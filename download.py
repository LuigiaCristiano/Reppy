#!/usr/bin/env python
#####################################################################################################
#
#               PLEASE add comments on changes on top with date and your initials.
#               makes traceback of changes much easier.
#
#####################################################################################################



##############
### 130819 (CW): small bugfixes on data consistency checks and length trimming of horizontals before rotation
#
### 130816 (CW): added option '--forcedl' to enforce redownload, does not check whether files are already present
#                default port reset to 18002
#
### 130815 (CW): 
###              VERSION working for ObsPy 0.8.4 ####
#                modified routeflag settings for italian data, only routeflag needs to be set in each loop, port can be same
#                italian data was available on port 18001, while 18002 is default ...
#
#                changed downsample function in d2fcts (using decimate) with new resample method
#                decimate introduces phase shift, resample method does not
#
####130430 :if condition on network name to download italian data using special port and route flag
### 130322 (CW): patch of temporary folder where obspy.arclink creates tmp files while downloading
#                previously tmp files were in /var/tmp but not removed automatically
#                in /tmp they will be automatically deleted once the stream object is read.
#
#                minor rearrangements and commenting
#
### 120705 (CW): included loadmat2trace function, improved reading of mat-files.
#                IMPORTANT: bug in save_stream_object for .mat files resolved
#                           only one component (commonly T, if rotation applied) was saved in three files
#
### 120704 (CW): stationname in input file may contain locationcode, e.g. RAYN.10, to pre-select specific instrument
#                does not apply to already downloaded data, which is loaded and converted.
#                some minor changes:
#                  - stream length is checked against (eventbased) requests, 
#                    if trace starts 5% (of duration) late or ends 20% early from request it is dropped
#                  - clients are opened before main loop
### 120529 (CW): added syncpath: folder to primary database e.g. /rayleighdata/data where eventfolders reside (2010/201012312359)
### 120524 (CW): comment lines in input file ('#') to skip them
### 120523 (CW): check start and endtimes of downloaded trace against requested times
#
### 120515 (CW): major rework:
#               -- introduced command line options:
#                        -r, -p, -s, -o (see USAGE)
#                       if given, they override input file
#               -- merged functions in single module d2fcts
#               -- output channels are determined by component field in input file,
#                       can be any combination of 'ZNERT'. '*' defaults to 'ZRT' if event info is given, else 'ZNE'
#                       only those channels will be saved
#               -- channels xH1 xH2 from IRIS are rotated to return xHN xHE (IRIS only! see retrieveWF.py for details)
#               -- stationfile (gmtlatlon....) is updated if a requested station is not present in the list
#
###############change 04092012 former  line 219
#		the  line was creating error message  #if Stationr != Station: del Location 
#		 this has been temporary changed with  if Stationr != Station: Location=[] to avoid the error message when "*" is used for the entry station
######################
#                 To do: update stationinfofile if we cut e and n components
####		change 04092012  previously the event information file was produced also when no data was available, now is produced only if the data are available.
##		succ_down and fail_download were created in the event folder. I moved in the parent folder. data folder was created even if no data is available. I create it only if data are available
##############################condition on  outpath. If is not given by command line it is the current folder
#################

### Check if command line is correct before loading all the modules
import sys
if len(sys.argv) <= 2:
        print ""
        print "USAGE: %s inputfile [-r -p -s sps -o path -syncpath --forcedl]" % sys.argv[0]
        print "options:"
        print "-r: raw data flag, save also data in raw format, no instrument correction"
        print "-p: print flag, print progress info on screen"
        print "-s sps: resample traces to 'sps' (samples per second), e.g. '-s 1' makes LH channels"
        print "-o path: override output data folder from input file to 'path'"
        print "-syncpath path: 'path' to primary database against which requests will be checked. default: output data folder"
        print "--forcedl: force download, does not check whether files are already present"
        print ""
        if len(sys.argv) < 2:
                sys.exit()
        else:
                print ""
                print "press any key to continue (x to exit)"
                dd = raw_input()
                if dd == 'x': sys.exit()
# default settings
import numpy as np
#import matplotlib.pyplot as plt
import csv
import time
import pickle, os, threading
import subprocess
import tempfile ## needed to set 'correct' path to /tmp directory where obspy.arclink creates tmp-files while downloading

# ObsPy modules
from obspy.arclink import Client
from obspy.iris import Client as iClient
from obspy.core import Stream, Trace, UTCDateTime,read
from obspy.signal.rotate import rotate_NE_RT
from obspy.core.util.geodetics import gps2DistAzimuth
from scipy.io import loadmat

# own modules, functions
from retrieveWF  import fetchwebdc, fetchiris
import d2fcts as d2f
from loadmat2trace import loadmat2trace


raw_data_flag = False # save raw data
print_flag = False  # print progress info, verbose mode
force_flag = False  # enforce download
outsps = None       # resample data
outpath = None                   #################where info files on download procedure are saved##############
syncpath = None                  #####################permanent database folder################

#OPEN INPUT FILE AND READ PARAMETERS
InputFile = str(sys.argv[1]); del sys.argv[1]

## loop over OPTIONAL COMMAND LINE ARGUMENTS
while len(sys.argv) > 1:
        option = sys.argv[1]
        del sys.argv[1]

        if option == '-r':
                raw_data_flag = True
#               raw_data_format = sys.argv[1]; del sys.argv[1]
        if option == '-p' or option == '--printonscreen':
                print_flag = True
        if option == '-s': # output samples per seconds (downsampling if necessary)
                outsps = int(sys.argv[1]); del sys.argv[1]
        if option == '-o': # outpath
                outpath = str(sys.argv[1]); del sys.argv[1]
        if option == '-syncpath': # path to primary database for check if request already exists
                syncpath = str(sys.argv[1]); del sys.argv[1]
                if not syncpath[-1] == '/': syncpath = syncpath + '/'
                if not os.path.isdir(syncpath):
                        print "SYNCPATH directory does not exist"
                        sys.exit()
        if option == '--forcedl': # force download
                force_flag = True

## check and set correct tmp-dir /tmp :: default /var/tmp will not remove created tmp-files!
if tempfile.tempdir != '/tmp':
        tempfile.tempdir = '/tmp'

# LOGFILE
### open Logfile, as 'InputFile basename'.log
Log_FileName=os.path.basename(os.path.splitext(InputFile)[0] + '.log')
Log_File=open(Log_FileName,'a')
## counter of lines in input file, i.e. no of requests
LineNo=0

### OPEN clients
host = "webdc.eu"
port = 18001 # 18001 is arclink proxy at GFZ (email Andres Heinloo, GFZ, 130816)
clientwebdc = Client(host, port, user="CAU Kiel")
clientiris = iClient(user="CAU Kiel")
#routeflag='True'
### OPEN InputFile
data_down_infos=csv.reader(open(InputFile,"r"), skipinitialspace=True) # allows for blanks in input line
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### LOOP 1: input file line by line
for line in data_down_infos: # line is list of elements in input line
        LineNo+=1
        if print_flag:
                print "+++   +++   +++   +++   +++   +++   +++   +++   +++"
                print 'Line No',LineNo
                print'time: ',time.strftime("%d.%m.%Y:%H:%M:%S",time.localtime())
        if line[0][0] == '#':
                if print_flag: print "skip:", ", ".join(line)
                continue
        else:
                if print_flag: print "line:", ", ".join(line)

        [torigin, event_lat, event_lon, depth, mag, path_data, dist_class, duration_class, toll_class, \
         rotation, correction, Network, Station, Location, Channel, sub, dataformat] = d2f.read_input_file(line)
        ### some 'dead' variables: rotation, 
        ### important variables:
                #toll: window should begin 'toll' seconds BEFORE origin time
###################italian data is available only through port 18001 and routeflag==False###########
        if Network in ('IV', 'MN', 'GU'):
                 routeflag=False
        else:
                routeflag=True
        Channel.sort(); Channel.reverse()
        reqcomp = "".join([i[-1] for i in Channel]) # requested components, e.g. 'ZTR', 'ZNE'


        Log_File.write('Data Downloading Start Time of This Event file\t')
        Log_File.write(time.strftime("%d.%m.%Y:%H:%M:%S",time.localtime())) 
        Log_File.write('\n') 

# determine (and create) data output folder
        # override path_data from input file by command line if given
        if outpath: path_data = outpath
        # check if path_data is empty -- CRITICAL, may happen if neither input file nor command line contains output folder
        if not path_data:
                path_data = "." # save in current folder
                #print "No output data folder given ... CRASH!"
                #sys.exit()
        if not outpath:
		outpath="."   
		###outpath=path_data
        # path_data must have trailing '/'
        if not path_data[-1] == '/': path_data = path_data + '/'

        # append event info to path_data
        if sub: # write data to subfolders of YYYY/YYYYMMDDHHMM
                path_data=str(path_data) + torigin.strftime("%Y/%Y%m%d%H%M/")
                #print path_data
                syncdir = path_data
        if syncpath:
                syncdir = syncpath + torigin.strftime("%Y/%Y%m%d%H%M/")

        if print_flag:
                print "syncfolder %s " % syncdir
                print "data output folder %s" % path_data

        if not os.path.exists(path_data):
                os.umask(0)
                os.makedirs(path_data, mode=0775)
                if print_flag:
                        print " ... newly created"

        print outpath
	succ_down=file("%ssucc_download.txt" %(outpath),'a')
        fail_down=file("%sfail_download.txt" %(outpath),'a')
#        succ_down=file("%ssucc_download.txt" %(path_data),'a')
#        fail_down=file("%sfail_download.txt" %(path_data),'a')


# and select proper client
        flag_IRIS = d2f.check_IRIS_net_list(Network)
        if flag_IRIS==0: 
                client = clientwebdc
        elif flag_IRIS==1:
                client = clientiris
        else:
                print "ERROR on flag_IRIS"
                sys.exit()


# select stations and their coordinates
        try:
                [Stationr,stalat,stalon] = d2f.rsi(Network, Station, flag_IRIS, torigin, max(duration_class),Channel,print_flag,routeflag,port) ## NEW, Stationr, stalat, stalon are lists!
        except:
                if 'No response' in str(sys.exc_info()[1]):
                        msg = 'StationNotFound'
                else:
                        msg = str(sys.exc_info()[1]).replace(" ", "")
                # keep line from input file
                fail_down.write("%s %s \n" % (", ".join(line), msg))
                Log_File.write("%s\t%s\t%s\t%s\t%d  failed station info not found\n" % (torigin.strftime("%Y\t%m\t%d\t%H\t%M\t%S"), Network, ",".join(Station),"_".join(Channel),LineNo))

                if print_flag:
                        print 'Station %s not found. msg = %s' % (Station, msg)

                continue


        #if Stationr != Station: del Location # remove preselected Locationcodes, if list(Station) included wildcards, i.e. was changed
        if Stationr != Station: Location=[]
        if print_flag:
                if len(Stationr) > 1: print "list of stations %s" % (", ".join(Stationr))

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### LOOP 2: for all stations in request line
        for inds in range(len(Stationr)): 
                t0 = time.time()
                if event_lat: # if event coords are given, determine distance class and window length (duration + pre-origin)
                        dist,az,baz=gps2DistAzimuth(event_lat,event_lon,stalat[inds],stalon[inds])
                        for ind_v in range(len(dist_class)):
                                if (dist/111000)<int(dist_class[ind_v]):
                                        laufzeit = duration_class[ind_v]
                                        toll_time = toll_class[ind_v]
                       #d2f.save_location_travel_time(path_data, Stationr[inds], stalat[inds], stalon[inds], event_lat, event_lon, depth, mag, dist, az, baz, torigin, toll_time, Channel)

                else:
                        laufzeit = max(duration_class)
                        toll_time = 0 # min(toll_class)
                        # remove eventual RT components from req. channels/components (since no event info given)
                        if 'R' in reqcomp: del Channel[[i[-1] for i in Channel].index('R')]
                        if 'T' in reqcomp: del Channel[[i[-1] for i in Channel].index('T')]
                        reqcomp = "".join([i[-1] for i in Channel]) # requested components, e.g. 'ZNE'

                ### start time of request (incl. pre-origin toll_time)
                tstart = torigin - toll_time
                tend = tstart + laufzeit

###############################################################################
                ### check if files to be downloaded are already available
                # - available: full path of available file (no further tasks)
                # - wrongfmt: full path of available files in wrong dataformat but correct channel (conversion needed)
                # - rotateme: full path of available NE files to be rotated
                reqchan = Channel[:]
                if outsps == 1: # check if 'LHx' channels are present
                        reqchan = ['L' + chn[1:] for chn in Channel]

                if not force_flag: # skip if download forced by option --forcedl
                        available, wrongfmt, rotateme = d2f.check_data_avail(tstart, Network, Stationr[inds], reqchan, laufzeit, correction, syncdir, dataformat)
                else:
                        available = []
                        wrongfmt = []
                        rotateme = []

                if print_flag:
                        if wrongfmt: 
                                print "wrongfmt files found in %s: " % syncdir
                                print "%s" % (", ".join([i.split("/")[-1] for i in wrongfmt]))
                        if available: 
                                print "available files found in %s: " % syncdir
                                print "%s" % (", ".join([i.split("/")[-1] for i in available]))

                # convert wrong format, will also downsample existing files, and add to available
                for infile in wrongfmt: 
                        infmt = infile.split(".")[-1]
                        if infmt in ['mseed', 'sac', 'gse']:
                                st = read(infile)
                        elif infmt == 'mat':
                                tr = loadmat2trace(infile)
                                st = Stream([tr])
                        else:
                                if print_flag: print "unrecognised fileformat ... skipping conversion of %s" % infile
                                continue

                        # decimate, if necessary
                        if outsps:
                                d2f.downsample(st, outsps)

                        print "........ converting %s ...." % infile
                        print st
                        print st[0].stats
                        print "........................."
                        d2f.save_stream_object(st, reqcomp, dataformat, path_data, correction, available)
                        del st

                # now we have available and rotateme to deal with

                # if we have for each channel one file available, all is done
                if len(Channel) == len(available): 
                        if print_flag: print "request available or solved by format conversion / downsampling" # in %s" % syncdir
                        Log_File.write("%s\t%s\t%s\t%s\t%d  converted existing file in %s\n" % (torigin.strftime("%Y\t%m\t%d\t%H\t%M\t%S"), Network, Stationr[inds],"_".join(Channel),LineNo, syncdir))
                        continue 

                if print_flag:
                        if rotateme: print "rotateme files: %s" % (", ".join([i.split("/")[-1] for i in rotateme]))

                # do we need to download additional files?
                dflag = True
                if available or rotateme:
                        if not available:
                                av_chan = rotateme[:] # copy instead of reference
                        if not rotateme:
                                av_chan = available[:] # copy instead of reference
                        else:
                                av_chan = available[:]
                                [av_chan.append(i) for i in rotateme]

                        av_comp = [i.split(".")[2].split("_")[0][-1] for i in av_chan] # available components incl. rotateme, e.g. ['Z', 'N']

                        Channel_copy = Channel[:] # test which channels are available/rotateable
                        for chn in Channel:
                                if not Channel_copy: break # all needed channels are available
                                if chn[-1] in av_comp: # channel is already available, done
                                        del Channel_copy[Channel_copy.index(chn)]
                                if (chn[-1] == 'R' or chn[-1] == 'T') and ('N' in av_comp and 'E' in av_comp): # rotation can be done
                                        del Channel_copy[Channel_copy.index(chn)] # remove R or T

                        if not Channel_copy: dflag = False  # all needed files are present, no download needed

                # initialize Stream objects
                st = Stream()
                straw = Stream()

                # download all channels with wildcard, if not all required channels are available
                if dflag:
                        t_down = time.time()
                        if print_flag:
                                print "downloading all channels for station: %s" % Stationr[inds]
                        dlchn = Channel[0][:2] + '*'
                        if Location and Location[inds] != "":
                                loclist = [Location[inds]]
                        else:
                                loclist = ['', '00', '10', '01']
                        for locationcode in loclist: # infinite loop: caution needs break at end otherwise will hang
                                try:
                                        if flag_IRIS==0:
                                                sttemp,strawtemp =fetchwebdc(client, Network, Stationr[inds], dlchn, locationcode, tstart, laufzeit,correction,True,routeflag)
                                        else:
                                                sttemp,strawtemp =fetchiris(client, Network, Stationr[inds], dlchn, locationcode, tstart, laufzeit,correction,saveraw=True)
                                        straw +=strawtemp
                                        st += sttemp

					## #########################SAVE INFOS ON THE EVENT LOCATION AND CALCULATE TRAVEL TIMES#######
                                        if sttemp and event_lat: # event location can only be saved if it is given, for example 'dist' is not defined if there is no event
                                                d2f.save_location_travel_time(path_data, Stationr[inds], stalat[inds], stalon[inds], event_lat, event_lon, depth, mag, dist, az, baz, torigin, toll_time, Channel)
	
                                        if print_flag:
                                                print "... successful %s, %s, \'%s\'" % (sttemp[0].stats.network, sttemp[0].stats.station, sttemp[0].stats.location)
                                        succ_down.write("%s\t%s\t%s\t%s\n" % ((torigin.strftime("%Y\t%m\t%d\t%H\t%M\t%S"), sttemp[0].stats.network, sttemp[0].stats.station, sttemp[0].stats.channel)))
                                        Log_File.write("%s\t%s\t%s\t%s\t%d  Downloaded Now \n" % (torigin.strftime("%Y\t%m\t%d\t%H\t%M\t%S"), sttemp[0].stats.network, sttemp[0].stats.station, "_".join([tr.stats.channel for tr in sttemp]), LineNo))
                                        break
                                except: # if no data, try different locationcode, else (e.g. timeout) stop
                                        if 'No waveform data' in str(sys.exc_info()[1]) or 'No data' in str(sys.exc_info()[1]): # IRIS, WebDC no data
                                                msg = 'NoData'
                                                if print_flag: print "No data for location \'%s\'" % locationcode
                                                if locationcode != loclist[-1]: # try with next location code
                                                        continue
                                        elif 'timeout' in str(sys.exc_info()[1]).lower(): # IRIS, WebDC no data
                                                msg = 'Timeout'
                                        else:
                                                msg = str(sys.exc_info()[1]).replace(" ", "")
                                        if print_flag: print "DL failed: %s" % msg #str(sys.exc_info()[1])
                                        line[15] = Stationr[inds]
                                        fail_down.write("%s %s \n" % (", ".join(line), msg))
                                        Log_File.write("%s\t%s\t%s\t%s\t%d  failed download\n" % (torigin.strftime("%Y\t%m\t%d\t%H\t%M\t%S"), Network, Stationr[inds],"_".join(Channel),LineNo))


## TODO old, keep until other than no data error is encountered and above procedure confirmed
#                                       if print_flag:
#                                               if 'No waveform data' in str(sys.exc_info()[1]) or 'No data' in str(sys.exc_info()[1]): # IRIS, WebDC no data
#                                                       print "No data for location \'%s\'" % locationcode
#                                               else:
#                                                       print "DL failed: %s" % str(sys.exc_info()[1])
#                                       if locationcode == loclist[-1]: # only when last locationcode fails
#                                               # catch exception
#                                               if 'No waveform data' in str(sys.exc_info()[1]) or 'No data' in str(sys.exc_info()[1]): # IRIS, WebDC no data
#                                                       msg = 'NoData'
#                                               elif 'timeout' in str(sys.exc_info()[1]).lower(): # IRIS, WebDC no data
#                                                       msg = 'Timeout'
#                                               else:
#                                                       msg = str(sys.exc_info()[1]).replace(" ", "")
#                                               # keep line from input file, replace station (may be wildcard) by current Stationr[inds]
#                                               line[15] = Stationr[inds]
#                                               fail_down.write("%s %s \n" % (", ".join(line), msg))
#                                               Log_File.write("%s\t%s\t%s\t%s\t%d  failed download\n" % (torigin.strftime("%Y\t%m\t%d\t%H\t%M\t%S"), Network, Stationr[inds],"_".join(Channel),LineNo))

                        if print_flag: print "download time: ", time.time() - t_down
	 
                # load rotateme - only if this avoids an additional download (i.e. Z is available, NE is rotateable), 
                #     else downloading 3 components is more consistent
                else: # RAW will be empty and not saved
                        if print_flag and rotateme:
                                print "loading files for rotation: %s" % (", ".join(rotateme))
                        for infile in rotateme:
                                infmt = infile.split(".")[-1]
                                if infmt in ['mseed', 'sac', 'gse']:
                                        st += read(infile)
                                elif infmt == 'mat':
                                        tr = loadmat2trace(infile)
                                        st += Stream([tr])


  
                ### should not happen, but to be sure
                if not st:
			if print_flag:
                        	print "No data found - either from download or files to be rotated .... skipping"
                        continue

###############################ROTATION
                ## do rotation and write data to file (should apply only to event data, not continuous "noise" recordings)
                # check start and endtimes of received data
                # starttime may be 5% of requested duration late from ORIGINtime(!)
                # endtime may be 20% of requested duration early
                # ... other traces are removed
                if event_lat: # and ('R' in reqcomp or 'T' in reqcomp): 
                        for tr in st: 
                                if tr.stats.starttime - torigin > laufzeit * .05: # earlier starttimes are acceptable
                                        #print tr.stats.starttime, torigin
                                        if print_flag: print "starttime of channel %s differs from request by %ss" % (tr.stats.channel, tr.stats.starttime - torigin)
                                        Log_File.write("%s\t%s\t%s\t%s\t%d  failed starttime too late\n" % (torigin.strftime("%Y\t%m\t%d\t%H\t%M\t%S"), tr.stats.network, tr.stats.station,tr.stats.channel,LineNo))
                                        print st
					st.remove(tr)
                                        continue
                                if tend - tr.stats.endtime > laufzeit * .2: # later endtimes are acceptable
                                        if print_flag: print "endtime of channel %s differs from request by %ss" % (tr.stats.channel, tend - tr.stats.endtime)
                                        Log_File.write("%s\t%s\t%s\t%s\t%d  failed endtime too early\n" % (torigin.strftime("%Y\t%m\t%d\t%H\t%M\t%S"), tr.stats.network, tr.stats.station,tr.stats.channel,LineNo))
                                        st.remove(tr)
                                        continue

                if len(st) == 0: # all traces had been removed
                        if print_flag: print " ... no traces left ... nothing left to do ... skipping ..."
                        Log_File.write("... no traces left ... skipping")
                        continue

                prescomp = [p.stats.channel[-1] for p in st] # available components
                if event_lat and ('N' in prescomp and 'E' in prescomp) and ('R' in reqcomp or 'T' in reqcomp): # then there should be as well nraw and eraw
                        for tr in st:
                                if tr.stats.channel[-1] == 'N': n = tr #.copy()
                                if tr.stats.channel[-1] == 'E': e = tr #.copy()
                        for tr in straw: # may be empty, if only rotateable files are read -> in the following protect straw with "if"
                                if tr.stats.channel[-1] == 'N': nraw = tr #.copy()
                                if tr.stats.channel[-1] == 'E': eraw = tr #.copy()

                        # cut the horizontal components if they havent  same length
                        if len(e)!=len(n): # same applies to nraw and eraw
                                starttrim=max(e.stats.starttime, n.stats.starttime)
                                endtrim=min(e.stats.endtime, n.stats.endtime)
                                if print_flag: 
                                        print "trim EN to equal length..."
                                        #print "Trimtimes start-end:", starttrim, endtrim
                                if endtrim > starttrim:
                                        e.trim(starttrim,endtrim)
                                        n.trim(starttrim,endtrim)
                                        if straw:
                                                eraw.trim(starttrim,endtrim)
                                                nraw.trim(starttrim,endtrim)


                                # cure remaining difference (there may still be a difference of 1 sample)
                                diff = e.stats.npts - n.stats.npts # diff should be zero, then nothing is done
                                if diff > 0: # crop samples at the end of E
                                        if print_flag: print "crop E by an extra %i samples at the end" % diff
                                        e.data = e.data[:-diff]
                                        if straw: eraw.data = eraw.data[:-diff]
                                if diff < 0: # crop samples at the end of N
                                        if print_flag: print "crop N by an extra %i samples at the end" % -diff
                                        n.data = n.data[:diff]
                                        if straw: nraw.data = nraw.data[:diff]

                        #print "Length E, N: ", len(e.data), len(n.data)
                        #print "compare lengths: ", len(e.data)==len(n.data)
                        #print "Starttime E, N: ", e.stats.starttime, n.stats.starttime
                        #print "Endtime E, N: ", e.stats.endtime, n.stats.endtime
                        r,trasv=rotate_NE_RT(n.data, e.data, baz)
                        #print "Length RAW E, N: ", len(eraw.data), len(nraw.data)
                        #print "compare RAW lengths: ", len(eraw.data)==len(nraw.data)
                        if straw: rraw,trasvraw=rotate_NE_RT(nraw.data, eraw.data, baz)

                        stats_r = tr.stats.copy()
                        stats_r.channel = stats_r.channel[:-1]+"R"
                        st += Stream([Trace(data=r, header=stats_r)])
                        if straw: straw += Stream([Trace(data=rraw, header=stats_r)])

                        stats_trasv = tr.stats.copy()
                        stats_trasv.channel = stats_trasv.channel[:-1]+"T"
                        st += Stream([Trace(data=trasv, header=stats_trasv)])
                        if straw: straw += Stream([Trace(data=trasvraw, header=stats_trasv)])

                # downsample if outsps is given
                if outsps:
                        d2f.downsample(st, outsps)

###########################################SAVE STREAM OBJECT#####################
                if print_flag:
                        print "writing Stream object:"
                        for tr in st:
                                if tr.stats.channel[-1] in reqcomp:
                                        print tr
                        if raw_data_flag:
                                print "... and RAW Stream object:"
                                for tr in straw:
                                        if tr.stats.channel[-1] in reqcomp:
                                                print tr
                if not os.path.exists(path_data):
                        os.umask(0)
                        os.makedirs(path_data, mode=0775)
                        if print_flag:
                                print " ... newly created"
 
                d2f.save_stream_object(st, reqcomp, dataformat, path_data, correction)
                if raw_data_flag:
                        d2f.save_stream_object(straw, reqcomp, dataformat, path_data, 'RAW')
###################################################################
                if print_flag: print "overall time for station %s: %i seconds" % (Stationr[inds], time.time() - t0)
###### LOOP 2 END

###### LOOP 1 END
fail_down.close()
succ_down.close()
Log_File.close()
