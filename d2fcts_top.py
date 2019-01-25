#!/bin/env python

import csv,sys,os
import subprocess
import glob
import commands
import os.path
import numpy as np

from scipy.io import savemat
from obspy.core import UTCDateTime, Stream, Trace, read #, UTCDateTime # potentially useful in rsi
from arclink import Client
from obspy.clients.fdsn import Client as iClient
from obspy.taup import TauPyModel
from obspy.signal.util import next_pow_2 as nextpow2
###################################SAVE STREAM OBJECTS THE ENCODING IS RELTED TO THE TYPE OF DATA WE SAVE, WE CANNOT SAVE DATA AFTER INSTRUMENT CORRECTION IN ANY OTHER FORMAT THAN ENC=5#############
def save_stream_object(st, reqcomp, dataformat, path, wtype, available=None):
        '''
        st: Stream object with traces to be saved
        reqcomp: string of requested components, e.g. 'ZTR'
        dataformat: string of requested data formats, e.g. 'mSM'
        path: string, path where data should be saved (includes trailing '/')
        wtype: string, type of seismogram, e.g. 'RAW', 'VEL'

        available (optional): list of filenames which is extended and returned by written files
        '''
        enc = 5 # float64
        if wtype == "RAW":
                enc = 'STEIM1' # Steim 1

        for tr in st:
                if not tr.stats.channel[-1] in reqcomp: continue # skip unrequested channels
                tracelen = int(tr.stats.endtime - tr.stats.starttime)
                if 'M' in dataformat: #mseed
                        filename = "%s%s.%s.%s_%s_%s_%i.mseed" % (path, tr.stats.network, tr.stats.station, tr.stats.channel, wtype, tr.stats.starttime.strftime("%Y-%2m-%2d.%2H-%2M-%2S"), tracelen)
                        tr.write(filename, format="MSEED", reclen=512, encoding=enc)
                if 's' in dataformat: #sac
                        filename = "%s%s.%s.%s_%s_%s_%i.sac" % (path, tr.stats.network, tr.stats.station, tr.stats.channel, wtype, tr.stats.starttime.strftime("%Y-%2m-%2d.%2H-%2M-%2S"), tracelen)
			if wtype=="RAW":
				if tr.stats.channel[-1]=="Z" or tr.stats.channel[-1]=="E" or tr.stats.channel[-1]=="N":
                                   tr.write(filename, format="SAC", reclen=512, encoding=enc)
			else:
                               tr.write(filename, format="SAC", reclen=512, encoding=enc)

                if 'g' in dataformat: #gse
                        filename = "%s%s.%s.%s_%s_%s_%i.gse" % (path, tr.stats.network, tr.stats.station, tr.stats.channel, wtype, tr.stats.starttime.strftime("%Y-%2m-%2d.%2H-%2M-%2S"), tracelen)
                        tr.write(filename, format="GSE2", reclen=512, encoding=enc)
                if 'm' in dataformat: #mat
                        filename = "%s%s.%s.%s_%s_%s_%i.mat" % (path, tr.stats.network, tr.stats.station, tr.stats.channel, wtype, tr.stats.starttime.strftime("%Y-%2m-%2d.%2H-%2M-%2S"), tracelen)
                        #for  ii, tr in enumerate(st):
                        mdict = dict([[j, str(k)] for j,k in tr.stats.iteritems()])
                        mdict['data'] = tr.data
                        #print "--> saving mat file: %s" % filename
                        #print mdict
                        savemat(filename,  mdict)
                # change permissions
                subprocess.call(["chmod", "664", filename])

                if isinstance(available, list):
                        available.append(filename)

        if isinstance(available, list):
                return available
        else:
                return

#======================================================================================================================
#======================================================================================================================

#########################################checks if the requested data are already present in the data base WE DISTINGUISH BETWEEN REQUEST WITH ROTATION AND WITH No ROTATION AS WELL BETWEEN COMPONENT="*" AND A LIST OF Components#####
def check_data_avail(treq, Network, Station, Channel, duration, correction, path, dataformat):
        '''
        treq: UTCDateTime object, requested starttime
        Network: str
        Station: str
        Channel: str, list of requested channels, e.g. ['BHZ', 'BHR', 'BHT']
        duration: int, length of requested tracce
        correction: str, output unit of request ('DIS', 'VEL', 'ACC')
        path: str, data path
        dataformat: str, one or more of 'msgM' (mat, sac, gse, MSeed)
        '''

        ### requested duration
        treqend = treq + duration # UTCDateTime object

#############################which format should be searched############
        reqfmt = []
        if 'M' in dataformat:
                reqfmt.append("mseed")
        if 'g' in dataformat:
                reqfmt.append("gse")
        if 's' in dataformat:
                reqfmt.append("sac")
        if 'm' in dataformat:
                reqfmt.append("mat")

##   sort Channel in descending order [ZTRNE]
        Channel.sort()
        Channel.reverse()

###### NEW SEARCH: we know the folder, find all components for event, station -> then check channels
        nome_trace_search = path + "*" + str(Network) + "." + str(Station) + '*' # select only network, station
        filefound = {'filename': [], 'fmt': [], 'comp': []}
        available = []
        wrongfmt = []
        rotateme = []
        for fullpath in glob.glob(nome_trace_search): ## e.g. full path /rayleighdata/data/2010/201010132243/GR.WET.BHZ_VEL_2010-10-13.22-39-03_2000.mat
                filename = fullpath.split('/')[-1]    ## e.g. filename GR.WET.BHZ_VEL_2010-10-13.22-39-03_2000.mat 
        	splt1, splt2,splt3,splt4 = filename.split("_") ### e.g. ['GR.WET.BHZ', 'VEL', '2010-10-13.22-39-03', '2000.mat']
                dur_found = int(splt4.split(".")[0])
                tfound = UTCDateTime(splt3.replace("-", ":").replace(".","T")) # starttime of found trace
                tfoundend = tfound + dur_found # UTCDateTime object
                # compare start and endtimes of request and found trace, times are within 10s corridor?
                if not (abs(treq - tfound) <= 10 and abs(treqend - tfoundend) <=10):
                        #print "file %s, times are more than 10s off ..." % filename
                        continue
                ### check output unit, continue also if raw seismogram is available
                if not (splt2 == correction or splt2 == 'RAW'): continue

                filefound['filename'].append(fullpath)

        ## sort found files in descending order [ZTRNE]
        filefound['filename'].sort()
        filefound['filename'].reverse()
        #print filefound, len(filefound['filename'])
        for f in filefound['filename']:
                filefound['fmt'].append(f.split(".")[-1])
                filefound['comp'].append(f.split(".")[-3][:3])

        ## possibly matching files are listed in filefound, including wrong formats and components
        for chn in Channel: # each requested channel
           for fmt in reqfmt:
                for i in range(len(filefound['comp'])):
                # 1. search for full match, i.e. component and fileformat matches
                        if filefound['comp'][i] == chn and filefound['fmt'][i] == fmt:
                                available.append(filefound['filename'][i])
                                break
                else: # no full match
                        for i in range(len(filefound['comp'])):
                                if filefound['comp'][i][0] == 'L' and filefound['comp'][i][0] != chn[0]: # LH is available but higher sampling rate (BH, HH, or other) is requested
                                        continue
                                elif filefound['comp'][i][-1] == chn[-1]: # right component, wrong format ## or just oversampled
                                        wrongfmt.append(filefound['filename'][i])
                                        break
                                # TODO: finds also BHN BHE if LHR LHT are requested
                                elif chn[-1] in 'RT' and filefound['comp'][i][-1] in 'NE': # wrong but useful channel (NE found for RT), format does not matter, file needs to be read anyway
                                        if not filefound['filename'][i] in rotateme: rotateme.append(filefound['filename'][i]) # avoid duplicate entries


        return available, wrongfmt, rotateme


#======================================================================================================================
#======================================================================================================================


def read_input_file(data):
        """
        read a text file containing on each line a  "request" of data download 
        """

# 120704 (CW): Station name may include locationcode, e.g. UOSS.00,
#              which is separated here in list(Station) and list(Location)

        #### origin time of event
        y=int(data[0])  #year
        mon=int(data[1]) #month
        d=int(data[2]) #day
        h=int(data[3]) #hour
        m=int(data[4]) #minute
        s=int(data[5])
        ms=0 #int(data[6])
        t = UTCDateTime(y, mon,d, h,m,s,ms/1000)
        ##### event parameters  # optional, may be empty
        event_lat=(data[7])
        if event_lat: event_lat=float(data[7])
        event_lon=(data[8])
        if event_lon: event_lon=float(data[8])
        depth=(data[10])
        if depth: depth=float(data[10])
        mag=(data[9])
        if mag: mag=float(data[9])

        dist_class=map(int, data[11].strip().split(';')) ### distance classes as int
        duration=map(int, data[12].strip().split(';')) ### request duration classes as int
        toll=map(int, data[13].strip().split(';')) ### preorigin time to be included (classes) as int

        ##### requested network, stations, channels
        Network=data[14].strip().upper()
        Station=data[15].strip().upper().split(';')
        ## NEW 120704: separate eventual locationcode from station, e.g. UOSS.00
        Location = list(Station) # copy of Station, overwrite below
        for i, sta in enumerate(Station):
                if "." in sta:
                        Station[i] = sta.split(".")[0].strip()
                        Location[i] = sta.split(".")[1].strip()
                else:
                        Station[i] = sta.strip()
                        Location[i] = ""

        channel=data[16].strip().upper() # e.g. BH or LH
        component = data[17].strip().upper() # should be any combination of 'ZNERT' or '*'
        component = component.translate(None, ';') # remove accidental ';', if e.g. Z;N;E is given (CW 131120)
        if component == '*': 
                if event_lat: 
                        component = 'ZRT'
                else:
                        component = 'ZNE'
        # Channel is list of channel + component, sorted alphabetically descending e.g. list ['BHZ', 'BHE']
        Channel = [channel + x for x in component]
        Channel.sort(); Channel.reverse() # inverse sorting

        ###### output data format 
        dataformat=data[18].strip()


        #### additional fields are optional but order must be kept
        ###### rotation
        if len(data) > 19: rotation=data[19].strip() 
        else: rotation = True
        ###### instrument correction to 'DIS', 'VEL' or 'ACC'
        if len(data) > 20: correction=data[20].strip().upper() 
        else: correction = 'VEL'
        ###### output path folder
        if len(data) > 21: root_path=data[21].strip() 
        else: root_path = None
        ###### creation of subfolders is default
        if len(data) > 22 and not data[22]: sub = False # if not present, disable the creation of subfolders
        else: sub = True

        return  t, event_lat, event_lon, depth, mag, root_path, dist_class, duration, toll, rotation, correction, Network, Station, Location, Channel, sub, dataformat

#======================================================================================================================
#======================================================================================================================

#def read_station_informations(Network, Station, Channel, flag_IRIS, t, duration):
def rsi(Network, Station, flag_IRIS, t, duration, Channel,print_flag=False,routeflag=True):
        '''
        for a given Network and list of Stations, e.g. [BFO, HLG] or ['*'], find station
        coordinates from static list, or download if not available
        '''

# 120515 (CW): Stations that are not listed in stationfile are added to stationfile!!
        Stationr=[] # output list of stationnames
        stalat=[] # output list of station latitudes
        stalon=[] # output list of station longitudes

        path = sys.argv[0].rpartition("/")[0]
        if not path: path = '.'
        #stationfile = path + "/gmtlatlon_2011_2011All.dat"
        stationfile = path + "/gmtlatlon_all.txt"
	#print stationfile
	os.path.isfile(stationfile)
	if os.path.isfile(stationfile)== True:
#### 131126 (CW): stationfile is read for each station request
         stafid = open(stationfile, 'ra')
         for ind in Station[:]: # loop over copy because elements of Station may be removed in the loop
                # create stationlist from stationfile
                stafid.seek(0) # reset stationfile to start
                # case wildcard in station
                if '*' in ind:
                        stationlist = [line.split() for line in stafid if line.split()[0] == Network]
                # case explicit station name (no wildcard)
                else: 
                        stationlist = [line.split() for line in stafid if line.split()[1] == ind and line.split()[0] == Network] 
                for entry in stationlist:
                        #stationstarttime = UTCDateTime(datetime.datetime.strptime(entry[2], '%d.%m.%Y')) ## potentially skip stations that started later than request date
                        stationnameline=str(entry[1])
                        stalatline=float(entry[5])
                        stalonline=float(entry[4])
                        if '*' in ind: # wildcard in stationname, collect all available stations that match
                                if '*' in ind and stationnameline[:ind.index('*')] == ind[:ind.index('*')]: # allows for wildcards a la 'C*' or 'CL*'
                                        Stationr.append(stationnameline)
                                        stalat.append(stalatline)
                                        stalon.append(stalonline)
                                        if ind in Station: Station.remove(ind) # at least one match was found for wildcard, remove from Station
                        else: # explicit station name ind
                                if stationnameline==ind and entry[0] == Network: # check also for correct network, else download
                                        stalat.append(stalatline)
                                        stalon.append(stalonline)
                                        Stationr.append(stationnameline)
                                        if ind in Station: Station.remove(ind) # remove station from request list

        # if list of return stations is allocated, i.e. stations were found in file, return
        # !!! file may not be complete, if wildcards are requested !!!
	if Stationr:
                if not Station: # all requests (min. 1 station per wildcard) fulfilled
                        return Stationr, stalat, stalon 
        #ELSE: #  continue and download station information from WebDC or IRIS also if in stationlist file the number of available stations is less than requested
	latlon_file=open(stationfile,'a')
	loclist = ['', '00', '10' ,'01','02']
        starttime = [] # needed for extending stationfile
        endtime = [] 
        flag = [] # WebDC or IRIS
        printline = []
        channelzero=Channel[1]
        channelreq=channelzero[:2]+"*"
        # download from WebDC
        if flag_IRIS == 0:
                client = Client(user=username,timeout=200) 
			##################change '*' with Station##########################
                try:
                        invent = client.getInventory(Network,'*','*', channel=str(channelreq), route=routeflag)
        		stationsgot = [value.code for key, value in invent.items() \
                    	  if key.startswith(Network + '.') and "code" in value and str(value.code) !=Network]
			stalatgot = [value.latitude for key, value in invent.items() \
                          if key.startswith(Network + '.') and "latitude" in value]
			stalongot = [value.longitude for key, value in invent.items() \
                          if key.startswith(Network + '.') and "longitude" in value]
			starttimegot=[value.start for key, value in invent.items() \
                          if key.startswith(Network + '.') and "start" in value]
 			endtimegot=[value.end for key, value in invent.items() \
                          if key.startswith(Network + '.') and "end" in value]
                except:
                        raise
		print Stationr
                # process stadict and assemble output lists
                d2="None      "
                d3="WebDC"
		#########change with one loop and check the presence with a grep on stationname##
		##time consuming send a request for each station###########
                if Station[0][-1] == '*': # station request with wildcard (network is missing in stationfile      
                        Stationr=stationsgot
			stalat=stalatgot
		        stalon=stalongot
                	for indstagot in range(len(stationsgot)):
	                      	
                                d1=starttimegot[indstagot].strftime("%d.%m.%Y")
				if endtimegot[indstagot]:
                                                d2=endtimegot[indstagot].strftime("%d.%m.%Y")
				
                                commandstation = ("awk \'($1==\"%s\"){print $2, $5, $6, $7}\' %s" % (Stationr[indstagot], stationfile))
                                stationfound = subprocess.check_output(commandstation, shell=True)
                                if not stationfound:

                                  printline.append([Network, Stationr[indstagot], d1, d2,  str("%.4f" %stalon[indstagot]),  str("%.4f" % stalat[indstagot]), d3])
                else: # stations are given explicitly
                        if print_flag==True:
                        for indstagot in range(len(stationsgot)):
                                if  stationsgot[indstagot] in Station:
                                        Stationr.append(stationsgot[indstagot])
                                        stalat.append(stalatgot[indstagot])
                                        stalon.append(stalongot[indstagot])
                                        d1=starttimegot[indstagot].strftime("%d.%m.%Y")
					if endtimegot[indstagot]:
						d2=endtimegot[indstagot].strftime("%d.%m.%Y")
					commandstation = ("awk \'($1==\"%s\"){print $2, $5, $6, $7}\' %s" % (stationsgot[indstagot],stationfile))
                	                stationfound = subprocess.check_output(commandstation, shell=True)
                        	        if not stationfound:

                                		        printline.append([Network, stationsgot[indstagot], d1, d2,  str("%.4f" % stalongot[indstagot]),  str("%.4f" % stalatgot[indstagot]), d3])

                if not printline:
                        raise Exception('No response information')


        # download from IRIS
        elif flag_IRIS == 1: 
                client = iClient("IRIS")
                try:
                        inv = client.get_stations(network=Network, station=",".join(Station), starttime=t, endtime=t+duration, level="station")
                except: 
                        raise

                d3=" IRIS"
                # write stations to inventory
                # select correct network (should be only one, but to be safe)
                for net in inv:
                        if net.code == Network: break
                for sta in net: 
                        # specific stations are requested (no wildcard); then skip all others
                        if Station[0][-1] != '*' and sta.code not in Station: continue 

                        Stationr.append(sta.code)
                        stalat.append(float(sta.latitude)) # sta.latitude is obspy.station.util.Latitude object
                        stalon.append(float(sta.longitude)) # sta.latitude is obspy.station.util.Latitude object

                        d1 = sta.start_date.strftime("%d.%m.%Y")
                        d2 = sta.end_date.strftime("%d.%m.%Y")
                        printline.append([Network, sta.code, d1, d2, str("%.4f" % stalon[-1]), str("%.4f" % stalat[-1]), d3])

        if printline:
                # open stationfile and append new stations
                stafid = open(stationfile, 'a')
		if print_flag==True:
                	print "adding lines to %s" % stationfile
                	print ". . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . "
                for i in range(len(printline)):
                        stafid.write("\t\t\t".join(printline[i]))
                        stafid.write("\n")
 			if print_flag==True:
                        	print "\t\t\t".join(printline[i])
                stafid.close()
		if print_flag==True:
                	print ". . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . "

  	return Stationr,stalat,stalon

#======================================================================================================================
#======================================================================================================================

def save_location_travel_time(path, stationname, slat, slon, event_lat, event_lon, depth, mag, dist, az, baz, torigin, toll_time, st):
##########################################CREATE TXT FILES CONTAINING EVENT INFOS AND EPICENTRAL DISTANCE FOR ALL THE STATIONS ASKED IN THE REQUEST, CALCULATES THE TRAVELTIMES AND SAVE THEM IN A TXT FILE####

        # origin time
        dateorigin = torigin.strftime("%Y%m%d%H%M")

        # start time of trace
        tstart = torigin - toll_time
        date = tstart.strftime("%Y%m%d%H%M")
        
        # event infos
	fo = file("%sevent_infos_%s.txt" %(path, dateorigin), 'w')
        fo.write("%3.4f\t%3.4f\t%2.1f\t%1.1f\t %s\t%s" % (event_lat, event_lon, depth, mag, dateorigin, torigin.strftime("%H %M %S %f")))
        fo.close()

        # station infos
        f1 = file("%sstation_infos_%s.txt" %(path, dateorigin), 'a')

        ### ONE channel should be sufficient?? Else loop over channels here
        #comp = Channel[0][-1]
        #f1.write("%s\t%2.2f\t%2.2f\t%d\t%2.2f\t%2.2f\t%s\t%s\t%d\n " % (stationname, slat, slon, dist, az, baz, comp, date, toll_time))
        ### ALL channels option
	
        for tr in st:
		starttime= tr.stats.starttime
                toll_time=torigin - starttime  ###to have it positive
                if tr.stats.channel[-1]=="Z" or tr.stats.channel[-1]=="E" or tr.stats.channel[-1]=="N": 
                 f1.write("%s\t%2.4f\t%2.4f\t%d\t%2.3f\t%2.3f\t%s\t%s\t%d\n " % (stationname, slat, slon, dist, az, baz, tr.stats.channel[-1], date, toll_time))
        #if depth<=1:
#		depth=1.3 
        model=TauPyModel("ak135")
        #if localvelmodfilename:
        #    model=TauPyModel("localfilename")
        tt=model.get_travel_times(source_depth_in_km=depth,distance_in_degree=dist/111000.,phase_list=["Pg","Pn","PmP","PgPg","PnPn","p","pP","P","pPn","pS","S","Sg","SgSg","SnSn","S","Sn","P^410P","P^660P","P^220P","PvmP","Pvmp","Pdiff","PKP","PKIKP","pS","sP","sS","PcP","PP","PKiKP","pPKiKP","sPKiKP","SS","ScS","sPn","PnS","SKS","PKKP","SKKS"])  ####phase,distance,time,purist_dist,ray_param,ray_param_index,name,purist_name,source_depth,takeoff_angle,incident_angle
                                             #tt = getTravelTimes(dist/111000, depth, model='ak135')
        command='grep -w \"'+ str(stationname)+'\"  '+str(path)+"station_infos_"+str(dateorigin)+'.txt | head -n 1'
        
       # print command
       
        # subprocess.check_output is available from Python2.7, older versions use commands.getoutput
        # TODO: not nice fix of version, is that needed here anyway? a is not used other than its presence checked.
	try:
                a = subprocess.check_output(command, shell=True)
        except:
                a = commands.getoutput(command)

        for linet in tt:
                label=linet.phase.name
                timet=linet.time
                if not a:
                        f1.write("%s\t%2.3f\n" %(label, timet))

        f1.close()
        return

#======================================================================================================================
#======================================================================================================================

def check_IRIS_net_list(Network):

#############CHECK IF A NETWORK  TO IRIS
	import csv
        path = sys.argv[0].rpartition("/")[0]
        if not path: path = '.'
	Inputfile=path + "/IRIS_Network_list"
	stationlist=csv.reader(open(Inputfile,"r"))
	IRIS_flag=0
	for stations in stationlist:
                if str(Network) in stations: 
                        IRIS_flag=1

    	return IRIS_flag

#======================================================================================================================
#======================================================================================================================

def downsample(st, outsps, debug=False):

        # change length of traces to nextpow2, else FFT in resample can be VERY slow
        for i,tr in enumerate(st):
                if tr.stats.sampling_rate < outsps: continue # do nothing if upsampling is requested

                end = tr.stats.endtime # keep endtime to trim trace at the end
		old_sampl=tr.stats.sampling_rate
                tr.data = np.append(tr.data, np.zeros(nextpow2(tr.stats.npts) - tr.stats.npts)) # append zeros to make 2^N length
                tr.resample(outsps)
                tr.trim(endtime=end)
                if round(tr.stats.sampling_rate, 2) == 20  or round(tr.stats.sampling_rate, 2) == 40: #rename channels to LHx
                        tr.stats.channel = 'B' + tr.stats.channel[1:]
                if round(tr.stats.sampling_rate, 2) == 1: # rename channels to LHx
                        tr.stats.channel = 'L' + tr.stats.channel[1:]


        return

