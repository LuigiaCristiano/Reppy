#!/usr/bin/env python

# make input file for download of continuous data

import sys, os
import datetime

from obspy.core import UTCDateTime
from obspy.arclink import Client

def string2UTCDateTime(otime):
        #### fix erroneous string interpretation in UTCDateTime objects
        if otime[-1] == ",": otime = otime[:-1] # remove trailing comma
        otime = otime.split(",") # list of values
        if len(otime) < 3: 
                print "time too short, needs year, mon, day"
                sys.exit()
        elif len(otime) == 3: # yyyy, mm, dd
                otime = datetime.datetime.strptime(" ".join(otime), "%Y %m %d")
        elif len(otime) == 4: # yyyy, mm, dd, hh
                otime = datetime.datetime.strptime(" ".join(otime), "%Y %m %d %H")
        elif len(otime) == 5: # yyyy, mm, dd, hh, mm
                otime = datetime.datetime.strptime(" ".join(otime), "%Y %m %d %H %M")
        elif len(otime) == 6: # yyyy, mm, dd, hh, mm, ss
                otime = datetime.datetime.strptime(" ".join(otime), "%Y %m %d %H %M %S")
        elif len(otime) == 7: # yyyy, mm, dd, hh, mm, ss, microsecond
                otime[6] = str("%06i" % (int(otime[6]) * 1000))
                otime = datetime.datetime.strptime(" ".join(otime), "%Y %m %d %H %M %S %f")
                print "!!! microseconds will be ignored in request line!!!"
        otime = UTCDateTime(otime)

        return otime

###### MAIN PROGRAM

call = sys.argv[:]

if len(sys.argv) <= 3:
        print 'prepare input file for download routine to download continuous data'
        print ""
        print "USAGE: %s -t time [-e endtime -chn channels]" % sys.argv[0]
        #print "file: name of file to be created"
        print "-t time: startday of request, comma-separated list of yyyy, mm, dd"
        print ""
        print "options:"
        print "-e endtime: endday of request, comma-separated list of yyyy, mm, dd" 
        print "            default: yesterday"
        print '-chn channels: type of required channels, DEFAULT: LH'
#       print "-e eventloc: event location, comma-separated list of lat, lon, depth, mag"
#       print "-c classes: distance classes for signal length and pre-origin time, comma-separated list dist, length, preorigin time"
#       print "              each class can be semicolon-separated list IN QUOTATION MARKS, e.g. '180;20;3, 7200;1200;600, 0;10;120'"
#       print "              DEFAULT: '180, 7200, 0'"
#       print "-r region: rectangular region, comma-separated list lat1, lon1, lat2, lon2 - NO NEGATIVE LONGITUDES"
#       print "-comp comp: components to be downloaded, DEFAULT '*'"
#       print "-fmt format: format of seismograms, any combination of msgM - M:mseed (DEFAULT), m:mat, s:sac, g:gse"
#       print "-n network: restrict station search to given network"
        print ""
        sys.exit()

#outfile = str(sys.argv[1]); del sys.argv[1];


# default variables:
chn = 'LH*' # channels for station search, e.g. BH
cls = '180, 86400, 0' # classes (dist, duration, preorigin time)
net = "*" # network for station search
reg = [-90, 0, 90, 360] # list of lat1, lon1, lat2, lon2 for regional station search
fmt = "m" # mat
comp = '*' # components to be requested
eloc = ",,," # event location
stationlist = list() # list of stations to be downloaded
# 
otime = "" # starttime
etime = UTCDateTime(UTCDateTime().date) - 1 # default endtime, yesterday 23:59:59


while len(sys.argv) > 1:
        option = sys.argv[1]
        del sys.argv[1]

        if option == '-t':
                while len(sys.argv) > 1 and sys.argv[1][0] != '-':
                        otime += str(sys.argv[1])
                        del sys.argv[1]
                otime = string2UTCDateTime(otime)
                otimestr = otime.strftime("%Y, %m, %d, %H, %M, %S, 0") # string, e.g. '2008, 12, 31, 00, 00, 00, 0', ignore microseconds

        if option == '-e':
                etime = ""
                while len(sys.argv) > 1 and sys.argv[1][0] != '-':
                        etime += str(sys.argv[1])
                        del sys.argv[1]
                etime = string2UTCDateTime(etime)
                etime = UTCDateTime(etime.year, etime.month, etime.day, 23, 59, 59) # adjust to end of selected day
                etimestr = etime.strftime("%Y, %m, %d, %H, %M, %S, 0") # string, e.g. '2008, 12, 31, 00, 00, 00, 0', ignore microseconds

        if option == '-chn':
                chn = str(sys.argv[1]) # TODO extend to allow for lists??
                del sys.argv[1]
                if chn[-1] != '*': chn = chn + '*' # append '*'

#       if option == '-net':
#               net = str(sys.argv[1]) # TODO extend to allow for lists??
#               del sys.argv[1]

#       if option == '-fmt':
#               fmt = str(sys.argv[1]) # TODO extend to allow for lists??
#               del sys.argv[1]

#       if option == '-comp':
#               comp = str(sys.argv[1]) # TODO extend to allow for lists??
#               del sys.argv[1]

#       if option == '-e':
#               eloc = ""
#               while len(sys.argv) > 1 and sys.argv[1][0] != '-':
#                       eloc += str(sys.argv[1])
#                       del sys.argv[1]
#               if eloc[-1] == ",": eloc = eloc[:-1] # remove trailing comma
#               eloc = eloc.replace(" ", "") #nice formatting
#               #if eloc[-1] != ",": eloc += "," # add trailing comma 
#               eloc = eloc.replace(",", ", ") # string, e.g. "44.814, 11.079, 9.6, 5.8,"

#       if option == '-c':
#               cls = ""
#               while len(sys.argv) > 1 and sys.argv[1][0] != '-':
#                       cls += str(sys.argv[1])
#                       del sys.argv[1]
#               if cls[-1] == ",": cls = cls[:-1] # remove trailing comma
#               if len(cls.split(",")) != 3:
#                       print "Error reading classes, QUOTATION marks are needed!"
#                       sys.exit()
                # cls is string, e.g. "180;20; 3, 7200; 1200;60, 0;10; 120" 
#               cls = cls.replace(" ", "")
#               cls = cls.replace(",", ", ")
                #if cls[-1] != ",": cls += "," # append trailing ","
                ##distclass, durclass, pretimeclass = cls.split(",") ## splitting of classes not needed ?
                ##distclass = distclass.split(";")

#       if option == '-r':
#               reg = []
#               while len(sys.argv) > 1 and sys.argv[1][0] != '-':
#                       if sys.argv[1][-1] == ",": val = sys.argv[1][:-1]
#                       else: val = sys.argv[1]
#                       del sys.argv[1]
#                       reg.append(int(val)) # list of [lon1, lon2, lat1, lat2]
#               if len(reg) != 4: 
#                       print "ERROR on region coordinates. DO NOT use negative values in longitude"
#                       sys.exit()
#               if reg[0]-reg[1] > 180: reg[0] -= 360

# READ stationlist from raw_input
raw = ""
while 1:
        raw = raw_input("Enter network, station(, location) or EOF to end: ")

        raw = raw.upper()
        if raw == "EOF": break
        raw = raw.split(",")
        if len(raw) < 2 or raw[1] == "": # station MUST NOT be empty
                print "... invalid - not written"
                continue
        #if len(raw) < 3: raw.append('') # add third field if no location is given
        if len(raw) == 3: # reformat locationcode to stationname as staname.loccode
                raw[1] = raw[1] + "." + raw.pop().strip()
		
        stationlist.append(raw)

	print stationlist

#print stationlist
outfile = "".join(['%s.%s_' % (n.strip(),s.strip()) for n,s in stationlist]) # 'IU.KEV_II.LVZ_'
outfile += otime.strftime("%Y%m%d") + "-" + etime.strftime("%Y%m%d") # 'IU.KEV_II.LVZ_20110520-20121231'
outfile += ".input" # append .input

### avoid long filenames
if len(stationlist) > 2:
        print ""
        print "many stations -> long filename %s" % outfile
        tstring = "_" + otime.strftime("%Y%m%d") + "-" + etime.strftime("%Y%m%d") + '.input'
        print "Enter alternative output filename (%s will be appended)" % tstring
        dum = raw_input('[ENTER keeps original]: ')
        if dum != "": outfile = dum + tstring
# 
if not otime:
        print "no starttime set ... QUIT"
        sys.exit()

# open client and get Inventory
#client = Client(user="CAU Kiel")
#inv = client.getInventory(net, starttime=otime, channel=chn, 
                          #min_latitude=min(reg[2], reg[3]), max_latitude=max(reg[2], reg[3]), 
                          #min_longitude=min(reg[0], reg[1]), max_longitude=max(reg[0], reg[1]))

#print inv.keys()
ofile = open(outfile, 'w')
ofile.write("#" + " ".join(call) + "\n")

while otime < etime:
        for sta in stationlist:
#for key in inv.keys():
        #if len(key.split(".")) == 2:
                #network, station = key.split(".")
                ##print "%s, %s, %s, %s, %s, %s, %s, %s" % (otimestr, eloc, cls, network, station, chn[:-1], comp, fmt)
                ofile.write("%s, %s, %s, %s, %s, %s, %s, %s\n" % (otimestr, eloc, cls, sta[0], sta[1], chn[:-1], comp, fmt))
        otime += 86400 # add one day
        otimestr = otime.strftime("%Y, %m, %d, %H, %M, %S, 0") # string, e.g. '2008, 12, 31, 00, 00, 00, 0', ignore microseconds

ofile.close()

## file postprocessing, sort unique (since lines may be appended to existing file)
#cmd = "sort -u " + outfile + " > tmp.x; " # sort unique
# move comment lines to top
#cmd += "awk '{if ((substr($0,1,1) == \"#\" && NR == FNR) || (substr($0,1,1) != \"#\" && NR > FNR)) print $0; }' tmp.x tmp.x > " + outfile + "; "
# remove tmp file
#cmd += "rm tmp.x;"
#print cmd
#os.system(cmd)
