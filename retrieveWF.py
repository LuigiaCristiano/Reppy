
# Daten von webdc.eu runterladen und abspeichern
"""
 retrieveWF.py

functions to retrieve (and restitute) waveforms from WebDC or IRIS DMC:
         client = Client() ## WebDC client       :: must be opened before calling following functions
         iclient = iClient() ## IRIS DMC client  :: must be opened before calling following functions
         st = fetchwebdc(client, network, station, channel, location, UTCDateTime, duration, restitute=True)
         ist = fetchiris(iclient, network, station, channel, location, UTCDateTime, duration, restitute=True)
         ## calls parsestaxml(staxml, time) to process metadata in StationXML format


recent changes (newest first, please include your initials, makes it easier to trace back problems):
        20130819: (CW) added demean in merge_comp
        20130620: (CW) minor bugfix in parsestaxml
        20120515: (CW)  -- channels xH1 xH2 are rotated to xHN xHE automatically (IRIS client only!, webdc does not provide azimuths :( )
                        returns NE but not 12.
                        -- added function merge_comp to this module
        20120413: (LC)  merge the traces if the requested data are not available as a continuos record
        20120111: (CW) forward Exception in fetch* functions on getWaveform call to parent script
                  add option 'saveraw' to fetch* functions to return both instr.corrected and raw Stream
       07092012  detrend function in the instrument correction was giving error when one trace is available but as length equal to 1. detrend has so not enough samples to work.
		 I introduced an if condition on the trace length before applying the instrument correction (line 386).
       15102012   route=False for each request on webDC in order to get Italian data
"""

import numpy as np
import matplotlib.pyplot as plt

from obspy.arclink import Client
from obspy.iris import Client as iClient
from obspy.core import Stream, Trace, UTCDateTime
from obspy.signal import seisSim
from obspy.core.util import AttribDict
from lxml import etree
import urllib # URL decode, e.g. %2A to *
#import re # regexp
#from scipy.io import savemat
from scipy.signal import detrend
#from obspy.signal.rotate import rotate_NE_RT, gps2DistAzimuth
#from merge_raw import merge_comp ## added in this module
import time
import pickle, os, threading,sys
import shutil


############### FETCH WEBDC S ####################################################
def fetchwebdc(client, network, station, channel, location, t, duration,restitute='VEL', saveraw=False,routeflag=True):
     """
     Fetch data from WebDC. Client connection MUST BE OPEN as client

     instrument correction is assigned via restitute
     restitute: 'DIS', 'VEL', 'ACC' or None
         default: 'VEL'

     saveraw: BOOLEAN, if True, returns raw seismogram in addition to instr.corrected

     """

     try:
                st = client.getWaveform(network, station, location, channel, t , t+duration, metadata=True,route=routeflag)
                if (len(st)>3 and str(channel[-1])=="*"  or len(st)>1 and str(channel[-1])!="*" or
                    len(st)>1 and st[0].stats.channel== st[1].stats.channel):
                        #last condition is when we ask for * but only one component is available
                        st = merge_comp(st.copy())

     except:
             raise


     if not restitute :
             return st

     ## RESTITUTION
     stout = Stream()
     for i, tr in enumerate(st): ### LOOP CHANNELS
             res = insdecon(tr, restitute)
             stcorr=Stream([Trace(data=res, header=tr.stats)])
             stout += stcorr

     if saveraw: 
             return stout, st
     else: 
             return stout

############### FETCH WEBDC E ####################################################

############### FETCH IRIS S ####################################################
def fetchiris(client, network, station, channel, location, t, duration, restitute='VEL', saveraw=False):
        """
        Fetch data from IRIS Client connection MUST BE OPEN as iclient!
                
        instrument correction is assigned via restitute
        restitute: 'DIS', 'VEL', 'ACC' or None
        default: 'VEL'

        saveraw: BOOLEAN, if True, returns raw seismogram in addition to instr.corrected
                
        """


## TIMING
        try:
             st = client.getWaveform(network, station, location, channel, t , t+duration)
	     if (len(st)>3 and str(channel[2:3])=="*"  or len(st)>1 and str(channel[2:3])!="*" or
                 len(st)>1 and st[0].stats.channel == st[1].stats.channel):
                #last condition is when we ask for * but only one component is available
                        st = merge_comp(st.copy())

        except:
             raise

## get RESP information
        if not restitute:
             lvl='sta'
        else:
             lvl='resp'

        staxml = client.station(network=network, station=station, location=location, channel=channel, starttime=t, endtime=t+duration, level=lvl)
        # accepted levels: net, sta, chan, resp -- all else will default to sta
        stadict = parsestaxml(staxml, t) ## parse XML station file to station dictionary

        # assign channel information (coords, uptime, azimuth, paz) to tr.stats
        # here we assume that network, station AND location are explicitly given - NO wildcards
#cw     stout = Stream()
        for i, tr in enumerate(st): ### LOOP CHANNELS
             ## final adjustments to stats
             tr.stats['coordinates'] = AttribDict({'latitude': stadict[network][station]['Lat'], 'elevation': stadict[network][station]['Elevation'], \
                                                'longitude': stadict[network][station]['Lon']})
             tr.stats['uptime'] = AttribDict({'StartDate': UTCDateTime(stadict[network][station]['StartDate']), 'EndDate': UTCDateTime(stadict[network][station]['EndDate'])})
             if lvl == 'resp':
                     chn = tr.stats.channel
                     #cw if chn[-1] == '1' or chn[-1] == '2':  # keep azimuth if channels xH1 xH2 are present
                     tr.stats['azimuth'] = stadict[network][station][chn]['Azimuth']
                     if chn not in stadict[network][station].keys():
                             chn = chn + '.' + location
                             if chn not in stadict[network][station].keys():
                                     print ".......... channel info not found", chn
                                     continue
                     ### channel info found in stadict
                     tr.stats['paz'] = AttribDict({'zeros': stadict[network][station][chn]['zeros'], 'sensitivity': stadict[network][station][chn]['sensitivity'], \
                                                'poles': stadict[network][station][chn]['poles'], 'name': network + '.' + station + '.' + location + '.' + chn, \
                                                   'gain': stadict[network][station][chn]['gain'], 'input': stadict[network][station][chn]['InputUnits']})


        # handle channels xH1 xH2, if xHN xHE are not present
        avcomp = [tr.stats.channel[-1] for tr in st]
        if not ('N' in avcomp and 'E' in avcomp): # NE not available, rotate xH1 xH2 channels to xHN xHE
                c1 = c2 = np.array([])
                for i, tr in enumerate(st):
                        if tr.stats.channel[-1] == '1':
                                c1 = tr.data.copy()
                                az1 = tr.stats.azimuth
                        if tr.stats.channel[-1] == '2':
                                c2 = tr.data.copy()
                                az2 = tr.stats.azimuth
                if c1.any() and c2.any():
                        if abs(az1 - az2) != 90 and abs(az1 - az2) != 270: # latter e.g. 350 and 80 degrees
                                raise ValueError("Components 1 and 2 are not orthogonal!")
                        n, e = rotate_12_ne(c1, c2, az1)

                        nstats = st[avcomp.index('1')].stats.copy()
                        nstats.azimuth = 0.0
                        nstats.channel = nstats.channel[:-1] + 'N'
                        st += Stream([Trace(data=n, header=nstats)])

                        estats = st[avcomp.index('2')].stats.copy()
                        estats.azimuth = 90.0
                        estats.channel = estats.channel[:-1] + 'E'
                        st += Stream([Trace(data=e, header=estats)])

                        # remove components 12
                        del st[[tr.stats.channel[-1] for tr in st].index('1')]
                        del st[[tr.stats.channel[-1] for tr in st].index('2')]

        # at last: restitution
        if lvl == 'resp':
                stout = Stream()
                for i, tr in enumerate(st):
                        res = insdecon(tr, restitute)
                        stout += Stream([Trace(data=res, header=tr.stats)])

        if lvl == 'resp' and saveraw:
             return stout, st
        elif lvl == 'resp':
             return stout
        else:
             return st

############### FETCH IRIS E ####################################################

############### PARSE STAXML S ####################################################
def parsestaxml(staxml, time):
     '''
        stadict = parsestaxml(staxml, time)

        parse StationXML information to station metadata
        stadict is multi-level dictionary of order: 
                if level = 'net': network
                if level = 'sta': network -> station
                if level = 'chan': network -> station -> channel
                    !!! if more than 3 channels are found (i.e. location contains wildcards), channel key contains location code, e.g. BHZ.00
                if level = 'resp': network -> station -> channel (inkl. PAZ)

     '''


     xml_doc = etree.fromstring(staxml) ## etree Element
     ns = etree.FunctionNamespace(xml_doc.nsmap[None]) ## namespace xml_doc.nsmap -> {None: 'http://www.data.scec.org/xml/station/', 'xsi': 'http://www.w3.org/2001/XMLSchema-instance'}
     ns.prefix = 'foo'

     url=urllib.unquote(xml_doc.xpath('string(foo:ModuleURI)')) # URI 
             #['http://www.iris.edu/ws/station/query?network=II', 'level=net', 'timewindow=2004-12-26%2C2004-12-26', 'station=%2A', 'location=%2A', 'channel=%2A']
     url=url.split('?', 1)[1].split('&') # list of URI parameters: 
             # ['network=II', 'level=resp', 'timewindow=2004-12-26,2004-12-26', 'station=BFO', 'location=*', 'channel=BH*']
     # extract level from URL
     for i, item in enumerate(url):
             k,v = item.split('=')
             if k == 'level':
                     level = v
                     break


     s = dict() # return dictionary
     # add new dict per network
     for element in xml_doc.xpath('//foo:Network'): # LOOP over networks
             s[element.attrib['net_code']]=dict()

     if level == 'net': # mission accomplished
             return s
     

     ## extract station information (name, lat, lon, elevation, starttime, endtime)
     for element in xml_doc.xpath('//foo:StationEpoch'): # LOOP over stations
             if UTCDateTime(element.xpath('string(foo:StartDate)')) < time and UTCDateTime(element.xpath('string(foo:EndDate)')) > time:
                     network = element.getparent().attrib['net_code']
                     station = element.getparent().attrib['sta_code']
                     # add new dict for station information
                     s[network][station]=dict()
                     for e in element.iterchildren():
                             if e.text.find('\n') >= 0: # drop tags with no content
                                     continue
                             tag=e.tag
                             if e.tag[:1] == "{": # remove namespace from tag
                                     tag=e.tag[1:].split("}", 1)[1]

                             if 'Date' in tag:  ## format conversions to date float and ints
                                     s[network][station][tag]=UTCDateTime(e.text)
                             elif 'Lat' in tag or 'Lon' in tag or 'Elevation' in tag:
                                     s[network][station][tag]=float(e.text)
                             elif 'Number' in tag: 
                                     s[network][station][tag]=int(e.text)
                             else:
                                     s[network][station][tag]=e.text

                     if level != 'sta':
                         ## LOOP over channels
                         for e2 in element.xpath('.//foo:Epoch'): # go into each epoch (usually == channel)
                             if UTCDateTime(e2.xpath('string(foo:StartDate)')) < time and UTCDateTime(e2.xpath('string(foo:EndDate)')) > time:
                                     chn = e2.getparent().attrib['chan_code']
                                     loc = e2.getparent().attrib['loc_code']
                                     if float(s[network][station]['SelectedNumberChannels']) > 3:
                                             chn = chn + '.' + loc
                                     # add new dict for channel information
                                     s[network][station][chn]=dict()
                                     s[network][station][chn]['Location'] = loc ## add Location
                                     s[network][station][chn]['EquipType'] = e2.xpath('string(foo:Sensor/foo:EquipType)') ## add SensorType
                                     for e in e2.iterchildren():
                                             if e.text.find('\n') >= 0: # drop tags with no content
                                                     continue
                                             tag=e.tag
                                             if e.tag[:1] == "{": # remove namespace from tag
                                                     tag=e.tag[1:].split("}", 1)[1]

                                             if 'Date' in tag:  ## format conversions to date float and ints
                                                     s[network][station][chn][tag]=UTCDateTime(e.text)
                                             elif 'Lat' in tag or 'Lon' in tag or 'Elevation' in tag or \
                                                             'Azimuth' in tag or 'ClockDrift' in tag or \
                                                             'Depth' in tag or 'Dip' in tag or 'SampleRate' in tag:
                                                     s[network][station][chn][tag]=float(e.text)
                                             elif 'Number' in tag: 
                                                     s[network][station][chn][tag]=int(e.text)
                                             else:
                                                     s[network][station][chn][tag]=e.text

                                     ## extract response information -- PolesAndZeros
                                     if level == 'resp':
                                             A0 = float(e2.xpath('string(foo:Response[@stage=1]/foo:PolesZeros/foo:NormalizationFactor)')) ## 'gain' of seismometer = A0 (in SEED notation)
                                             s[network][station][chn]['gain'] = A0 ## in obspy 'gain' is A0 normalization of seismometer
                                             s[network][station][chn]['sensitivity'] = float(e2.xpath('string(foo:InstrumentSensitivity/foo:SensitivityValue)'))
                                             s[network][station][chn]['InputUnits'] = e2.xpath('string(foo:Response[@stage=1]/foo:PolesZeros/foo:InputUnits)')
                                             re = [float(i) for i in e2.xpath('foo:Response[@stage=1]/foo:PolesZeros/foo:Pole/foo:Real/text()')]
                                             im = [float(i) for i in e2.xpath('foo:Response[@stage=1]/foo:PolesZeros/foo:Pole/foo:Imaginary/text()')]
                                             s[network][station][chn]['poles'] = list()
                                             for i in range(0,len(re)):
                                                     s[network][station][chn]['poles'].append(complex(re[i], im[i]))
                     
                                             re = [float(i) for i in e2.xpath('foo:Response[@stage=1]/foo:PolesZeros/foo:Zero/foo:Real/text()')]
                                             im = [float(i) for i in e2.xpath('foo:Response[@stage=1]/foo:PolesZeros/foo:Zero/foo:Imaginary/text()')]
                                             s[network][station][chn]['zeros'] = list()
                                             for i in range(0,len(re)):
                                                     s[network][station][chn]['zeros'].append(complex(re[i], im[i]))
                         ## LOOP over channels END
                                     
             

     # LOOP over stations END

     return s

############### PARSE STAXML E ####################################################

############### CHKZEROS S ####################################################
def chkzeros(izeros, restitute):
        '''
        check number of zeros in PAZ to comply with requested restitution
        DIS, VEL, ACC need
          3,   2,   1 zero
        [0j, 0j, 0j], [0j, 0j], [0j]

        return: izeros is updated list with zeros 
        '''

        # tout is no of zeros that need to be returned
        if restitute == 'DIS':
                tout = 3
        elif restitute == 'VEL':
                tout = 2
        elif restitute == 'ACC':
                tout = 1
        else:
                #print ".... ERROR: restitute not recognized", restitute
                return [0j, 0j]

        tin = izeros.count(0j) # tin is number of zeros in input

        if tin != tout:
                while izeros.count(0j) > tout:
                        izeros.remove(0j)
                        #print " ... zero removed"
                while izeros.count(0j) < tout:
                        izeros.append(0j)
                        #print " ... zero added"
        return
############### CHKZEROS E ####################################################
     
############### INSDECON E ####################################################
def insdecon(tr, restitute):
        '''
        Instrument deconvolution wrapper for Trace objects

        tr: Trace
        restitute: 'DIS', 'VEL', 'ACC'
        '''

        if not restitute in ['DIS', 'VEL', 'ACC']:
                #print "ERROR on restitute ", restitute
                return 0

        ## RESTITUTION
        #print ".......... restitute channel", tr.stats.channel, " to ", restitute
        #if restitute == 'DIS': # apply prefilter for DIS correction
        #             tlen = (tr.stats.npts - 1) * tr.stats.delta
        #             filt = [4. / tlen, 6. / tlen, tr.stats.sampling_rate / 2., tr.stats.sampling_rate / 1.]
        #             #print "............. prefilter %ds -> %ds, %.1fHz -> %.1fHz" % (1/filt[0], 1/filt[1], filt[2], filt[3])
        #else: 
        filt = None
        ## ensure correct number of zeros in PAZ
        chkzeros(tr.stats.paz.zeros, restitute)
        if len(tr.data)>2:
          res = seisSim(detrend(tr.data), tr.stats.sampling_rate, tr.stats.paz, pre_filt=filt, nfft_pow2=True, remove_sensitivity=True)
        return res
############### INSDECON E ####################################################


############### ROTATE_12_NE S ####################################################
def rotate_12_ne(c1, c2, az):
        '''
        rotates Ch1 and Ch2 components of seismogram to North-East.

        angle is given as azimuth to Ch1.

        :param c1: data of Ch1
        :param c2: data of Ch2
        :param az: azimuth of Ch1, angle between North and orientation of Ch1
        :return: North and East component of seismogram
       
        '''

        if len(c1) != len(c2):
                raise TypeError("Components 1 and 2 have different length!?!")

        n = c1 * np.cos(az * 2 * np.pi / 360) - c2 * np.sin(az * 2 * np.pi / 360)
        e = c1 * np.sin(az * 2 * np.pi / 360) + c2 * np.cos(az * 2 * np.pi / 360)

        return n, e
############### ROTATE_12_NE E ####################################################

############### MERGE COMPONENTS S ####################################################
#merge the traces
def merge_comp(sttempf):
# 130819 (CW): demean traces before merging with fill value 0

                          z = sttempf.select(component="Z")
                          n = sttempf.select(component="N")
                          e = sttempf.select(component="E")
#cw                       z=Stream()
#cw                       e=Stream()
#cw                       n=Stream()
#cw                       for ind in range(len(sttempf)):
#cw                             if sttempf[ind].stats.channel[2:3]=="Z":
#cw                                  z +=sttempf[ind]
#cw                             if sttempf[ind].stats.channel[2:3]=="E":
#cw                                     e +=sttempf[ind]
#cw                             if sttempf[ind].stats.channel[2:3]=="N":
#cw                                     n +=sttempf[ind]
                          sttemp=Stream()
                          z.sort(['starttime'])
                          e.sort(['starttime'])
                          n.sort(['starttime'])
                          if len(z)>1: 
                                  z.detrend(type='demean')
                                  z.merge(method=1,fill_value=0, interpolation_samples=0)
                          if len(e)>1: 
                                  e.detrend(type='demean')
                                  e.merge(method=1, fill_value=0, interpolation_samples=0)
                          if len(n)>1: 
                                  n.detrend(type='demean')
                                  n.merge(method=1,fill_value=0, interpolation_samples=0)

                          sttemp +=z
                          sttemp +=e
                          sttemp +=n
			
                          return sttemp

############### MERGE COMPONENTS E ####################################################
