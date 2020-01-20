from scipy.io import loadmat
from obspy.core import AttribDict, UTCDateTime, Trace,Stream
import numpy as np

# scipy version handling
from scipy import __version__ as sv
scipyver = int(sv.split(".")[1])

def loadmat2trace(infile, debug=False):
    '''
    wrapper function for scipy.io.loadmat
    reformats seismogram read from .mat file to ObsPy Trace object
    including all available tr.stats entries
    '''
    ## conversion libraries
    tofloat = ['azimuth', 'calib', 'delta', 'sampling_rate']
    toUTCDateTime = ['starttime']
    toremove=['__header__','__globals__', '__version__', 'endtime'] # endtime is read-only
    st = loadmat(infile, squeeze_me=True)

    if debug:
        print "1) after loadmat"
        print st

    # remove unwanted keys
    for i in st.keys():
        if i in toremove: del st[i]

    # remove Unicode from strings
    for i in st.keys():
        if scipyver == 10: # works on scipy 0.10 on ifg machines
            if i[:2] != "__" and st[i].dtype.kind == "U":  st[i] = st[i].item().encode()
        elif scipyver >= 12: # works on scipy 0.12 on rzcluster
            if type(st[i]) in [unicode, str]:
                st[i] = st[i].encode() # Strings
            elif type(st[i]) is np.ndarray:
                st[i] = st[i].astype(np.float) # NumPy Arrays
    if debug:
        print ""
        print "2) after unicode removal"
        print st

    # convert strings to AttribDicts, floats, UTCDateTime
    for i in st.keys():
        if st[i][:10] == "AttribDict" or i in ['processing']: # AttribDict or list
            st[i] = eval(st[i])
        elif i in tofloat:
            st[i] = float(st[i])
        elif i in toUTCDateTime:
            st[i] = UTCDateTime(st[i])
        elif i == 'location' and len(st[i]) == 0: # location == ''
            st[i] = ''

    if debug:
        print ""
        print "3) after data type assignment"
        print st

    ## create Trace
    tr = Trace(data=st['data']) # adds default header, starttime 1970-1-1-0:0:0, sampling rate 1
    for i in st.keys():
        if i == 'data': continue
        tr.stats[i] = st[i]

    if debug:
        print ""
        print "4) new trace"
        print tr.stats
        print tr.data

    return tr #, st :: is that needed?

if __name__ == '__main__':

    #tr = loadmat2trace(infile, debug=True)
    tr = loadmat2trace('/rayleighdata/data/2009/200909301016/GR.GRA1.BHZ_VEL_2009-09-30.10-16-08_9600.mat')
    print tr
    print tr.stats
