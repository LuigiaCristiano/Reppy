###from obspy web page

from obspy import read

from obspy.io.xseed import Parser

from obspy.signal import PPSD
import d2fcts_mod as d2f
import sys
reload(sys)

from loadmat2trace import loadmat2trace
#OPEN INPUT FILE AND READ PARAMETERS
InputFile = str(sys.argv[1]); del sys.argv[1]

#####we read the list of stations we want to process and calculate the ppsd for the given number of days
######Input1 station list
######Input 2 number of days
data_down_infos=csv.reader(CommentStripper(open(InputFile,"r")), skipinitialspace=True) # allows for blanks in input line
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Read data and select a trace with the desired station/channel combination:

### LOOP 1: input file line by line
for line in data_down_infos: # line is list of elements in input line
[torigin, event_lat, event_lon, depth, mag, path_data, dist_class, duration_class, toll_class, \
rotation, correction, Network, Station, Location, Channel, sub, dataformat] = d2f.read_input_file(line)
######search for the  data in the database
available, wrongfmt, rotateme = d2f.check_data_avail(tstart, Network, Stationr[inds], reqchan, laufzeit, correction, path_data, dataformat)
    ppsd = PPSD(tr.stats, metadata=parser)

    ppsd.add(st)
    print("number of psd segments:", len(ppsd.times))
    ppsd.plot()
    ppsd.plot("/tmp/ppsd.png")  
    ppsd.plot("/tmp/ppsd.pdf")  
