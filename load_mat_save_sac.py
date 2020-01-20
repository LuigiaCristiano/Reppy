
###load .mat file and save them in SAC format
###(Requirements: Python 3 and Obspy 1)

#path is the first input parameter
#filename is the second input parameter
import sys
reload(sys)
path=str(sys.argv[1])
filename=str(sys.argv[2])
filename_new=filename.replace(".mat",".sac")
from loadmat2trace import loadmat2trace

st=loadmat2trace([path+filename])
st.write(filename_new,format="SAC"
