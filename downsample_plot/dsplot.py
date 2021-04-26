import numpy as np
import matplotlib.pyplot as plt
import os
import struct
#from tqdm import tqdm
import argparse
import sys

def plot_downsamp(avgfac,samps,nchan,data):
	avgsampuse = int(np.floor(samps/avgfac))
	timmy = np.zeros((nchan,avgsampuse))
	for tavg in np.arange(avgsampuse):
		subdata = data[:,tavg*avgfac : (tavg+1)*avgfac]
		subavg = np.reshape(np.sum(subdata,axis=1), (nchan,1)) / avgfac
		timmy[:,tavg] = subavg[:,0]
	return timmy
	
def grab_data(FILE,STARTSAMP,NUMSAMP,NCHAN,DTYPE,FIL):
    with open(FILE,'r') as F:
        if not FIL:
            headlen = 4096
            #print "hi there, I'm assuming the DADA header is 4096 bits"
        else:
            #print "hi there, I'm reading a filterbank header now"
            thehead = header(F)
            headlen = len(thehead)
        #F.seek(4096+NCHAN*STARTSAMP)
        F.seek(headlen+NCHAN*STARTSAMP)
        #data = np.fromfile(F,dtype=np.uint8,count=int(NCHAN*NUMSAMP))
        data = np.fromfile(F,dtype=DTYPE,count=int(NCHAN*NUMSAMP))
        data = np.reshape(data,(-1,NCHAN)).T
    return data

def header(afile):
    inread = ""
    while True:
        tmp = afile.read(1)
        inread = inread + tmp 
        flag = inread.find('HEADER_END')
        if flag != -1: 
            break
    return inread

def get_headparam(head,parlist):
    # grabbing header values
    how2parse={
        'nchans': ('i',4),
        'tsamp': ('d',8),
        'foff': ('d',8),
        'fch1': ('d',8),
        'tstart': ('d',8),
        'ibeam': ('i',4),
        'nbits': ('i',4)
    }   
    n=0 
    for i in parlist:
        i1 = head.find(i)
        i2 = i1 + len(i)
        nbytes = how2parse[i][1]
        cstr = how2parse[i][0]

        val = struct.unpack(cstr, head[i2:i2+nbytes])[0]
        parlist[n] = val 
        n += 1
        return parlist

def ds_n_plot(f,ds,nchan,start_counter,endsamp,nsamps,tsamp,outname,show):
	# set counter and read-in block size
	counter = start_counter 
	blocksize=int(np.ceil(ds))
	# initiate outplot array
	totdown = nsamps/ds
	outdat = np.zeros((nchan,totdown))	
	# set progres bar
	bar_length = 30
	# loop through data and downsample
	for i in range(np.shape(outdat)[1]):
		start_tmp = counter
		end_tmp = np.minimum(start_tmp + blocksize,endsamp)
		n_tmp = end_tmp - start_tmp
		dat = grab_data(f,start_tmp,n_tmp,nchan,np.uint8,True)
		pltdat = plot_downsamp(ds,np.shape(dat)[1],nchan,dat)
		outdat[:,i] = pltdat[:,0]
		pc = 100.*i/np.shape(outdat)[1]
		sys.stdout.write('\r')
		sys.stdout.write("Progress: [{:{}}] {:>3}%"
					.format('='*int(pc/(100./bar_length)),
					bar_length,int(pc)))
		sys.stdout.flush()
		counter += blocksize
	fig,ax = plt.subplots()
	ax.imshow(outdat,aspect='auto',interpolation='none',cmap='afmhot_r')
	ax.set_ylabel('Channel')
	ttt = np.round(ds*tsamp,2)
	ax.set_xlabel("Sample (downsamp factor {}, 1 sample={} sec)".format(ds,ttt))
	ax.set_title("{}".format(os.path.basename(f)))
	fig.tight_layout()
	if show:
		plt.show()
	plt.savefig(outname)


# unfortunately get rid of tqdm, maybe make a shitty progress bar instead
# add options
# 	input file
#	downsampling factor 
#	start, stop, all file
#	here or in a separate script, FFT a single channel



if __name__ == "__main__":
	desc = """
		downsample and plot filterbank
			"""
	parser = argparse.ArgumentParser(description=desc)
	parser.add_argument('--fil',type=str,help='filterbank file')
	parser.add_argument('--downsamp',type=int,default=64,
			help="downsampling factor (power of 2), default=64")
	parser.add_argument('--start',type=float,default=0.,
			help="starting time in seconds, default=0")
	parser.add_argument('--end',type=float,default=np.inf,
			help="end time in seconds, default=end of file")
	parser.add_argument('--out',type=str,default=os.getcwd(),
			help="plot directory,default=current dir")
	parser.add_argument('--show',action='store_true',default=False,
			help="show plot")
	args = parser.parse_args()

	# find start/end of plot
	st = np.maximum(0,args.start)
	en = np.minimum(np.inf,args.end)

	# outfile path and name
	basename = os.path.splitext(os.path.basename(args.fil))[0]
	outname_base = "{}_DS_{}.png".format(basename,args.downsamp)
	outname = os.path.join(args.out,outname_base)
	
	# get file params
	with open(args.fil,'r') as F:
		head = header(F)
	tsamp = get_headparam(head,['tsamp'])[0]
	nchan = get_headparam(head,['nchans'])[0]

	# size, samples, length of obs, downsampled samples
	si = os.path.getsize(args.fil)
	si -= len(head)
	totsamps = si/nchan
	obslen = totsamps * tsamp

	# start/stop conditions
	starty = int(st/tsamp)
	endy = en/tsamp
	if endy > totsamps:
		endy = int(np.floor(totsamps))
	else:
		endy = int(endy)
	numsamps = endy - starty

	ds_n_plot(args.fil,args.downsamp,nchan,starty,endy,numsamps,tsamp,outname,args.show)





