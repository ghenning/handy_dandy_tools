import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
import struct
from tqdm import tqdm
import random

def header(afile):
	# read the header
	inread = ""
	while True:
		tmp = afile.read(1)
		inread = inread + tmp
		flag = inread.find('HEADER_END')
		if flag != -1:
			break
	return inread

def grab_data(FILE,STARTSAMP,NUMSAMP,NCHAN,DTYPE):
	# grabbing data
	with open(FILE,'r') as F:
		thehead = header(F)
		headlen = len(thehead)
		F.seek(headlen+NCHAN*STARTSAMP)
		data = np.fromfile(F,dtype=DTYPE,count=int(NCHAN*NUMSAMP))
		data = np.reshape(data,(-1,NCHAN)).T
	return data

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

def scaleit(DATA,MEAN):
	# divide each channel by median of data (or 1 if median is 0)
	###bandpass_shape = np.median(DATA,1).reshape((-1,1))
	bandpass_shape = np.mean(DATA,1).reshape((-1,1))
	# take mean of current mean and "running" mean
	bp_mean = np.mean(np.array([MEAN,bandpass_shape]),axis=0)
	###tmpdat = DATA/np.maximum(bandpass_shape,1)
	tmpdat = DATA/np.maximum(bp_mean,1)
	# divide everything by data mean and standard deviation
	tmpdat -= tmpdat.mean()
	tmpdat /= tmpdat.std()
	# clip data between -2 and 4 and re-scale to 8 bit
	tmpdat = np.clip(tmpdat,-2.,4.) + 2.
	outdat = np.uint8(tmpdat * 256 /6)
	return outdat, bp_mean

def read_write_scale(FIL,OUTFIL,NCHAN,TSAMP,START,END):
	# grab header of original file
	with open(FIL,'r') as X:
		print "working on file {}".format(FIL)
		thehead = header(X)
		headlen = len(thehead)
	# write header of new file
	with open(OUTFIL,'w') as F:	
		print "writing header to new fil {}".format(OUTFIL)
		F.write(thehead)		
	print "scaling..."
	# loopy loop through the file and scale
	blocksize = 2000
	seeker = headlen + START
	# grab random part of data to start "running smoothing mean"
	where = int(np.random.uniform(.2,.8)*END/NCHAN)
	randdata = grab_data(FIL,where,blocksize,NCHAN,np.uint8)
	randmean = np.mean(randdata,1).reshape((-1,1))
	with tqdm(total=END-seeker) as pbar:
		while seeker<END:
			# grabbing data and sending it to scaleit
			start_tmp = seeker/NCHAN
			end_tmp = np.minimum(start_tmp + blocksize,END/NCHAN)
			n_tmp = end_tmp - start_tmp
			dataIn = grab_data(FIL,start_tmp,n_tmp,NCHAN,np.uint8)
			seeker += len(dataIn.flatten())
			dataIn, randmean = scaleit(dataIn,randmean)
			dataIn = dataIn.T
			with open(OUTFIL,'a') as F:
				F.write(dataIn.flatten())
			pbar.update(len(dataIn.flatten()))

if __name__ == "__main__":
	desc = """
		An attempt (work in progress) to scale data that
		have an uneven bandpass. Currently: only works
		for filterbanks. 
		Scale data or plot a dynamic spectrum window to
		take a look at the data. 
		To scale, use --fil, to scale a slice of the 
		original data use --start --end.
		To plot, use --fil --show --time, and optionally
		use --dur to set the size of the window.
		Lots of TODOs, will update when I feel like it.
		"""
	parser = argparse.ArgumentParser(description=desc)
	parser.add_argument('--fil',type=str,help='filterbank file')
	parser.add_argument('--show',action='store_true',
		help='plot dynspec window',default=False)
	parser.add_argument('--time',type=float,
		help='plot window start (in sec). default=0',default=0)
	parser.add_argument('--dur',type=float,
		help='how big a plot window (in sec). default=0.1',default=0.1)
	parser.add_argument('--start',type=float,
		help='where to start scaling file (in sec). default=0',default=0)
	parser.add_argument('--end',type=float,
		help='where to end scaling file (in sec). default=end',default=np.inf)
	args = parser.parse_args()

	# find start/end of scaling
	st = np.maximum(0,args.start)
	en = np.minimum(np.inf,args.end)

	# set outfile path and name
	# TODOs: add start and stop time in name if start/end flags are used
	basename = os.path.splitext(os.path.basename(args.fil))[0]
	basename = "{}_scaled.fil".format(basename)
	thedir = os.path.dirname(args.fil)
	totalout = os.path.join(thedir,basename)

	# get file parameters
	with open(args.fil,'r') as F:
		head = header(F)
	tsamp = get_headparam(head,['tsamp'])[0]
	nchan = get_headparam(head,['nchans'])[0]

	# find size and observing length of file
	si = os.path.getsize(args.fil)
	si -= len(head)
	print "data size {} bytes".format(si)
	totsamps = si / nchan	
	print "total samples {}".format(totsamps)
	obslen = totsamps * tsamp
	print "obs len {:.2f} seconds".format(obslen)

	# will probably split this into two py scripts:
	#	one for scaling
	#	one for plotting
	#		- dynspec
	#		- data histograms
	#		- compare two files (e.g. scaled/not scaled)
	#	add axes labels
	# plot dynamic spectrum
	if args.show:
		startsamp = int(args.time / tsamp)
		nsamp = int(args.dur / tsamp)
		pltdata = grab_data(args.fil,startsamp,nsamp,nchan,np.uint8)
		fig,ax = plt.subplots()
		ax.imshow(pltdata,aspect='auto',interpolation='hamming',cmap='jet')
		plt.show()
	else:
		# check start/stop conditions
		starty = int(st / tsamp)
		endy = en / tsamp
		if endy > totsamps:
			endy = int(np.floor(totsamps))
		else:
			endy = int(endy)
		starty *= nchan
		endy *= nchan
		# scale the data
		read_write_scale(args.fil,totalout,nchan,tsamp,starty,endy)
	
