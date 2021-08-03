import numpy as np
import psrchive
import matplotlib.pyplot as plt
import argparse
import time
from subprocess import Popen,PIPE
import os

# Setting up the font and fontsize
rc = {"font.family" : "serif","mathtext.fontset" : "stix"}
plt.rcParams.update(rc)
plt.rcParams["font.serif"] = ["Times New Roman"] + plt.rcParams["font.serif"]
FS = 20
FS2 = 22
FS3 = 12

# Get all the archive info
def tst(A,DS,SB,W,C,R):
	# Load the archive, convert to Stokes, remove baseline
	ar = psrchive.Archive_load(A)
	if ar.get_state()=='Coherence':
		ar.convert_state('Stokes')
	ar.remove_baseline()

	# Downsample, subband, find zapped channels
	if DS > 1:
		ar.bscrunch(DS)
	if SB > 1:
		ar.fscrunch(SB)
	if W:
		w = ar.get_weights()
	else:
		w = np.ones(ar.get_nbin())

	# Archive parameters
	bw = ar.get_bandwidth()
	nbin = ar.get_nbin()
	nchan = ar.get_nchan()
	fchan = abs(bw)/ar.get_nchan()
	flo = ar.get_centre_frequency() - abs(ar.get_bandwidth()/2.)
	fhi = ar.get_centre_frequency() + abs(ar.get_bandwidth()/2.)
	length = ar.integration_length()
	
	# for pulsars - optionize this
	###boo = ar.get_Integration(0)
	###length = boo.get_folding_period()

	# move data to center
	if C:
		ar.centre_max_bin()
	### add option to rotate 
	if R != 0:
		ar.rotate(R)

	# find the edges of the burst 
	# how to get this to work?
	#print ar.find_max_phase()
	#print ar.find_min_phase()
	##print ar.find_peak_edges(0,10)
	#print ar.find_transitions(80,10,3)

	# grab the archive data
	thedat = ar.get_data()

	# grab Stokes params
	I = thedat[0][0][:][:]
	Q = thedat[0][1][:][:]
	U = thedat[0][2][:][:]
	V = thedat[0][3][:][:]

	# Apply zaps
	# optionize this (default True)
	I = (I.T * w).T
	Q = (Q.T * w).T
	U = (U.T * w).T
	V = (V.T * w).T

	# Testing some time series
	# It's a little bit different than what the ar time series is
	#I_TS = np.average(I,axis=0,weights=np.squeeze(w))
	#TSI = np.sum(I,axis=0)/np.count_nonzero(w)
	#TSI /= np.amax(TSI)
	#TSQ = np.sum(Q,axis=0)/np.count_nonzero(w)
	#TSQ /= np.amax(TSQ)
	#TSU = np.sum(U,axis=0)/np.count_nonzero(w)
	#TSU /= np.amax(TSU)
	#TSV = np.sum(V,axis=0)/np.count_nonzero(w)
	#TSV /= np.amax(TSV)
	#TSL = np.sqrt(TSQ**2 + TSU**2)
	

	# Reduce data to time series
	# I don't know whether this takes into account zapped channels
	# perhaps manually reduce I,Q,U,V that have been zapped 
	ar.fscrunch_to_nchan(1)
	TSI = ar.get_data()[0][0][0]
	TSQ = ar.get_data()[0][1][0]
	TSU = ar.get_data()[0][2][0]
	TSV = ar.get_data()[0][3][0]
	TSL = np.sqrt(TSQ**2 + TSU**2)

	return TSI,TSV,TSL,I,Q,U,V,flo,fhi,fchan,TSQ,TSU,nbin,bw,nchan,length


if __name__=="__main__":
	desc = "trying to use something other than psrplot"
	parser = argparse.ArgumentParser(description=desc)
	optional = parser._action_groups.pop()
	required = parser.add_argument_group('required args')
	required.add_argument('--data',type=str,help='Burst archive')
	optional.add_argument('--todo',type=str,help='TBD')
	optional.add_argument('--downsample',type=int,default=1,
		help="Downsampling factor, default=1 (no downsampling)")
	optional.add_argument('--subband',type=int,default=1,
		help="Subbanding factor, default=1 (no subbanding)")
	optional.add_argument('--weights',action='store_false',default=True,
		help="Apply pazi zap info. Default=True")
	optional.add_argument('--cmap',type=str,default='afmhot',
		help="Plot colormap. Default=afmhot")
	optional.add_argument('--zapcol',type=str,default='turquoise',
		help="Color of zapped channels. Default=turquoise")
	optional.add_argument('--center_max',action='store_false',default=True,
		help="Center the max bin. Default=True")
	#optional.add_argument('--print_stats',action='store_true',default=False,
		#help="Print pol stats. Default=False")
	optional.add_argument('--blabel',type=str,default='tst',
		help="Burst number")
	optional.add_argument('--notime',action='store_true',default=False,
		help="No time (x-axis) label")
	optional.add_argument('--nofreq',action='store_true',default=False,
		help="No frequency (y-axis) label")
	optional.add_argument('--nopa',action='store_true',default=False,
		help="No PA (PA y-axis) label")
	optional.add_argument('--rotate',type=float,default=0,
		help="Rotate profile by X seconds")
	optional.add_argument('--Ilim',type=float,default=0.4)
	optional.add_argument('--PAout',type=str,
		help="Directory to save PA info")
	optional.add_argument('--Iout',type=str,
		help="Directory to save stokes I")
	optional.add_argument('--special',action='store_true',default=False)
	optional.add_argument('--Nticks',type=int,default=6,
		help="yaxis.set_major_locator(plt.MaxNlocator(Nticks)),default=6")
	optional.add_argument('--savepan',type=str,
		help='Directory to save panel plot')
	optional.add_argument('--speciallab',action='store_true',default=False)
	optional.add_argument('--shiftPApos',action='store_false',default=True,
		help="Shift negative PA values to positive, then subtract 90 deg")
	optional.add_argument('--xlim',action='store_true',default=False,
		help="set xlim to 50 ms for longer ars")
	optional.add_argument('--print_stats',type=str,
		help='Directory to save pol stats')
	optional.add_argument('--allout',type=str,
		help="Directory to save all stokes params")
	optional.add_argument('--dynV',type=str,
		help="Directory to save Stokes V time-freq")
	optional.add_argument('--phi1',type=float,
		help="start phase")
	optional.add_argument('--phi2',type=float,
		help="end phase")
	# manually select on/off pulse 
	# save PA as a file
	# choose to display axis labels
	parser._action_groups.append(optional)
	args = parser.parse_args()

	A = args.data
	# I 	: Intensity time series
	# Q		: Stokes Q time series
	# U		: Stokes U time series
	# L 	: Linear polariztion time series
	# V 	: Circular polarization time series 
	# dyn	: Frequency-time intensity data
	# flo	: Frequency lower bound
	# fhi	: Frequency upper bound
	# fch	: Channel bandwidth
	# nbin	: Number of time bins
	# more stuff ... get dynspec of Q,U,V
	I,V,L,dyn,dynQ,dynU,dynV,flo,fhi,fch,QQ,UU,nbin,bw,nchan,length = tst(A,args.downsample,args.subband,args.weights,args.center_max,args.rotate)
	# flip the band if bandwidth is positive
	if bw > 0: 
		dyn = np.flipud(dyn)
		dynQ = np.flipud(dynQ)
		dynU = np.flipud(dynU)
		dynV = np.flipud(dynV)

	tstplots = False

	# some testing on not plotting edges
	#fhi -= 26
	f0 = int(26*fch)
	#flo += 31
	f1 = 512 - int(49*fch)
	###dyn[512-50:] = 0

	# set up Hi-res plot
	fig,ax = plt.subplots()

	# summin each channel 
	# was using this to inspect stuff 
	sums = np.sum(abs(dyn),axis=1)/np.amax(dyn)
	sums -= np.median(sums)
	###dyn[np.where(sums<(-np.std(sums)))] = 0

	# plot the data
	ax.imshow(dyn,aspect='auto',interpolation='none',cmap=args.cmap)#,extent=(0,1,flo,fhi))
	# color in the zapped channels
	zaps = np.squeeze(np.where(np.sum(dyn,axis=1)==0))
	for z in zaps:
		ax.axhspan(z-.5,z+.5,color=args.zapcol)

	# setting up the multi-panel plot
	#fig,ax = plt.subplots(3,1,sharex=True,gridspec_kw={'height_ratios':[1,3,8]},figsize=(6.4,4.8))
	fig,ax = plt.subplots(3,1,sharex=True,gridspec_kw={'height_ratios':[2,4,9]})
	#fig.subplots_adjust(hspace=0)

	# plot dynamic spectrum
	ax[2].imshow(dyn,aspect='auto',interpolation='none',cmap=args.cmap)#,extent=(0,1,flo,fhi))
	# color in zapped channels
	zaps = np.squeeze(np.where(np.sum(dyn,axis=1)==0))
	for z in zaps:
		ax[2].axhspan(z-.5,z+.5,color=args.zapcol)

	#phi = .5 * np.arctan(UU/QQ)
	phi = .5 * np.arctan2(UU,QQ)
	if args.shiftPApos:
		##ophi = phi[on_phi]
		##neg = np.where(ophi < 0)	
		##ophi[neg] = ophi[neg] + np.pi
		##ophi -= np.pi/2
		neg = np.where(phi < 0)
		phi[neg] = phi[neg] + np.pi
		phi -= np.pi/2
	phi = np.rad2deg(phi)
	P2 = np.square(QQ) + np.square(UU)
	#######################################
	# need a better way to determine on/off
	# doing it manually for now
	#######################################
	#if args.phi1:
	#	off1 = np.arange(0,int(args.phi1*nbin/args.downsample))
	#	off2 = np.arange(int(args.phi2*nbin/args.downsample),nbin/args.downsample)
	#else:
	#	off1 = np.arange(0,int(.25*nbin/args.downsample))
	#	off2 = np.arange(int(.75*nbin/args.downsample),nbin/args.downsample)
	if args.phi1:
		off1 = np.arange(0,int(args.phi1*nbin))
		off2 = np.arange(int(args.phi2*nbin),nbin)
	else:
		off1 = np.arange(0,int(.25*nbin))
		off2 = np.arange(int(.75*nbin),nbin)
	stdI = np.std(I[np.r_[off1,off2]])
	Inorm = I/np.amax(I)
	stdInorm = np.std(Inorm[np.r_[off1,off2]])
	Ilim = 6
	Ilim = args.Ilim * np.max(Inorm)/stdInorm
	#Ilim = .6 * np.max(Inorm)/stdInorm
	on_phi = np.where(Inorm/stdInorm > Ilim)
	off_phi = np.where(Inorm/stdInorm < Ilim)

	# unbiased L
	tmpL = np.copy(L)
	yes = np.where(tmpL/stdI >= 1.57)
	no = np.where(tmpL/stdI < 1.57)
	tmpL[yes] = stdI * np.sqrt( (tmpL[yes]/stdI)**2 -1)
	tmpL[no] = 0

	# plotting Stokes I,L,V
	ax[1].plot(I,color='k',lw=.7)
	#ax[1].plot(L,color='r',lw=.7)
	ax[1].plot(tmpL,color='r',lw=.7)
	#ax[1].axvline(off1[-1])
	#ax[1].axvline(off2[0])
	ax[1].plot(V,color='b',lw=.7)

	dQ = np.std(QQ[off_phi])
	dU = np.std(UU[off_phi])
	dphi = np.sqrt((np.abs(QQ)/P2)*dQ**2 + (np.abs(UU)/P2)*dU**2)
	dphi = np.rad2deg(dphi)
	ax[0].errorbar(np.arange(len(phi))[on_phi],phi[on_phi],yerr=dphi[on_phi],fmt='.',elinewidth=.8,ms=1,color='k')
	for tick in ax[0].yaxis.get_major_ticks():
		tick.label.set_fontsize(FS)
	for tick in ax[0].xaxis.get_major_ticks():
		tick.label.set_fontsize(FS)
	for tick in ax[1].yaxis.get_major_ticks():
		tick.label.set_fontsize(FS)
	for tick in ax[1].xaxis.get_major_ticks():
		tick.label.set_fontsize(FS)
	for tick in ax[2].yaxis.get_major_ticks():
		tick.label.set_fontsize(FS)
	for tick in ax[2].xaxis.get_major_ticks():
		tick.label.set_fontsize(FS)
	ax[1].tick_params(axis='y',which='both',labelleft=False,left=False)
	ax[1].tick_params(axis='x',which='both',bottom=False)
	if args.notime:
		print 'no xlabel'
	else:
		ax[2].set_xlabel('Time (ms)',fontsize=FS2)
	if args.nofreq:
		print 'no ylabel'
		ax[2].tick_params(axis='y',which='both',labelleft=False)
		#ax[2].tick_params(axis='y',which='both',labelcolor='white')
		#ax[2].set_ylabel('Frequency (MHz)',fontsize=FS2,color='white')
	else:
		#ax[2].set_ylabel('Frequency (MHz)',fontsize=FS2)
		ax[2].set_ylabel('Freq. (MHz)',fontsize=FS2)
	ax[1].set_xlim([0,len(I)-1])
	#print on_phi
	PAmed = int(np.median(phi[on_phi]))
	#PAmed = int(np.mean(phi[on_phi]))
	#PAmed = int(np.median(phi))
	#ax[0].set_ylim([PAmed-30,PAmed+30])
	ax[0].set_ylim([PAmed-20,PAmed+20])
	#ax[0].set_yticks([PAmed-15,PAmed,PAmed+15])
	ax[0].set_yticks([PAmed-10,PAmed,PAmed+10])
	ax[0].set_yticklabels(['',PAmed,''])
	if args.nopa:
		print 'no PA label'
	else:
		ax[0].set_ylabel('PA',fontsize=FS2)

	# get rid of band edges
	# add this as an option later
	# grab original and new number of channels to fix this
	topcut = 38.5 # for 512 chan
	botcut = 462.5 # for 512 chan
	#topcut /= args.subband
	#botcut /= args.subband
	topcut /= (512./nchan)
	botcut /= (512./nchan)
	tst = .5 * np.floor(2*topcut)
	tstt = .5 * np.ceil(2*botcut)
	botcut -= botcut%.5
	topcut += topcut%.5
	#ax[2].set_ylim([0.90332*nchan,0.075195*nchan])
	ax[2].set_ylim([botcut,topcut])

	# change from channel number to frequency
	ax[2].yaxis.set_major_locator(plt.MaxNLocator(args.Nticks))
	ticks = ax[2].get_yticks()
	print "orig yticks {}".format(ticks)
	print "fhi {} fch {}".format(fhi,fch)
	print "nchan {} bw {}".format(nchan,bw)
	nyticks = (fhi - fch * ticks).astype(int)
	print "new yticks {}".format(nyticks)
	ax[2].set_yticklabels(nyticks)

	# change from time bin to time in ms
	ticks = ax[2].get_xticks()
	ticktime = int((10e-3 * nbin)/length)
	print "bins per 10 ms {} nbins {}".format(ticktime,nbin)
	cbin = int(nbin/2)
	if args.xlim:
		ax[2].set_xlim([cbin-2.5*ticktime,cbin+2.5*ticktime])
	#nticks = np.arange(cbin,nbin,ticktime)
	nticks = np.arange(cbin,cbin+2*ticktime+1,ticktime)
	#nticks2 = np.sort(np.arange(cbin,0,-ticktime))
	nticks2 = np.sort(np.arange(cbin,cbin-2*ticktime-1,-ticktime))
	nxticks = np.unique(np.r_[nticks2,nticks])
	ax[2].set_xticks(nxticks)
	nticklabs = np.round(1e3*(nxticks/float(nbin))*length,1)
	nticklabs -= nticklabs[len(nticklabs)/2]
	nticklabs = np.round(nticklabs,0)
	ax[2].set_xticklabels(nticklabs.astype(int))
	if args.special:
		nticks = np.arange(cbin,cbin+5*ticktime+1,ticktime)
		nticks2 = np.sort(np.arange(cbin,cbin-5*ticktime-1,-ticktime))
		nxticks = np.unique(np.r_[nticks2,nticks])
		ax[2].set_xticks(nxticks)
		nticklabs = np.round(1e3*(nxticks/float(nbin))*length,1)
		nticklabs -= nticklabs[len(nticklabs)/2]
		nticklabs[0] = -50
		nticklabs[1] = -40
		nticklabs[-1] = 50
		nticklabs[-2] = 40
		nticklabs = np.round(nticklabs,0)
		print nticklabs
		ax[2].set_xticklabels(nticklabs.astype(int))
		ax[2].set_xlim([int(.26*nbin),int(.75*nbin)])
		ax[2].set_xlim([int(.23*nbin),int(.78*nbin)])
		#ax[2].set_xlim([int(.15*nbin),int(.95*nbin)])
		ax[0].set_ylim([PAmed-60,PAmed+60])
		ax[0].set_ylim([5,85])
		#ax[0].set_yticks([PAmed-10,PAmed,PAmed+10])
		ax[0].set_yticks(np.arange(5,86,20))
		#ax[0].set_yticklabels(['',PAmed,''])
		ax[0].set_yticklabels(['','','45','',''])
		

	# time and freq res
	tres = float(length)/nbin
	#fres = float(bw)/nbin
	fres = fch
	treslab = r"{:.0f} $\mu$s".format(tres*1e6)
	freslab = r"{:.2f} MHz".format(fres)
	bothlab = r"{:.0f} $\mu$s{}{:.2f} MHz".format(tres*1e6,'\n',fres)
	burstlab = "$\#${}".format(args.blabel)
	#burstlab = "$\#${}".format('tst')
	#ax[1].text(.01,.75,treslab,fontsize=FS,transform=ax[1].transAxes)
	#ax[1].text(.01,.5,freslab,fontsize=FS,transform=ax[1].transAxes,bbox=dict(facecolor='white',alpha=.5,edgecolor='none'),ha='left',va='top')

	if args.speciallab:
		ax[2].text(.015,.95,bothlab,fontsize=FS,transform=ax[2].transAxes,bbox=dict(facecolor='white',alpha=.9,edgecolor='none'),ha='left',va='top')
	else:
		ax[1].text(.015,.92,bothlab,fontsize=FS,transform=ax[1].transAxes,bbox=dict(facecolor='white',alpha=.7,edgecolor='none'),ha='left',va='top')
		
	ax[1].text(.99,.92,burstlab,fontsize=FS,transform=ax[1].transAxes,ha='right',va='top')
	if args.special:
		blab2 = "$\#${}".format(10)
		ax[1].text(.48,.92,blab2,fontsize=FS,transform=ax[1].transAxes,ha='right',va='top')
	

	#ax[0].set_xlim([0,len(phi)])
	#fig.tight_layout(h_pad=-1)
	fig.tight_layout(h_pad=-.5)
	#fig.tight_layout()
	if args.savepan:
		bname = os.path.splitext(os.path.basename(args.data))[0]
		oname = "{}.png".format(bname)
		outpath = os.path.join(args.savepan,oname)
		print outpath
		fig.savefig(outpath,dpi=200)

	if not (args.savepan or args.PAout or args.print_stats or args.Iout):
		plt.show()

	# I profile and dynamic spectrum (for hires plots)
	# add as option
	fig,ax = plt.subplots(2,1,sharex=True,gridspec_kw={'height_ratios':[1,3]})
	fig.subplots_adjust(hspace=0)
	# dynamic spectrum
	ax[1].imshow(dyn,aspect='auto',interpolation='none',cmap=args.cmap)
	# color in zapped channels
	zaps = np.squeeze(np.where(np.sum(dyn,axis=1)==0))
	for z in zaps:
		ax[1].axhspan(z-.5,z+.5,color=args.zapcol)
	# intensity
	ax[0].plot(I,color='k',lw=.7)
	ax[1].set_xlim([0,len(I)-1])

	ax[0].tick_params(axis='y',which='both',labelleft=False,left=False)
	ax[0].tick_params(axis='x',which='both',bottom=False)
	#ax[1].set_ylim([0.90623*nchan,0.07422*nchan])
	ax[1].set_ylim([botcut,topcut])
	ax[1].set_xlabel('Time (ms)',fontsize=FS2)
	ax[1].set_ylabel('Frequency (MHz)',fontsize=FS2)

	# change x and y ticks
	ax[1].yaxis.set_major_locator(plt.MaxNLocator(7))
	ax[1].set_xticks(nxticks)
	ax[1].set_xticklabels(nticklabs)
	ax[1].set_yticklabels(nyticks)

	for tick in ax[1].yaxis.get_major_ticks():
		tick.label.set_fontsize(FS)
	for tick in ax[1].xaxis.get_major_ticks():
		tick.label.set_fontsize(FS)

	#ax[0].text(.01,.75,treslab,fontsize=FS,transform=ax[0].transAxes)
	#ax[0].text(.01,.5,freslab,fontsize=FS,transform=ax[0].transAxes)
	#ax[0].text(.9,.75,burstlab,fontsize=FS,transform=ax[0].transAxes)
	ax[0].text(.015,.92,bothlab,fontsize=FS,transform=ax[0].transAxes,bbox=dict(facecolor='white',alpha=.7,edgecolor='none'),ha='left',va='top')
	ax[0].text(.99,.92,burstlab,fontsize=FS,transform=ax[0].transAxes,ha='right',va='top')

	fig.tight_layout(h_pad=-1)
	if not (args.savepan or args.PAout or args.print_stats or args.Iout):
		plt.show()

	# Circular polarization spectrum
	fig,ax = plt.subplots(2,1,sharex=True,gridspec_kw={'height_ratios':[1,3]})
	fig.subplots_adjust(hspace=0)
	# dynamic spectrum
	ax[1].imshow(dynV,aspect='auto',interpolation='none',cmap=args.cmap)
	# color in zapped channels
	zaps = np.squeeze(np.where(np.sum(dyn,axis=1)==0))
	for z in zaps:
		ax[1].axhspan(z-.5,z+.5,color=args.zapcol)
	# intensity
	ax[0].plot(V,color='b',lw=.7)
	ax[1].set_xlim([0,len(I)-1])

	ax[0].tick_params(axis='y',which='both',labelleft=False,left=False)
	ax[0].tick_params(axis='x',which='both',bottom=False)
	#ax[1].set_ylim([0.90623*nchan,0.07422*nchan])
	ax[1].set_ylim([botcut,topcut])
	ax[1].set_xlabel('Time (ms)',fontsize=FS2)
	ax[1].set_ylabel('Frequency (MHz)',fontsize=FS2)

	# change x and y ticks
	ax[1].yaxis.set_major_locator(plt.MaxNLocator(7))
	ax[1].set_xticks(nxticks)
	ax[1].set_xticklabels(nticklabs)
	ax[1].set_yticklabels(nyticks)

	for tick in ax[1].yaxis.get_major_ticks():
		tick.label.set_fontsize(FS)
	for tick in ax[1].xaxis.get_major_ticks():
		tick.label.set_fontsize(FS)

	#ax[0].text(.01,.75,treslab,fontsize=FS,transform=ax[0].transAxes)
	#ax[0].text(.01,.5,freslab,fontsize=FS,transform=ax[0].transAxes)
	#ax[0].text(.9,.75,burstlab,fontsize=FS,transform=ax[0].transAxes)
	ax[0].text(.015,.92,bothlab,fontsize=FS,transform=ax[0].transAxes,bbox=dict(facecolor='white',alpha=.9,edgecolor='none'),ha='left',va='top')
	ax[0].text(.99,.92,burstlab,fontsize=FS,transform=ax[0].transAxes,ha='right',va='top')

	fig.tight_layout(h_pad=-1)
	if not (args.savepan or args.PAout or args.print_stats or args.Iout):
		plt.show()

	# bandpass plot
	if tstplots:
		dyn2 = np.flipud(dyn)
		#print np.shape(dyn2)
		#print(dyn[500,:])
		#print np.count_nonzero(dyn[500,:])
		#print np.sum(dyn[500])
		#print np.sum(dyn[128,:])
		plt.plot(sums)
		plt.hlines(-np.std(sums),0,512,color='red')
		plt.show()

	# PA plot
	#phi = .5 * np.arctan(UU/QQ)
	#phi = np.rad2deg(phi)
	#P2 = np.square(QQ) + np.square(UU)
	#dQ = np.std(QQ)
	ddQ = np.std(QQ[np.r_[off1,off2]])
	#dU = np.std(UU)
	ddU = np.std(UU[np.r_[off1,off2]])
	print "dq {} ddq {}".format(dQ,ddQ)
	print "du {} ddu {}".format(dU,ddU)
	#dphi = np.sqrt((QQ/P2)*dQ**2 + (UU/P2)*dU**2)
	#dphi = np.rad2deg(dphi)
	ddphi = np.sqrt((np.abs(QQ)/P2)*ddQ**2 + (np.abs(UU)/P2)*ddU**2)
	ddphi = np.rad2deg(ddphi)
	stdInorm = np.std(Inorm[np.r_[off1,off2]])
	#print np.shape(phi)
	#print np.shape(dphi)
	#plt.plot(phi[500:2000],'.')
	plt.plot(phi,'.')
	plt.errorbar(np.arange(len(phi[on_phi])),phi[on_phi],yerr=dphi[on_phi],fmt='.',elinewidth=.8)
	#plt.errorbar(np.arange(len(phi[on_phi])),phi[on_phi],yerr=ddphi[on_phi],fmt='.',color='m')
	if not (args.savepan or args.PAout or args.print_stats or args.Iout):
		plt.show()

	# I plot 
	fig,ax = plt.subplots()
	Inorm = I/np.amax(I)
	stdInorm = np.std(Inorm[np.r_[off1,off2]])
	ax.plot(Inorm/stdInorm,color='m')
	ax.axhline(Ilim)
	if not (args.savepan or args.PAout or args.print_stats or args.Iout):
		plt.show()

	# another PA plot
	if tstplots:
		on_phi = np.where(Lnorm/stdLnorm > .5)
		on_phi = np.where(Lnorm > .1)
		#print on_phi
		#on_phi = np.where(Lnorm/stdLnorm > 4)
		#print np.shape(on_phi)
		#phi_on = phi[on_phi]
		#dphi_on = dphi[on_phi]
		#plt.errorbar(np.arange(len(phi)),phi,yerr=dphi,fmt='.',color='r')
		plt.errorbar(np.arange(len(phi))[on_phi],phi[on_phi],yerr=dphi[on_phi],fmt='.')
		plt.xlim([0,len(phi)])
		plt.show()

	#pipe = Popen("psrstat -j FTDp {} |grep on:start".format(A),shell=True,stdout=PIPE).stdout
	#out_start = pipe.read()
	#pipe = Popen("psrstat -j FTDp {} |grep on:end".format(A),shell=True,stdout=PIPE).stdout
	#out_end = pipe.read()

	#ee = np.array(out_end.split()[-1].split(",")).astype(int)
	#ss = np.array(out_start.split()[-1].split(",")).astype(int)

	#idx_on = np.hstack([np.arange(ss[i],ee[i],1) for i in range(len(ee))])
	#idx_on_ds = np.arange(int(idx_on[0]/32.),(int(idx_on[-1]/32.)))
	#print idx_on
	#print idx_on_ds
	#print on_phi

	#print out_start
	#print out_end
	#print ee
	#print ss

	#print "std I {}".format(np.std(I))
	#print "max I {}".format(np.amax(I))
	#II = I/np.amax(I)
	#fwhm = np.where(II > .5)
	#print "fwhm I {}".format(fwhm)
	#print "std I/Imax {}".format(np.std(II))
	#plt.plot(II)
	#plt.show()

	### statistics ###
	# add option whether to print or not
	if args.print_stats:
		I_on = np.std(I[on_phi])
		L_on = np.std(L[on_phi])
		V_on = np.std(V[on_phi])
		I_off = np.std(I[off_phi])
		L_off = np.std(L[off_phi])
		V_off = np.std(V[off_phi])

		I_sum = np.sum(I[on_phi])
		L_sum = np.sum(L[on_phi])
		V_sum = np.sum(V[on_phi])
		V_sum_abs = np.sum(abs(V[on_phi]))
		Q_sum = np.sum(abs(QQ[on_phi]))
		U_sum = np.sum(abs(UU[on_phi]))

		# sqrt(Q^2 + U^2 + V^2)/I
		pfrac = np.sum(np.sqrt(QQ[on_phi]**2 + UU[on_phi]**2 + V[on_phi]**2))/np.sum(I[on_phi])
		# sqrt(L^2 + V^2)/I
		pfrac2 = np.sum(np.sqrt(L[on_phi]**2 + V[on_phi]**2))/np.sum(I[on_phi])

		lfrac = L_sum/I_sum

		print "std V on {} off {}".format(V_on,V_off)
		print "std V on/off {}".format(V_on/V_off)
		#print "std I on {} off {}".format(I_on,I_off)
		#print "std L on {} off {}".format(L_on,L_off)
		#print "Lin pol frac {}".format(L_sum/I_sum)
		print "Cir pol frac {}".format(V_sum/I_sum)
		print "|Cir| pol frac {}".format(V_sum_abs/I_sum)

		n = os.path.splitext(os.path.basename(args.data))[0]
		nn = "{}.stats".format(n)
		nnn = os.path.join(args.print_stats,nn)

		head = "#stdV_on\tstdV_off\tstdV(on/off)\tVfrac\tabsVfrac\tLfrac\tpfrac\tpfrac2\n"
		with open(nnn,'w') as F:
			F.write(head)
			tab = "\t"
			stats = "{:.4f}{}{:.4f}{}{:.4f}{}{:.4f}{}{:.4f}{}{:.4f}{}{:.4f}{}{:.4f}\n".format(V_on,tab,V_off,tab,V_on/V_off,tab,V_sum/I_sum,tab,V_sum_abs/I_sum,tab,lfrac,tab,pfrac,tab,pfrac2)
			F.write(stats)

	if args.PAout:
		n = os.path.splitext(os.path.basename(args.data))[0]
		nn = "{}.PAinfo".format(n)
		p = os.path.join(args.PAout,nn)
		print p
		np.savetxt(p,[phi[on_phi],dphi[on_phi]])

	if args.Iout:
		n = os.path.splitext(os.path.basename(args.data))[0]
		nn = "{}.F{:.0f}_F{:.0f}_T{:.0f}_StokesI".format(n,flo,fhi,length*1e3)
		p = os.path.join(args.Iout,nn)
		print p
		print np.shape(dyn)
		np.savetxt(p,dyn)
		#np.savetxt(p,[phi[on_phi],dphi[on_phi]])

	if args.allout:
		n = os.path.splitext(os.path.basename(args.data))[0]
		nn = "{}.allStokes".format(n)
		p = os.path.join(args.allout,nn)
		#np.savetxt(p,(dyn,QQ,UU,V))
		np.savetxt(p,np.column_stack((I,QQ,UU,V)))

	if args.dynV:
		n = os.path.splitext(os.path.basename(args.data))[0]
		nn = "{}.StokesV".format(n)
		p = os.path.join(args.dynV,nn)
		np.savetxt(p,dynV)

	fig,ax = plt.subplots(2,1,sharex=True)
	ax[0].plot(np.sqrt(QQ**2+UU**2+V**2)/I)
	#print np.sqrt(np.sum(QQ[on_phi]**2 + UU[on_phi]**2 + V[on_phi]**2))/np.sum(I[on_phi])
	#print np.sqrt(np.sum(QQ[on_phi]**2 + UU[on_phi]**2 + V[on_phi]**2))
	#print np.sqrt(np.sum(QQ[on_phi]**2) + np.sum(UU[on_phi]**2) + np.sum(V[on_phi]**2))
	#print np.sum(I[on_phi])
	#print np.sum(np.sqrt(QQ[on_phi]**2 + UU[on_phi]**2 + V[on_phi]**2))
	#print np.sum(np.sqrt(L[on_phi]**2 + V[on_phi]**2))
	ax[0].set_ylim(-1,2)
	ax[1].plot(I,color='k')
	ax[1].plot(L,color='r')
	ax[1].plot(V,color='b')
	ax[0].hlines(1,xmin=0,xmax=len(I),color='k')
	if not (args.savepan or args.PAout or args.print_stats or args.Iout):
		plt.show()
	
	phi = .5 * np.arctan2(UU,QQ)
	# shift negative PAs to positive
	ophi = phi[on_phi]
	neg = np.where(ophi < 0)	
	ophi[neg] = ophi[neg] + np.pi
	ophi -= np.pi/2
	ophi_deg = np.rad2deg(ophi)
	
