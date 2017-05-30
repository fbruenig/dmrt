#!/usr/bin/env python

from __future__ import print_function,division
import diffTools_func as dtf
import unitools as uni
import numpy as np
from scipy import integrate
from os import listdir
from os.path import isfile, join

import pydmrt as pydmrt_module

def pydmrt(inputVec, start=-2.0, interval=0.1, end=2.0, mode="rtcross", verb=False):
    dmrtTms, dmrtCts, dmrtUpts = pydmrt_module.dmrtInp(inputVec, start, interval, end, mode,verb)
    return dmrtTms, dmrtCts, dmrtUpts

def loadEvalTxt(folder,filestring,mode,optionalstring="",verbose=0, recompute = False, rtCorrect=True):
	finalTms=[]
	finalCts=[]
	if uni.checkForPkl(folder,optionalstring+"_tms")==True and recompute==False:
		tms =   uni.openTemp(folder+filestring+mode+"_"+optionalstring+"_tms.pkl")
		cts =   uni.openTemp(folder+filestring+mode+"_"+optionalstring+"_cts.pkl")
		dists = uni.openTemp(folder+filestring+mode+"_"+optionalstring+"_dists.pkl")
		return tms,cts,dists
	for f in listdir(folder):
		if f.startswith(filestring) and mode in f and optionalstring in f and "cts" in f and not ".pkl" in f:
			if verbose:
				print("Loading", f)
			ar = np.genfromtxt(folder+f)
			finalCts.append(ar)
		if f.startswith(filestring) and mode in f and optionalstring in f and "tms" in f and not ".pkl" in f:
			ar = np.genfromtxt(folder+f)
			finalTms.append(ar)
	cts = finalCts[0]
	for ar in finalCts[1:]:
		cts=np.add(cts,ar)
	tms = finalTms[0][:-1,:]
	dists = finalTms[0][-1,:]
	for ar in finalTms[1:]:
		tms=np.add(tms,ar[:-1,:])
	if not mode.startswith("tftp"):
		tms = tms/cts
	if mode.startswith("tftp"):
		dists=dists-0.5*(dists[1]-dists[0])
	if "rt" in mode and rtCorrect:
		tms = tms.T + tms
		tms = +np.tril(tms)-np.triu(tms)
	uni.safeTemp(tms,folder+filestring+mode+"_"+optionalstring+"_tms.pkl")
	uni.safeTemp(cts,folder+filestring+mode+"_"+optionalstring+"_cts.pkl")
	uni.safeTemp(dists,folder+filestring+mode+"_"+optionalstring+"_dists.pkl")
	return tms,cts,dists


def getSolRdfFromTxt(folder,filestring,optionalstring, mfpt="mftp",cross="cross",verbose=0, recompute=False):
		finalRdf=[]
		if uni.checkForPkl(folder,optionalstring+"_hist")==True and recompute==False:
			rdf =   uni.openTemp(folder+filestring+"_"+optionalstring+"_hist.pkl")
			return rdf

		for f in listdir(folder):
			if f.startswith(filestring) and optionalstring in f and "hist" in f and not ".pkl" in f:
				if verbose:
					print("Loading", f)
				ar = np.genfromtxt(folder+f)
				finalRdf.append(ar)
		rdf=finalRdf[0]
		for r in finalRdf[1:]:
			rdf[:,1]=np.add(rdf[:,1],r[:,1])
		if cross=="cross":
			rdf[:-1,0]=(rdf[1:,0]+rdf[:-1,0])/2
		uni.safeTemp(rdf,folder+filestring+"_"+optionalstring+"_hist.pkl")
		return rdf


def fullAnalysis(dmrtMat, countsMat,dists, rdf, rmin = -1.5, rmax= 1.5,minVal = 0.2, maxVal = 1.0,mfpt="mftp",cross="cross", dim=1, verbose=0, smoothwidth=0):
    downsampleD=1
    print("Resolution: ", countsMat.shape,"Max counts: ", countsMat.max())
    dmrtMat,dists = dtf.filterDmrtMat(dmrtMat,dists ,rmin,rmax)
    dmrtMat=dmrtMat[::downsampleD,::downsampleD]
    dists=dists[::downsampleD]
    if np.isnan(dmrtMat[:,-1]).any():
        print("Warning some MFPTs did not reach final Qf after filtering!")
    if mfpt=="mftp":
        dmrtMat = -dmrtMat
    if smoothwidth==0:
        ddmrt,ddists = dtf.diffDmrtMat(dmrtMat,dists)
    else:
        ddmrt,ddists = dtf.diffDmrtMatYann(dmrtMat,dists,smoothwidth)
    ddmrt = ddmrt[:,:-1]
    np.fill_diagonal(ddmrt,0.0)
    ddmrt[ddmrt == 0] = uni.epsilon

    rdf[:,1]=rdf[:,1]/np.sum(rdf[:,1])

    if dim ==1:
        f = -np.log(rdf[:,1])
    elif dim ==3:
        -np.log(r[:, 1])-2*np.log(r[:, 0])
    print(np.exp(-f).shape,rdf[0].shape)
    integral = integrate.cumtrapz(np.exp(-f), rdf[:,0], initial=0.0)
    if mfpt=="mftp":
        d = np.exp(f) * integral
    else:
        d = np.exp(f) * integral[-1]
    interp = np.interp(ddists, rdf[:,0], d)

    diff = np.einsum('ij,i->ij', ddmrt, 1/interp)
    invdiff = 1/diff

    filtMat = dtf.filterMatrix(minVal,maxVal,dists)

    finalDmrts = dtf.averageSelFilter(diff,filtMat[1:,1:])

    final =1/finalDmrts

    return final,ddists, ddmrt, diff
