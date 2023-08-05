#!/usr/bin/env python

from __future__ import print_function,division
import diffTools_func as dtf
import numpy as np
from scipy import integrate
from os import listdir
from os.path import isfile, join
import pickle
import gc

import pydmrt as pydmrt_module



class DiffTools():

    def __init__(self,mode=None,rtt=True,mfpt=False,lfpt=False,ptpx=False,cross=True,bins=False,dist=False):
        if mode is not None:
            self.decodeMode(mode)
        else:
            self.rtt=rtt
            self.mfpt=mfpt
            self.lfpt=lfpt
            self.ptpx=ptpx
            self.bins=bins
            self.cross=cross
            self.dist=dist
            self.encodeMode()
        return

    def calcTimes(self,dmrtTms,dmrtCts,dmrtVars,rtt=False):
        tms = np.array(dmrtTms)[:-1,:]
        cts = np.array(dmrtCts)
        wars = np.array(dmrtVars)
        dists = np.array(dmrtTms)[-1,:]
        if rtt:
            tms = tms.T + tms
            cts = cts.T + cts
            tms = +np.tril(tms)-np.triu(tms)
            #tms = tms*2 #make them RTT again, but this is not necessary anymore
            wars = wars.T + wars
        tms = tms/cts
        errs = np.sqrt(wars)/cts
        return dists,tms,cts,errs

    def calcPTPR(self,dmrtTms,dmrtCts):
        normal0 = np.array(dmrtTms)[:-1,0]
        normal1 = np.array(dmrtTms)[:-1,1]
        ptpx0= np.array(dmrtCts)[:,0]
        ptpx1= np.array(dmrtCts)[:,1]
        dists = np.array(dmrtTms)[-1,:]
        total0 = normal0 + ptpx0
        total1 = normal1 + ptpx1
        total = normal0 + normal1 + ptpx0 + ptpx1
        ptpx = (ptpx0 + ptpx1)/total
        return dists, ptpx, total, ptpx0/total0, ptpx1/total1

    def compute(self,data,start=-2.0, interval=0.1, end=2.0, mode=None, verb=False, radii=None, errorCalc = None):
        if mode is not None:
            self.decodeMode(mode)
        if data.shape[1]!=2:
            print("More than 1 dimension parsed. Evaluating each dim. seperately!")
            if errorCalc == "bootstrap":
                print("Error is bootstrapped.")
            dmrtTmss, dmrtCtss, dmrtUptss, dmrtVarss, dmrtDists, dmrtTPDists, dmrtErrss = [],[],[],[],[],[],[]
            for i in range(data.shape[1]-1):
                if radii is not None:
                  dmrtTms, dmrtCts, dmrtUpts, dmrtVars, dmrtDist, dmrtTPDist = pydmrt_module.dmrtInpRadii(self.mode,int(verb),data[:,[0,i+1]],radii)
                else:
                  dmrtTms, dmrtCts, dmrtUpts, dmrtVars, dmrtDist, dmrtTPDist = pydmrt_module.dmrtInp(data[:,[0,i+1]], start, interval, end,self.mode,int(verb))
                dmrtTmss.append(np.array(dmrtTms))
                dmrtCtss.append(np.array(dmrtCts))
                dmrtUptss.append(np.array(dmrtUpts))
                if errorCalc == "bootstrap":
                    dmrtErrss.append(np.array(dmrtTms)[:-1,:]/np.array(dmrtCts))
                dmrtVarss.append(np.array(dmrtVars))
                dmrtDists.append(dmrtDist)
                dmrtTPDists.append(dmrtTPDist)
            dmrtTms = np.sum(dmrtTmss,axis=0)
            dmrtTms[-1,:] = dmrtTmss[0][-1,:]
            dmrtCts = np.sum(dmrtCtss,axis=0)
            dmrtUpts = np.sum(dmrtUptss,axis=0)
            dmrtVars = np.mean(dmrtVarss,axis=0)
            if errorCalc == "bootstrap":
                dmrtErrsBootstrap = np.std(dmrtErrss,axis=0)/np.sqrt(data.shape[1]-1)
            dmrtDists = []
            dmrtTPDist = []
            #dmrtDists = [[[]] for i in tp for tp in tp2 for tp2 in dmrtDist]

        else:
            if radii is not None:
                dmrtTms, dmrtCts, dmrtUpts, dmrtVars, dmrtDist, dmrtTPDist = pydmrt_module.dmrtInpRadii(self.mode,int(verb),data,radii)
            else:
                dmrtTms, dmrtCts, dmrtUpts, dmrtVars, dmrtDist, dmrtTPDist = pydmrt_module.dmrtInp(data, start, interval, end,self.mode,int(verb))
        gc.collect()
        if self.rtt or self.mfpt or self.lfpt:
            ret1,ret2,ret3,ret4 = self.calcTimes(dmrtTms,dmrtCts,dmrtVars,rtt=self.rtt)
            if errorCalc == "bootstrap" and data.shape[1]!=2:
                ret3 = dmrtErrsBootstrap
        elif self.ptpx:
            return self.calcPTPR(dmrtTms,dmrtCts)
        if self.dist:
            return ret1,ret2,ret3,ret4,np.array(dmrtUpts),dmrtDist,dmrtTPDist
        else:
            return ret1,ret2,ret3,ret4,np.array(dmrtUpts)

    def decodeMode(self,mode):
        if mode.startswith("rt"):
            self.rtt,self.mfpt,self.lfpt,self.ptpx=True,False,False,False
        elif mode.startswith("mfpt"):
            self.rtt,self.mfpt,self.lfpt,self.ptpx=False,True,False,False
        elif mode.startswith("lfpt"):
            self.rtt,self.mfpt,self.lfpt,self.ptpx=False,False,True,False
        elif mode.startswith("ptpx"):
            self.rtt,self.mfpt,self.lfpt,self.ptpx=False,False,True,False
        if "bins" in mode:
            self.cross,self.bins=False,True
        elif "cross" in mode:
            self.cross,self.bins=True,False
        if "dist" in mode:
            self.dist=True
        else:
            self.dist=False
        self.mode=mode

    def encodeMode(self):
        if self.rtt:
            mode="rt"
        elif self.mfpt:
            mode="mfpt"
        elif self.lfpt:
            mode="lfpt"
        elif self.ptpx:
            mode="ptpx"
        if self.cross:
            mode=mode+"cross"
        elif self.bins:
            mode=mode+"bins"
        if self.dist:
            mode=mode+"dist"
        self.mode=mode

    # DEPRECATED
    def pydmrt(self, inputVec, start=-2.0, interval=0.1, end=2.0, mode="rtcross", verb=False,radii=None):
        print("WARNING: this wrapper is deprecated! Use DiffTools.compute(args) instead.")
        if radii is not None:
            dmrtTms, dmrtCts, dmrtUpts, dmrtVars, dmrtDist, dmrtTPDist = pydmrt_module.dmrtInpRadii(self.mode,int(verb),inputVec,radii)
        else:
            dmrtTms, dmrtCts, dmrtUpts, dmrtVars, dmrtDist, dmrtTPDist = pydmrt_module.dmrtInp(inputVec, start, interval, end,self.mode,int(verb))
        return dmrtTms, dmrtCts, dmrtUpts, dmrtDist, dmrtTPDist




def loadEvalTxt(folder,filestring,mode,optionalstring="",verbose=0, recompute = False, rtCorrect=True):
	finalTms=[]
	finalCts=[]
	if checkForPkl(folder,optionalstring+"_tms")==True and recompute==False:
		tms =   openTemp(folder+filestring+mode+"_"+optionalstring+"_tms.pkl")
		cts =   openTemp(folder+filestring+mode+"_"+optionalstring+"_cts.pkl")
		dists = openTemp(folder+filestring+mode+"_"+optionalstring+"_dists.pkl")
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
	safeTemp(tms,folder+filestring+mode+"_"+optionalstring+"_tms.pkl")
	safeTemp(cts,folder+filestring+mode+"_"+optionalstring+"_cts.pkl")
	safeTemp(dists,folder+filestring+mode+"_"+optionalstring+"_dists.pkl")
	return tms,cts,dists


def getSolRdfFromTxt(folder,filestring,optionalstring, mfpt="mfpt",cross="cross",verbose=0, recompute=False):
		finalRdf=[]
		if checkForPkl(folder,optionalstring+"_hist")==True and recompute==False:
			rdf =   openTemp(folder+filestring+"_"+optionalstring+"_hist.pkl")
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
		safeTemp(rdf,folder+filestring+"_"+optionalstring+"_hist.pkl")
		return rdf


def fullAnalysis(dmrtMat, countsMat,dists, rdf, rmin = -1.5, rmax= 1.5,minVal = 0.2, maxVal = 1.0,mfpt="mfpt",cross="cross", dim=1, verbose=0, smoothwidth=0):
    downsampleD=1
    print("Resolution: ", countsMat.shape,"Max counts: ", countsMat.max())
    dmrtMat,dists = dtf.filterDmrtMat(dmrtMat,dists ,rmin,rmax)
    dmrtMat=dmrtMat[::downsampleD,::downsampleD]
    dists=dists[::downsampleD]
    if np.isnan(dmrtMat[:,-1]).any():
        print("Warning some MFPTs did not reach final Qf after filtering!")
    if mfpt=="mfpt":
        dmrtMat = -dmrtMat
    if smoothwidth==0:
        ddmrt,ddists = dtf.diffDmrtMat(dmrtMat,dists)
    else:
        ddmrt,ddists = dtf.diffDmrtMatYann(dmrtMat,dists,smoothwidth)
    ddmrt = ddmrt[:,:-1]
    np.fill_diagonal(ddmrt,0.0)
    ddmrt[ddmrt == 0] = 2.22044604925e-16 # epsilon

    rdf[:,1]=rdf[:,1]/np.sum(rdf[:,1])

    if dim ==1:
        f = -np.log(rdf[:,1])
    elif dim ==3:
        -np.log(r[:, 1])-2*np.log(r[:, 0])
    print(np.exp(-f).shape,rdf[0].shape)
    integral = integrate.cumtrapz(np.exp(-f), rdf[:,0], initial=0.0)
    if mfpt=="mfpt":
        d = np.exp(f) * integral
    else:
        d = np.exp(f) * integral[-1]
    interp = np.interp(ddists, rdf[:,0], d)

    diff = np.einsum('ij,i->ij', ddmrt, 1/interp)
    invdiff = 1/diff

    filtMat = dtf.filterMatrix(minVal,maxVal,dists)

    finalDmrts = dtf.averageSelFilter(diff,filtMat[1:,1:])

    final =1/finalDmrts
    print(finalDmrts)

    return final,ddists, ddmrt, diff

