#!/usr/bin/env python

from __future__ import print_function,division

import numpy as np

epsilon = np.finfo(float).eps

def averageSelFilter(dmrts,filt):
    avg = []
    for i,d in enumerate(dmrts):
        avg.append(np.nanmean(d[filt[i,:]==1]))
    return np.array(avg)

def diffDmrtMat(dmrtMat,dists):
		ddmrt = np.einsum('ij,i->ij', np.diff(dmrtMat,axis=0), 1.0/np.diff(dists))
		dr = dists[1] - dists[0]
		dists = dists[:-1] + dr * 0.5
		return ddmrt,dists

def diffDmrtMatYann(dmrtMat,dists,width):
		ddmrt = np.zeros((dmrtMat.shape[0]-1,dmrtMat.shape[1]))
		for j,val in enumerate(dists[:-1]):
			sel = np.logical_and(dists<val+width/2, dists>val-width/2)
			for i,d in enumerate(dmrtMat.T[:-1,:]):
				sell = np.logical_and(np.logical_and(sel, ~np.isnan(d)),d!=0.0)
				if sum(sell==True)>2:
					ist = dists[sell]
					dd = d[sell]
					ddmrt[j,i], c = np.linalg.lstsq(np.vstack([ist, np.ones(len(ist))]).T, dd)[0]
		dists = dists[:-1]
		return ddmrt,dists

def filterDmrtMat(dmrtMat,dists,minVal,maxVal):
		if (dists[0]> minVal):
				print("Warning: data range is samller than considered minR")
		if (dists[-1]< maxVal):
				print("Warning: data range is smaller than considered maxR")
		newdists = dists[dists>=minVal]
		newdists = newdists[newdists<=maxVal]
		newDmrt = dmrtMat[np.where(dists==newdists[0])[0][0]:np.where(dists==newdists[-1])[0][0]+1,np.where(dists==newdists[0])[0][0]:np.where(dists==newdists[-1])[0][0]+1]
		return newDmrt, newdists

def filterMatrix(minVal,maxVal,dists):
    matDistToEscape = np.zeros((dists.shape[0],dists.shape[0]))
    for i in range(dists.shape[0]):
        matDistToEscape[i,:]=abs(dists-dists[i])
    matDistToEscape[matDistToEscape<=minVal]=0.0
    matDistToEscape[matDistToEscape>=maxVal]=0.0
    matDistToEscape[matDistToEscape!=0.0]=1.0
    return matDistToEscape

