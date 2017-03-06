#!/usr/bin/python

import pydmrt
import diffTools as dt
import numpy as np
import matplotlib.pyplot as plt
from pylab import *

print("Starting pydmrt test")

testtraj = np.loadtxt("./test_traj_doublewell_od.txt")

#&inputVec1, &start, &interval,&end,&mode,&verb
dmrtTms, dmrtCts, dmrtUpts = pydmrt.dmrtInp(testtraj, -1.7, 0.1, 1.7, "rtcross",True)

refTms = np.loadtxt("./reference_Tms.txt")
refCts = np.loadtxt("./reference_Cts.txt")
refUpts = np.loadtxt("./reference_Upts.txt")

diffTms = np.sum(np.sum(dmrtTms-refTms))
diffCts = np.sum(np.sum(dmrtCts-refCts))
diffUpts = np.sum(np.sum(dmrtUpts-refUpts))
if (diffTms + float(diffCts) + float(diffUpts) != 0.0):
  raise("WARNING: Differences to reference tables: tms %f, cts %i, upts %i "%(diffTms, diffCts, diffUpts))

tms = np.array(dmrtTms)[:-1,:]
cts = np.array(dmrtCts)
dists = np.array(dmrtTms)[-1,:]

tms = tms/cts
tms = tms.T + tms
tms = +np.tril(tms)-np.triu(tms)


for i in range(0,tms.shape[0]-1):
  plt.plot(dists,tms[:,i])
title("RT times of overdamped doublewell Langevin sim.")
plt.show()


hist,bins = np.histogram(testtraj[:,1], bins=dists)
rdf = np.zeros((hist.shape[0],2))
rdf[:,0] = bins[:-1]
rdf[:,1] = hist


final,ddists, ddmrt, diff = dt.fullAnalysis(tms, cts, dists, rdf, rmin = -1.5 , rmax= 1.5 ,minVal = -1., maxVal = 1.,
                mfpt="rt",cross="cross", dim=1, verbose=True, smoothwidth=0.0)

plt.plot(ddists,final)
title("diffusivity profile of overdamped Langevin simulation")
plt.show()
