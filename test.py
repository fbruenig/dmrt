#!/usr/bin/python

from diffTools import *
import numpy as np
import matplotlib.pyplot as plt

def plotTimes(dists,tms):
    for i in range(0,tms.shape[0]-1):
      plt.plot(dists,tms[:,i])
    plt.title("RT times of overdamped doublewell Langevin sim.")
    plt.show()

print("Starting pydmrt test")

testtraj = np.loadtxt("./test_traj_doublewell_od.txt")
#testtraj = np.loadtxt("./test_traj_short.txt")

dt=DiffTools()

#&inputVec1, &start, &interval,&end,&mode,&verb
dmrtTms, dmrtCts, dmrtUpts, dmrtDist, dmrtTPDist = dt.pydmrt(testtraj, -1.7, 0.1, 1.7, "rtcross",True)

refTms = np.loadtxt("./reference_Tms.txt")
refCts = np.loadtxt("./reference_Cts.txt")
refUpts = np.loadtxt("./reference_Upts.txt")

diffTms = np.sum(np.sum(dmrtTms-refTms))
diffCts = np.sum(np.sum(dmrtCts-refCts))
diffUpts = np.sum(np.sum(dmrtUpts-refUpts))
if (diffTms + float(diffCts) + float(diffUpts) != 0.0):
  print("WARNING: Differences to reference tables: tms %f, cts %i, upts %i "%(diffTms, diffCts, diffUpts))
else:
  print("Reference data is reproduced!")

dists,tms,cts,errs = dt.calcTimes(dmrtTms,dmrtCts,True)
plotTimes(dists,tms)

hist,bins = np.histogram(testtraj[:,1],bins=dists)
rdf=np.zeros((hist.shape[0],2))
rdf[:,0]=bins[:-1]
rdf[:,1]=hist

final,ddists, ddmrt, diff = fullAnalysis(tms, cts, dists, rdf, rmin = -1.5 , rmax= 1.5 ,minVal = -1., maxVal = 1.,
                mfpt="rt",cross="cross", dim=1, verbose=True, smoothwidth=0.0)

plt.plot(ddists,final)
plt.plot(ddists, np.ones(ddists.shape),'--')
plt.title("diffusivity profile of overdamped Langevin simulation")
plt.show()


print("Short test: radii mode")
radii = [-2,-1.5,0.0,1.5,2]
dists,tms,cts,upts,vars,rtDist,rtTPDist = dt.compute(testtraj,radii=radii, mode="rtcrossdist",verb=True)
print(tms)
hist,bins = np.histogram(np.concatenate([np.array(rtDist[-2][1]),np.array(rtDist[1][-2])]))
plt.plot(bins[:-1],hist)
plt.title("MFPT distribution of doublewell barrier hopping.")
plt.show()


print("Short test: radii and bins mode")
dists,tms,cts,upts,vars = dt.compute(testtraj,radii=radii, mode="rtbins",verb=True)
print(tms)
