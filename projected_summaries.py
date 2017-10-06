import dask.dataframe as dd
import numpy as np
import pandas as pd

nlag = 30
xnames = []
for j in range(0, nlag + 1):
    xnames.append("Speed[" + str(-j) + "]")
for j in range(0, nlag + 1):
    xnames.append("FcwRange[" + str(-j) + "]")
colnames = list(xnames)
colnames += ["Brake","meandir0","cd0","cd1"]

df = dd.read_csv("/scratch/stats_flux/luers/smproj_8pc_*.txt",
                 usecols =  colnames,
                 assume_missing=True)

mdbins = [0.75, 1.25, 1.75, 2.25, 2.75, 3.25]
cd0bins = [-5.0, -2.5, -1, -0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1.0,1.2,1.4]
print("bins for mean direction: ")
print(mdbins)
print("bins for cd0: ")
print(cd0bins)

df = df.assign(md_interval = lambda x: np.digitize(x.meandir0, mdbins))
df = df.assign(cd0_interval = lambda x: np.digitize(x.cd0, cd0bins))
mddict = {}
for xn in xnames:
    mddict[xn] = 'mean'
mddict['meandir0'] = 'count'
cd0dict = {}
for xn in xnames:
    cd0dict[xn] = 'mean'
cd0dict['cd0'] = 'count'

dfmd = df.groupby('md_interval').agg(mddict)
dfcd0 = df.groupby('cd0_interval').agg(cd0dict)
dfmd = dfmd.compute().rename(columns = {"meandir0": "nsegments"})
dfcd0 = dfcd0.compute().rename(columns = {"cd0": "nsegments"})

dfmd.to_csv("md0_binmeans.txt")
dfcd0.to_csv("cd0_binmeans.txt")
