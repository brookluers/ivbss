import dask.dataframe as dd
import numpy as np
import pandas as pd

dirfile = pd.read_csv('data/directions_0pc.txt',
    #'/scratch/stats_flux/luers/directions_0pc.txt',
                      usecols=['varname'])
xnames = dirfile.varname.values
colnames = list(xnames)
colnames += ["Brake_1sec", "Driver"]
# colnames += ["Brake","meandir0","cd0","cd1"]

df = dd.read_csv("data/lagdat_small_*",
    #"/scratch/stats_flux/luers/smproj_8pc_*.txt",
                 usecols =  colnames,
                 assume_missing=True)

mddict = {}
for xn in xnames:
    mddict[xn] = ['mean','np.min','np.max','count', 'np.std','np.median']

dfa = df.groupby('Brake_1sec').agg([np.mean, np.min, np.max, np.std]) # .agg(mddict)
dfa = dfa.compute().stack().reset_index()
dfcount = df.groupby(["Driver","Brake_1sec"]).count().compute()
dfcount.to_csv("count_y.txt")
dfa.to_csv("marg_sumstats.txt")
