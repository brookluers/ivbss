import dask.dataframe as dd
import numpy as np
import pandas as pd

dirfile = pd.read_csv('/scratch/stats_flux/luers/directions_0pc.txt',
    #'data/directions_0pc.txt',
                      usecols=['varname'])
xnames = dirfile.varname.values
colnames = list(xnames)
colnames += ["Brake_1sec", "Driver"]

df = dd.read_csv("/scratch/stats_flux/luers/lagdat_0pc_*.txt",
    #"data/lagdat_small_*",
                 usecols =  colnames,
                 assume_missing=True)

dfa = df.groupby(['Driver', 'Brake_1sec']).agg([np.mean])
dfa = dfa.compute().stack().reset_index()
dfa.to_csv("/scratch/stats_flux/luers/driver_means.txt")
