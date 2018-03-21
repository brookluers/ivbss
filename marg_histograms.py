import dask.dataframe as dd
import dask.array as da
import numpy as np
import pandas as pd
#import numexpr
#numexpr.set_nthreads(1)
df = dd.read_csv("/scratch/stats_flux/luers/lagdat_*",
                 usecols = ['Speed[0]','FcwRange[0]', 'Brake_1sec'],
                 assume_missing=True)
sp_h0, sp_bins0  = da.histogram(df[df.Brake_1sec==0]['Speed[0]'].values,
                                bins=50, range=[7, 42])
sp_h1, sp_bins1  = da.histogram(df[df.Brake_1sec==1]['Speed[0]'].values,
                                bins=50, range=[7, 42])
r_h0, r_bins0 = da.histogram(df[df.Brake_1sec==0]['FcwRange[0]'].values,
                             bins=50, range=[0, 125])
r_h1, r_bins1 = da.histogram(df[df.Brake_1sec==1]['FcwRange[0]'].values,
                             bins=50, range=[0, 125])
print("Speed[0] histogram, Brake==0")
print(sp_bins0)
print(sp_h0.compute())
print("\n-------\nSpeed[0] histogram, Brake==1")
print(sp_bins1)
print(sp_h1.compute())
print("\n-------\nRange[0] histogram, Brake==0")
print(r_bins0)
print(r_h0.compute())
print("\n-------\nRange[0] histogram, Brake==1")
print(r_bins1)
print(r_h1.compute())


