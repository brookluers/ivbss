import dask.dataframe as dd
import dask.array as da
import numpy as np
import pandas as pd

dirnames = ['meandir','cd1','cd2']
cols = ['Driver','Brake']
cols.extend(dirnames)

df = dd.read_csv('/scratch/stats_flux/luers/scores_8pc_py_*.txt', 
                 dtype={'Driver': str},
                 usecols = cols,
                 assume_missing=True)

minmax = df.groupby('Driver').agg(['min','max']).compute()

maxID = 108

b_lwr = [] # bun upper boundaries
b_upr = [] # bin lower boundaries
h = [] # histogram counts
driver = [] #  driver ids
vname = [] #  variable names
brake_ind = []
nbins = 100
for i in range(1, maxID + 1):
    driver_padded = "{0:03d}".format(i)
    driver_nopad = "{0:d}".format(i)

    for cur_dir in dirnames:
        rg = minmax.loc[driver_nopad, cur_dir]
        ch, cb = da.histogram(df[df.Driver==driver_nopad][cur_dir].values,
                              bins = nbins, range = [rg['min'], rg['max']],
                              density=True)
        ch1, cb1 = da.histogram(df[ ( df.Driver==driver_nopad) & (df.Brake > 0)][cur_dir].values,
                                bins = nbins, range = [rg['min'], rg['max']],
                                density=True)
        ch0, cb0 = da.histogram(df[ ( df.Driver==driver_nopad) & (df.Brake < 0.5 )][cur_dir].values,
                                bins = nbins, range = [rg['min'], rg['max']],
                                density=True)
        ch1 = ch1.compute()
        ch0 = ch0.compute()
        ch = ch.compute()
        h.extend(ch)
        h.extend(ch1)
        h.extend(ch0)
        b_lwr.extend(cb[0:(len(cb) - 1)])
        b_lwr.extend(cb1[0:(len(cb1) - 1)])
        b_lwr.extend(cb0[0:(len(cb0) - 1)])
        b_upr.extend(cb[1:len(cb)])
        b_upr.extend(cb1[1:len(cb1)])
        b_upr.extend(cb0[1:len(cb0)])
        vname.extend([cur_dir] * (len(ch) + len(ch1) + len(ch0)))
        driver.extend([i] * (len(ch) + len(ch1) + len(ch0)))
        brake_ind.extend([-99] * len(ch))
        brake_ind.extend([1.0] * len(ch1))
        brake_ind.extend([0.0] * len(ch0))
        

d = {"Driver": driver, "Dir": vname, "bin_lwr": b_lwr, "bin_upr": b_upr, "hdens": h, "Brake": brake_ind}
res = pd.DataFrame(d)
res.to_csv("/scratch/stats_flux/luers/proj_histograms.txt", index=False)

    
    
