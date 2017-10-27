import dask.dataframe as dd
import dask.array as da
import numpy as np
import pandas as pd

dirnames = ['meandir0','cd0','cd1']
cols = ['Driver']
cols.extend(dirnames)

df = dd.read_csv('psmall_*.txt', 
                 dtype={'Driver': str},
                 usecols = cols,
                 assume_missing=True)

minmax = df.groupby('Driver').agg(['min','max']).compute()

maxID = 3

b_lwr = [] # bun upper boundaries
b_upr = [] # bin lower boundaries
h = [] # histogram counts
driver = [] #  driver ids
vname = [] #  variable names
nbins = 100
for i in range(1, maxID + 1):
    driver_padded = "{0:03d}".format(i)
    driver_nopad = "{0:d}".format(i)

    for cur_dir in dirnames:
        rg = minmax.loc[driver_nopad, cur_dir]
        ch, cb = da.histogram(df[df.Driver==driver_nopad][cur_dir].values,
                              bins = nbins, range = [rg['min'], rg['max']],
                              density=True)
        ch = ch.compute()
        h.extend(ch)
        b_lwr.extend(cb[0:(len(cb) - 1)])
        b_upr.extend(cb[1:len(cb)])
        vname.extend([cur_dir] * len(ch))
        driver.extend([i] * len(ch))

d = {"Driver": driver, "Dir": vname, "bin_lwr": b_lwr, "bin_upr": b_upr, "hdens": h}
res = pd.DataFrame(d)
res.to_csv("proj_histograms.txt", index=False)

    
    
