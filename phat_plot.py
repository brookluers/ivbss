import matplotlib.pyplot as plt
import dask.dataframe as dd
import numpy as np
import pandas as pd

df = dd.read_csv('phat_*.txt', dtype={'Driver': str}, assume_missing=True)
#print(df.dtypes)
#df.set_index('Driver')
gdf = df.groupby('Driver')

plotdir = ['meandir0', 'cd0']
maxID = 2
for cdir in plotdir:
    fig, ax = plt.subplots()
    for name in [str(i) for i in range(1,maxID+1)]:
        cg = gdf.get_group(name)
        ax.plot(cg[cdir], cg[cdir + '_phat'])
        #label=name)
        #plt.show()
        plt.title(cdir)
        plt.savefig(cdir + '_phat.png')
        #plt.close(fig)
