import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import dask.dataframe as dd
import numpy as np
import pandas as pd

df = dd.read_csv('/scratch/stats_flux/luers/phat_*.txt', dtype={'Driver': str}, assume_missing=True)
#print(df.dtypes)
#df.set_index('Driver')
gdf = df.groupby('Driver')

plotdir = ['meandir0', 'cd0', 'pc0','cd1','pc1']
plotdirlab = {'meandir0': 'Standardized difference in mean vectors',
              'cd0': 'First DOC eigenvector',
              'cd1': 'Second DOC eigenvector',
              'pc0': 'First principal component',
              'pc1': 'Second principal component'}
maxID = 108
minID = 1
plotIDs = [str(i) for i in range(minID, maxID+1)]
for cdir in plotdir:
    fig, ax = plt.subplots()
    for name in plotIDs:
        cg = gdf.get_group(name)
        ax.plot(cg[cdir], cg[cdir + '_phat'], alpha=0.6, color='#e34a33')
        #label=name)
        #plt.show()
        plt.title(plotdirlab[cdir])
        plt.xlabel('Projected value')
        plt.ylabel('Conditional braking probability')
        plt.savefig(cdir + '_phat_108.png')
        #plt.close(fig)
