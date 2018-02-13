import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

maxID = 2

for i in range(1, maxID + 1):
    istr = '{:03d}'.format(i)
    d = pd.read_csv('/scratch/stats_flux/luers/lagdat_' + idstr + '.txt',
                    usecols=['Brake_1sec', 'cd0', 'cd1', 'cd2'])
    fig = plt.figure(figsize=(4,2))
    ax = fig.add_subplot(121)
    ax2 = fig.add_subplot(122, sharex=ax)
    ax.plot(d.cd1, d.cd0, color=d.Brake_1sec, marker=',', lw=0,linestyle="")
    ax.xlabel('first cov. direction')
    ax.ylabel('mean direction')
    ax2.plot(d.cd1, d.cd2, color=d.Brake_1sec, marker=',', lw=0, linestyle="")
    ax2.xlabel('first cov. direction')
    ax2.ylabel('second cov. direction')
    fig.savefig('projscatter_' + idstr + '.png', dpi=100)

    
