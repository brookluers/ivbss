import pandas as pd
import matplotlib
import matplotlib.patches as mpatches
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.tight_layout()

def main():
    maxID = 3
    colors = {'1': 'r', '0': 'b'}
    nrow =1
    ncol = 3
    fig, axs = plt.subplots(nrows=1, ncols=ncol,figsize=(15,10))
    min0 = -17
    max0 = 68
    min1 = -46
    max1 = 47
    axs_reshape = axs.reshape(-1)
    red_patch = mpatches.Patch(color='red', label='y=1')
    blue_patch = mpatches.Patch(color='blue', label='y=0')
    for i in range(1, maxID + 1):
        ax = axs_reshape[i-1]
        idstr = '{:03d}'.format(i)
        d = pd.read_csv('data/ldproj_0pc_small_' + idstr + '.txt',
                        usecols=['Driver', 'lda0', 'lda1','Brake_1sec'])
        ax.scatter(d['lda0'], d['lda1'], c=d['Brake_1sec'].apply(lambda x: colors[str(round(x))]),
                    marker=',', lw=0, s=1, alpha=0.5)
        ax.set_xlim([min0,max0])
        ax.set_ylim([min1,max1])
        ax.set_title("Driver " + str(i))
        if i==1:
            ax.legend(handles=[red_patch, blue_patch])
    
    fig.tight_layout(rect=[0, 0.05,1,0.95])
    fig.savefig('ldaplot_0pc.png', dpi=400)

if __name__ == "__main__":
    main()
