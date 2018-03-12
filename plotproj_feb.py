import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.tight_layout(w_pad=0.75)

def main():
    maxID = 3
    
    for i in range(1, maxID + 1):
        idstr = '{:03d}'.format(i)
        d = pd.read_csv('data/projdat_0pc_' + idstr + '.txt',
        # '/scratch/stats_flux/luers/lagdat_' + idstr + '.txt',
                    usecols=['Brake_1sec', 'doc0', 'doc1', 'doc2'])
        d0 = d.query('Brake_1sec==0')
        d1 = d.query("Brake_1sec==1")
        fig = plt.figure()
        # fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2)
        fig, ((ax1, ax2)) = plt.subplots(nrows=1, ncols=2, figsize=(8.5,4))
        fig.suptitle("Driver " + str(i))
        ax1.scatter(d0.doc0, d0.doc1, color='blue', s=1,marker=',',alpha=0.5)
        ax1.scatter(d1.doc0, d1.doc1, color='red',  s=1, marker=',',alpha=0.5)
        ax1.set_ylabel('first cov. direction')
        ax1.set_xlabel('mean direction')
        ax2.scatter(d0.doc1, d0.doc2, color='blue', s=1, marker=',',alpha=0.5)
        ax2.scatter(d1.doc1, d1.doc2, color='red', s=1,marker=',',alpha=0.5)
        ax2.set_xlabel('first cov. direction')
        ax2.set_ylabel('second cov. direction')
        # ax3.hist(d.
        fig.savefig('projscatter_' + idstr + '.png', dpi=250)

if __name__ == "__main__":
    main()
