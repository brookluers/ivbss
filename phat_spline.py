import pandas as pd
from scipy import interpolate
import numpy as np

rr = pd.read_csv('proj_ranges.txt', header=[0,1])

# collapse multiindex
rr.columns = ['_'.join(c) for c in rr.columns.values]

md0min = rr.meandir0_min.max() 
md0max = rr.meandir0_max.min()
cdmin = rr.cd0_min.max()
cdmax= rr.cd0_max.min()

mdgrid, cdgrid = np.mgrid[md0min:md0max:100j, cdmin:cdmax:100j]

maxID = 2

for i in range(1, maxID+1):
    idstr = '{:03d}'.format(i)
    d = pd.read_csv('smproj_8pc_' + idstr + '.txt',
                    usecols=['Brake','meandir0','cd0'])
    uvs_md = interpolate.UnivariateSpline(np.array(d.meandir0), np.array(d.Brake))
    uvs_cd = interpolate.UnivariateSpline(np.array(d.cd0), np.array(d.Brake))
    #bvs = interpolate.SmoothBivariateSpline(np.array(d.meandir0), np.array(d.cd0), np.array(d.Brake))
    tck = interpolate.bisplrep(np.array(d.meandir0), np.array(d.cd0), np.array(d.Brake))
    yhat_bvs = interpolate.bisplev(mdgrid[:,0], cdgrid[0,:], tck)
    #yhat_bvs = bvs(mdgrid[:,0], cdgrid[0,:])
    yhat_md_uvs=  uvs_md(mdgrid[:,0])
    yhat_cd_uvs = uvs_cd(cdgrid[0,:])

    out = np.vstack((np.ravel(yhat_bvs), np.ravel(mdgrid), np.ravel(cdgrid))).T
    out_1d = np.vstack((mdgrid[:,0], cdgrid[0,:], yhat_md_uvs, yhat_cd_uvs)).T
    out = np.hstack((out, np.reshape([i] * out.shape[0], (out.shape[0], 1))))
    out_1d = np.hstack((out_1d, np.reshape([i] * out_1d.shape[0], (out_1d.shape[0], 1))))

    np.savetxt('hmap_spl_' + idstr + '.txt', out, delimiter=',',
               header='yhat,md0,cd0,Driver',
               comments='')
    np.savetxt('phat1d_spl_' + idstr + '.txt', out_1d, delimiter=',',
               comments='',
               header='md0,cd0,yhat_md0,yhat_cd0,Driver')



