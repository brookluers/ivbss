import numpy as np
from sklearn import neighbors

def cartesian(arrays, out=None):
    """
    Generate a cartesian product of input arrays.

    Parameters
    ----------
    arrays : list of array-like
        1-D arrays to form the cartesian product of.
    out : ndarray
        Array to place the cartesian product in.

    Returns
    -------
    out : ndarray
        2-D array of shape (M, len(arrays)) containing cartesian products
        formed of input arrays.

    Examples
    --------
    >>> cartesian(([1, 2, 3], [4, 5], [6, 7]))
    array([[1, 4, 6],
           [1, 4, 7],
           [1, 5, 6],
           [1, 5, 7],
           [2, 4, 6],
           [2, 4, 7],
           [2, 5, 6],
           [2, 5, 7],
           [3, 4, 6],
           [3, 4, 7],
           [3, 5, 6],
           [3, 5, 7]])

    """

    arrays = [np.asarray(x) for x in arrays]
    dtype = arrays[0].dtype

    n = np.prod([x.size for x in arrays])
    if out is None:
        out = np.zeros([n, len(arrays)], dtype=dtype)

    m = n / arrays[0].size
    out[:,0] = np.repeat(arrays[0], m)
    if arrays[1:]:
        cartesian(arrays[1:], out=out[0:m,1:])
        for j in xrange(1, arrays[0].size):
            out[j*m:(j+1)*m,1:] = out[0:m,1:]
    return out

maxID = 2
n_neighbors = 1000

for curID in range(1, maxID+1):
    dname = "/scratch/stats_flux/luers/smproj_pca_" + '{:03d}'.format(curID) + ".txt"
    d = np.loadtxt(dname, skiprows=1, delimiter=",")
    print("Python: Computing heatmap using file " + dname)
    Xm = d[:,1:3] # mean direction and first covariance direction
    Xcd = d[:, 2:4] # first two covariance directions
    y = d[:,4] # brake indicators
    mins_m = np.min(Xm, axis=0)
    maxs_m = np.max(Xm, axis=0)
    mins_cd = np.min(Xcd, axis=0)
    maxs_cd = np.max(Xcd, axis=0)

    grid1_m = np.linspace(mins_m[0], maxs_m[0], 100)
    grid2_m = np.linspace(mins_m[1], maxs_m[1], 100)
    grid1_cd = np.linspace(mins_cd[0], maxs_cd[0], 100)
    grid2_cd = np.linspace(mins_cd[1], maxs_cd[1], 100)
    
    predX_m = cartesian((grid1_m, grid2_m))
    predX_cd = cartesian((grid1_cd, grid2_cd))

    knn_m = neighbors.KNeighborsRegressor(n_neighbors, weights='uniform')
    knn_cd = neighbors.KNeighborsRegressor(n_neighbors, weights='uniform')

    ypred_m = knn_m.fit(Xm, y).predict(predX_m)
    ypred_cd = knn_cd.fit(Xcd, y).predict(predX_cd)
    fname_cd = "/scratch/stats_flux/luers/hmap_pca_cd1_cd2_" + '{:03d}'.format(curID) + ".txt"
    fname_m = "/scratch/stats_flux/luers/hmap_pca_meandir_cd1_" + '{:03d}'.format(curID) + ".txt"
    np.savetxt(fname_cd, np.column_stack((predX_cd, ypred_cd)), fmt="%.10f", delimiter=",",header="cd1,cd2,ypred",
           comments="")
    np.savetxt(fname_m, np.column_stack((predX_m, ypred_m)), fmt="%.10f", delimiter=",", header="meandir,cd1,ypred",
             comments="")


