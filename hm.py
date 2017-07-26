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


n_neighbors = 100
d = np.loadtxt("smproj_001.txt", skiprows=1, delimiter=",")
X = d[:,0:2] # mean direction and first covariance direction
y = d[:,3] # brake indicators
print(X.shape)
print(y.shape)
print(X)
print(y)
mins = np.min(X, axis=0)
maxs = np.max(X, axis=0)
#print(mins)
#print(maxs)
#print(mins.shape)
#print(maxs.shape)
grid1 = np.linspace(mins[0], maxs[0], 100)
grid2 = np.linspace(mins[1], maxs[1], 100)
predX = cartesian((grid1, grid2))
#print(predX.shape)
#print(predX)
knn = neighbors.KNeighborsRegressor(n_neighbors, weights='uniform')
ypred = knn.fit(X, y).predict(predX)
np.savetxt("hmap_001.txt", np.column_stack((predX, ypred)), fmt="%.10f", delimiter=",",header="meandir,cd1,ypred",
           comments="")

