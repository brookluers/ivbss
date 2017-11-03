import numpy as np
import pandas as pd

nlag = 30
xnames = []
for j in range(0, nlag + 1):
    xnames.append("Speed[" + str(-j) + "]")
for j in range(0, nlag + 1):
    xnames.append("FcwRange[" + str(-j) + "]")
colnames = list(xnames)
colnames += ["Brake"]

npc = 8

maxid = 108

for cid in range(1, maxid + 1):
    id_pad = '{:03d}'.format(cid)
    dir = pd.read_csv('data/directions_' + str(npc) +  'pc.txt',
                      index_col=0,
                      header=0)
    dcur = pd.read_csv('data/smproj_8pc_' + id_pad + '.txt',
                       header=0,   
                       usecols=colnames, 
                       engine='c')
# Inner product of kinematics and direction loadings, add brake column
    dcur = dcur[xnames].dot(dir.loc[xnames]).assign(Brake=dcur.Brake)
    dcur.to_csv('data/scores_' + str(npc) + 'pc_py_' + id_pad + '.txt',
                index = False)
