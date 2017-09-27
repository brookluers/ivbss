import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import dask.dataframe as dd
import numpy as np
import pandas as pd

df = dd.read_csv('/scratch/stats_flux/luers/smproj_multi_*.txt', 
                 dtype={'Driver': str}, 
                 usecols = ['Driver', 'Brake'],
                 assume_missing=True)
ycounts = df.groupby('Driver')['Brake'].aggregate(['sum','count'])
ycounts = ycounts.compute()
ycounts.to_csv('ycounts.txt')
