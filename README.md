# `fit_all.go`
Obtains the mean direction and first two covariance directions from DOC using data from a specified number of drivers. The drivers are treated exchangeably when computing the covariance matrices for Brake==1 and Brake==0. The data for each driver are then projected onto these three directions and the results are written to /scratch, one file per driver. 
Written files are named `smproj_multi_###.txt` for each driver.
Must run `summary_trip1.R` before running this script.

## Variables used in DOC fit:

   * `OutsideTemperature`: 10 hz temperature
   * `OnStudyElapsed`: time (in days) between the start of the first trip and the current 10 hz measurement
   * `IvbssEnable`: binary indicator of whether the Ivbss sensors are turned on
   * `LaneWidth`
   * `FcwRange`: lagged range values, 3-second window
   * `Speed`: lagged speed, 3-second window, meters/second
   * `AccelPedal`: lagged acceleration pedal (in percent), 3-second window
  
There are about 100 variables total.

## Metadata/filtering

   * `data_###.txt: FcwValidTarget==1` indicator that the range values are valid
   * Total trip distance greater than zero
   * Speed at end of 3-second window is at least 7 m/s
   * Remove all rows where `OutsideTemperature` is not available

# `fit_all_pca.go`
Performs the same DOC computation as `fit_all.go` but first projects the data onto the first 10 PC directions. Requires two streams through the data (one for marginal covariance, one for projecting data and computing DOC).
Written files named `smproj_pca_###.txt` with one row per 3-second window and columns for Driver ID, mean direction score, first covariance direction score, second covariance direction score, and `y` (braking indicator)

# `hm.py`, `hm_pc.py`
Use the projected data (`smproj_` files) for each driver to compute 2d heatmaps using K-nearest neighbors regression.
The estimated conditional probability of braking is computed on a 100 by 100 grid within the range of the projected points.
Writes each computed heatmap to /scratch. Two files for each driver: mean direction/first covariance direction map, first covariance direction/second covariance direction map.
Written files prefixed `hmap_multi_meandir_cd1_` and `hmap_multi_cd1_cd2_`

# `hmap_knn.R`, `hmap_knn_pca.R`
Plot the estimated conditional braking probability heat maps, save as pdf files.

# Data files
`data_###.txt`, where \#\#\# is the driver ID, contain 10 hz measurements, including speed, range, lane width, brake indicator, and acceleration pedal
`data2_###` contain 10 hz temperature
`summary.txt` contains trip-level summary information, including the calendar start time of each trip, total trip distance, `IvbssEnable`
Files are joined using an implementation of [left joins for Dstreams](https://github.com/brookluers/dstream/blob/master/dstream/leftjoin.go).

# Todo
- [ ] Reduce memory usage of fits using all drivers. Possible culprits: `dstream.Regroup`; `LeftJoin` uses `AddCol` which copies entire column into float slice.
- [ ] Finish adding PC projections to Dstream pacakge
- [ ] `RoadtypeEvents.txt` contains start/end timestamps for windows when the road type changes. Road types are coded as 0 (unknown), 1 (freeway), 3 (major surface), 4 (minor surface), 5 (local), 6 (ramps). Need to confirm this coding. These are my guesses based on frequency tables in published Ivbss documents.
- [ ] Compute DOC directions on each driver separately. Compare spectra across drivers
