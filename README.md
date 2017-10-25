See <a href="http://www.brookluers.com/research/ivbss_midas_poster.pdf">here</a> for a description of the IVBSS data and difference in covariances (DOC) methodology.

# `fit_poster2.go`

Must run `summary_trip1.R` before running this script.

Convert the raw kinematic data to three-second segments with lagged 10 hz speed and range and a binary braking indicator at the end of the segment.

Project these 62-dimensional vectors onto the first *q* (e.g. 8) principal components. Then obtain the mean direction and first two covariance directions from DOC, treating the drivers exchangeably (pooling their PC-projected driving segments).

The direction loadings are saved to `/scratch` for several values of *q*. The three-second segments are projected against the PCA+DOC directions and the projected data are saved to `/scratch`. The projected data stored in files prefixed with `smproj_`, one file per driver.

## projected data (`smproj_`) file format
Comma delimited, one row per three-second driving segment.  
Variable names: `Driver`, `Trip`, `Time`, `Speed[0]`...`Speed[-30]`, `FcwRange[0]`...`FcwRange[-30]`, `pc0`, `pc1`,... `meandir0`, `cd0`, `cd1`, ...  

The scores on each principal component are stored in `pc0`, ...   
`meandir0` is the score on the DOC mean direction
`cd0` is the score on the first DOC eigenvector

## Variables in `data_###.txt`

   * `Driver`
   * `Trip`: trip number. A trip is a single ignition cycle.
   * `Time`: hundredths of a second
   * `FcwRange`: 10 hertz distance to lead vehicle (range)
   * `Speed`: 10 hz speed (meters per second)
   * `IvbssEnable`: binary indicator of whether the Ivbss warnings are active
   * `LaneWidth`
   * `AccelPedal`: acceleration pedal (percent)
   * ...

## filtering

   * `data_###.txt: FcwValidTarget==1` indicator that the range values are valid
   * Total trip distance greater than zero
   * Speed at end of 3-second window is at least 7 m/s (`Speed[0]`)

# Data files
`data_###.txt`, where \#\#\# is the driver ID, contain 10 hz measurements, including speed, range, lane width, brake indicator, and acceleration pedal
`data2_###` contain 10 hz temperature
`summary.txt` contains trip-level summary information, including the calendar start time of each trip, total trip distance, `IvbssEnable`
Files are joined using an implementation of [left joins for Dstreams](https://github.com/brookluers/dstream/blob/master/dstream/leftjoin.go).

