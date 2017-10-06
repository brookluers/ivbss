package main

import (
	"fmt"
	"math"
	"github.com/gonum/matrix/mat64"
	"bufio"
	"strconv"
	"github.com/brookluers/dstream/dstream"
	"os"
	//"sort"
	"time"
)

const (
	maxSpeedLag int     = 30 //30 samples = 30 * 100 milliseconds = 3 seconds
	maxRangeLag int     = 30
	minCount    int     = 100
	nrow        int     = 100
	ncol        int     = 100
	minSpeed    float64 = 7.0
	chunkSize   int     = 10000
)

// eliminates rows where the variable is less than 0
// or is NaN
func selectGtNotNaN(w float64) dstream.FilterFunc {
	f := func(x interface{}, ma []bool) bool {
		anydrop := false
		z := x.([]float64)
		for i, v := range z {
			if v < w {
				ma[i] = false
				anydrop = true
			} else if math.IsNaN(v) {
				ma[i] = false
				anydrop = true
			}
		}
		return anydrop
	}
	return f
}

// remove rows with math.NaN values
func notNaN(x interface{}, ma []bool) bool {
	anydrop := true
	z := x.([]float64)
	for i, v := range z {
		if math.IsNaN(v) {
			ma[i] = false
			anydrop = true
		}
	}
	return anydrop
}

// selectGt returns a funvtion that can be used with Filter to retain
// only rows where a given variable is greater than a provided value
func selectGt(w float64) dstream.FilterFunc {
	f := func(x interface{}, ma []bool) bool {
		anydrop := false
		z := x.([]float64)
		for i, v := range z {
			if v < w {
				ma[i] = false
				anydrop = true
			}
		}
		return anydrop
	}
	return f
}

func selectBt(lwr, upr float64) dstream.FilterFunc {
	f := func(x interface{}, ma []bool) bool {
		anydrop := true
		z := x.([]float64)
		for i, v := range z {
			if v < lwr || v > upr {
				ma[i] = false
				anydrop = true
			}
		}
		return anydrop
	}
	return f
}

//selectEq creates a FilterFunc that will drop any
// rows where the variable is not equal to w
func selectEq(w float64) dstream.FilterFunc {
	f := func(x interface{}, ma []bool) bool {
		anydrop := true
		z := x.([]float64)
		for i, v := range z {
			if v != w {
				ma[i] = false
				anydrop = true
			}
		}
		return anydrop
	}
	return f
}

//pminFunc returns an ApplyFunc that
// computes the minimum "horizontally" across
// a set of variables defined by vnames (variable names)
func pminFunc(vnames []string) dstream.ApplyFunc {
	nvar := len(vnames)
	data := make([][]float64, nvar)
	f := func(v map[string]interface{}, z interface{}) {
		for i, vn := range vnames {
			data[i] = v[vn].([]float64)
		}
		y := z.([]float64)
		for i := 0; i < len(y); i++ {
			y[i] = data[0][i]
			for j := 1; j < nvar; j++ {
				if data[j][i] < y[i] {
					y[i] = data[j][i]
				}
			}
		}

	}
	return f
}

// diffCols returns an ApplyFunc that computes the
// difference between column c1 and column c2
// c1 - c2
func diffCols(c1, c2 string) dstream.ApplyFunc {
	f := func(v map[string]interface{}, z interface{}) {
		v1 := v[c1].([]float64)
		v2 := v[c2].([]float64)
		res := z.([]float64)
		for i := 0; i < len(res); i++ {
			res[i] = v1[i] - v2[i]
		}
	}
	return f
}

// sumCols returns an ApplyFunc that computes the
// sum of column c1 and column c2
func sumCols(c1, c2 string) dstream.ApplyFunc {
	f := func(v map[string]interface{}, z interface{}) {
		v1 := v[c1].([]float64)
		v2 := v[c2].([]float64)
		res := z.([]float64)
		for i := 0; i < len(res); i++ {
			res[i] = v1[i] + v2[i]
		}
	}
	return f
}

// applyIdent is an ApplyFunc that leaves a column
// unchanged (use for renaming columns)
func applyIdent(vname string) dstream.ApplyFunc {
	f := func(v map[string]interface{}, z interface{}) {
		vdat := v[vname].([]float64)
		y := z.([]float64)
		for i := range vdat {
			y[i] = vdat[i]
		}
	}
	return f
}

// driverTripTimeId is an ApplyFunc that creates
//  an ID value for each 10 hz measurement for each driver-trip
func driverTripTimeId(v map[string]interface{}, z interface{}) {
	dr := v["Driver"].([]float64)
	tr := v["Trip"].([]float64)
	ti := v["Time"].([]float64)
	//var dri, tri, tii uint64
	res := z.([]float64)
	for i := range dr {
		//dri = uint64(dr[i])
		//tri = uint64(tr[i])
		//tii = uint64(ti[i])
		res[i] = dr[i]*1000000000 + tr[i]*1000000 + ti[i]
	}
}

//driverTripId is an ApplyFunc that creates
// an ID value for each driver-trip combination
func driverTripId(v map[string]interface{}, z interface{}) {
	dr := v["Driver"].([]float64)
	tr := v["Trip"].([]float64)
	res := z.([]float64)
	for i := range dr {
		res[i] = dr[i]*1000000000 + tr[i]*1000000
	}
}

//fbrake populates z.([]float64) with zeroes
// at each position corresponding to either no braking
// or the first row where Brake==1
// all other positions of z are populated with ones
func fbrake(v map[string]interface{}, z interface{}) {

	b := v["Brake"].([]float64)
	y := z.([]float64)

	y[0] = 0
	for i := 1; i < len(y); i++ {
		y[i] = 0
		if (b[i] == 1) && (b[i-1] == 1) {
			y[i] = 1
		} else if b[i] == 1 {
			y[i] = 0
		}
	}
}

func diagMat(v []float64) *mat64.Dense {
	n := len(v)
	m := mat64.NewDense(n, n, nil)
	for i := 0; i < n; i++ {
		m.Set(i, i, v[i])
	}
	return m
}

func main() {
	maxDriverID := 108
	fnames := make([]string, maxDriverID)
	fmt.Println("File names:")
	for i := 1; i <= maxDriverID; i++ {
		fnames[i-1] = fmt.Sprintf("/nfs/turbo/ivbss/LvFot/data_%03d.txt", i)
		fmt.Println(fnames[i-1])
	}

	// reader for main data files
	mrdr := dstream.NewMultiReadSeek0(fnames, true)

	ivb := dstream.FromCSV(mrdr).SetFloatVars([]string{"Driver", "Trip",
		"Time", "Speed", "Brake",
		"FcwValidTarget", "FcwRange"}).HasHeader().Done()

	// ID column for driver-trip
	ivb = dstream.Apply(ivb, "DriverTrip", driverTripId, "float64")
	ivb = dstream.Convert(ivb, "DriverTrip", "uint64")
	ivb = dstream.Apply(ivb, "DriverTripTime", driverTripTimeId, "float64")
	// ID column for 10-hz measurement within each driver-trip
	ivb = dstream.Convert(ivb, "Driver", "uint64")

	srdr, err := os.Open("/nfs/turbo/ivbss/LvFot/summary.txt")
	if err != nil {
		panic(err)
	}
	summary := dstream.FromCSV(srdr).SetFloatVars([]string{"Driver", "Trip", "StartTime", "EndTime", "Distance", "TODTripStart"}).SetChunkSize(5000).HasHeader().Done()

	summary = dstream.Apply(summary, "DriverTrip", driverTripId, "float64")
	summary = dstream.Convert(summary, "DriverTrip", "uint64")
	summary = dstream.Segment(summary, []string{"DriverTrip"})
	summary = dstream.Apply(summary, "SummaryDistance", applyIdent("Distance"), "float64")

	ivb = dstream.Regroup(ivb, "DriverTrip", false)
	start := time.Now()
	ivb = dstream.LeftJoin(ivb, summary, []string{"DriverTrip", "DriverTrip"}, []string{"SummaryDistance"})
	elapsed := time.Since(start).Minutes()
	fmt.Printf("--- Finished joining trip summary data to main data stream ---\n")
	fmt.Printf("---\tElapsed time: %v minutes ---\n", elapsed)

	fmt.Printf("--- Starting data transformations and DOC fit --- \n")
	start = time.Now()
	// Divide into segments with the same trip and fixed time
	// deltas, drop when the time delta is not 100 milliseconds
	ivb.Reset()
	ivb = dstream.Regroup(ivb, "DriverTrip", false)
	ivb = dstream.DiffChunk(ivb, map[string]int{"Time": 1})
	ivb = dstream.Segment(ivb, []string{"DriverTrip", "Time$d1"})
	ivb = dstream.Filter(ivb, map[string]dstream.FilterFunc{"Time$d1": selectEq(10)})

	// lagged variables within the current chunks,
	// where the time deltas are the same
	// lagging removes first m=max(lags) rows removed from each chunk
	ivb = dstream.LagChunk(ivb, map[string]int{"Speed": maxSpeedLag,
		"FcwRange": maxRangeLag})

	// Drop consecutive brake points after the first,
	// require FCW to be active and minimum speed of 7 meters per second
	// total distance for trip > 0
	ivb = dstream.Apply(ivb, "brake2", fbrake, "float64")
	ivb = dstream.Filter(ivb, map[string]dstream.FilterFunc{"brake2": selectEq(0),
	      			  				"SummaryDistance": selectGtNotNaN(0)})
	
		//"FcwValidTarget": selectEq(1), "Speed[0]": selectGt(7),

	// keep driver, trip, time
	//ivb = dstream.DropCols(ivb, []string{"DriverTrip", "DriverTripTime", "Time$d1", "FcwValidTarget", "brake2", "SummaryDistance"})
	ivb = dstream.DropCols(ivb, []string{"DriverTripTime", "Time$d1", "brake2", "SummaryDistance"})

	ivb.Reset()
	ivb = dstream.Segment(ivb, []string{"DriverTrip"})
	fmt.Printf("Variable names after transformations: %v\n", ivb.Names())
	
	var regxnames []string
	for j := maxSpeedLag; j >= 0; j-- {
		regxnames = append(regxnames, fmt.Sprintf("Speed[%d]", -j))
	}
	for j := maxRangeLag; j >= 0; j-- {
		regxnames = append(regxnames, fmt.Sprintf("FcwRange[%d]", -j))
	}	
	
	pdim := len(regxnames)
	
	
	resfile, err := os.Create("/scratch/stats_flux/luers/segment_summaries.txt")
	if err != nil {
	   panic(err)
	}
	defer resfile.Close()
	w := bufio.NewWriter(resfile)

	_, err = w.WriteString("Driver\tTrip\tfcwvalid\tspeed_small\tspeed_large\tBrake\tcount\t")
	if err != nil {
	   panic(err)
	}
	for _, na := range regxnames {
	    w.WriteString(na + "_sum\t")
	}
	w.WriteString("\n")
	w.Flush()
	
	total_invalid := 0
	total_bigSpeed_b0 := 0
	total_bigSpeed_b1 := 0
	total_smallSpeed_b1 := 0
	total_smallSpeed_b0 := 0
	drivertrips_all := make(map[uint64]int)
	drivertrips_hasvalid := make(map[uint64]int)
	for d := 1; d <= maxDriverID; d++{
	    drivertrips_all[uint64(d)] = 0
	    drivertrips_hasvalid[uint64(d)] = 0
	}
	fmt.Printf("\n---Streaming through three-second segments to compute summary statistics---\n")

	start = time.Now()
	for ivb.Next() {
	    hasvalid := false
	    cDriver := ivb.Get("Driver").([]uint64)[0]
	    drivertrips_all[cDriver]++
	    cTrip := ivb.Get("Trip").([]float64)[0]
	    fcwvalid := ivb.Get("FcwValidTarget").([]float64)
	    braking := ivb.Get("Brake").([]float64)
	    cn := len(braking)
	    cn_invalid := 0
	    xMat := make([][]float64, pdim)
	    speed0 := ivb.Get("Speed[0]").([]float64)
	    for j, na := range regxnames {
	    	xMat[j] = ivb.Get(na).([]float64)
	    }
	    n_bigSpeed_b1 := 0
	    n_smallSpeed_b1 := 0
	    n_bigSpeed_b0 := 0
	    n_smallSpeed_b0 := 0
	    sums_smallSpeed_b1 := make([]float64, pdim) // sum of each variable for the low-speed excluded segments
	    sums_bigSpeed_b1 := make([]float64, pdim) // sum of each variable for the included
	    sums_smallSpeed_b0 := make([]float64, pdim)
	    sums_bigSpeed_b0 := make([]float64, pdim)
	    for i := 0; i < cn; i++ {
	    	if fcwvalid[i] > 0 {      // Lead vehicle detected
		   hasvalid = true
		   if speed0[i] > minSpeed {     // meets speed requirement
		      if braking[i] >0 { // and brake==1
		      	 for j := 0; j < pdim; j++ {
			     sums_bigSpeed_b1[j] += xMat[j][i]
			 }
		      	 n_bigSpeed_b1++
		      } else {    //brake==0
		      	for j := 0; j < pdim; j++ {
			    sums_bigSpeed_b0[j] += xMat[j][i]
			}
		      	n_bigSpeed_b0++
		      }
		   } else {     // small speed, EXCLUDED from analysis
		     if braking[i] > 0 { // valid range, small speed, brake==1
		     	for j := 0; j < pdim; j++ {
		     	    sums_smallSpeed_b1[j] += xMat[j][i]
			}
		     	n_smallSpeed_b1++
		     } else{ //valid range, small speed, brake==0
		       for j := 0; j < pdim; j++{
		       	   sums_smallSpeed_b0[j] += xMat[j][i]
		       }
		       n_smallSpeed_b0++
		     }
		   }

		} else {    // no lead vehicle detected, just count these segments
		  cn_invalid++ 
		}

	    }
	    
	    if hasvalid {
	       drivertrips_hasvalid[cDriver]++
	    }


	    dtripstr := strconv.FormatUint(cDriver, 10) + "\t" + strconv.FormatFloat(cTrip, 'g', 0, 64) + "\t"
	    w.WriteString(dtripstr) // driver and trip ID
	    w.WriteString("TRUE\tTRUE\tFALSE\t1\t") // valid FCW, small speed, brake==1
	    w.WriteString(strconv.Itoa(n_smallSpeed_b1))
	    sumstr := ""
	    for j := range regxnames {
	    	sumstr += "\t" + strconv.FormatFloat(sums_smallSpeed_b1[j], 'f', -1, 64)
	    }
	    w.WriteString(sumstr + "\n") 
	    w.Flush()

	    w.WriteString(dtripstr)
	    w.WriteString("TRUE\tTRUE\tFALSE\t0\t") // valid FCW, small speed, brake==0
	    w.WriteString(strconv.Itoa(n_smallSpeed_b0))
	    sumstr = ""
	    for j := range regxnames {
	    	sumstr += "\t" + strconv.FormatFloat(sums_smallSpeed_b0[j], 'f', -1, 64)
	    }
	    w.WriteString(sumstr + "\n") 
	    w.Flush()

	    w.WriteString(dtripstr)
	    w.WriteString("TRUE\tFALSE\tTRUE\t1\t") // valid FCW, large speed, brake==1
	    w.WriteString(strconv.Itoa(n_bigSpeed_b1))
	    sumstr = ""
	    for j := range regxnames {
	    	sumstr += "\t" + strconv.FormatFloat(sums_bigSpeed_b1[j], 'f', -1, 64)
	    }
	    w.WriteString(sumstr + "\n") 
	    w.Flush()

	    w.WriteString(dtripstr)
	    w.WriteString("TRUE\tFALSE\tTRUE\t1\t") // valid FCW, large speed, brake==0
	    w.WriteString(strconv.Itoa(n_bigSpeed_b0))
	    sumstr = ""
	    for j := range regxnames {
	    	sumstr += "\t" + strconv.FormatFloat(sums_bigSpeed_b0[j], 'f', -1, 64)
	    }
	    w.WriteString(sumstr + "\n") 
	    w.Flush()

	    w.WriteString(dtripstr)
	    w.WriteString("FALSE\tNA\tNA\tNA\t") // INVALID FCW
	    w.WriteString(strconv.Itoa(cn_invalid))
	    sumstr = ""
	    for j := 0; j < pdim; j++ {
	    	sumstr += "\tNA"
	    }
	    w.WriteString(sumstr + "\n") 
	    w.Flush()

	    total_invalid += cn_invalid
	    total_bigSpeed_b0 += n_bigSpeed_b0
	    total_bigSpeed_b1 += n_bigSpeed_b1
	    total_smallSpeed_b0 += n_smallSpeed_b0
	    total_smallSpeed_b1 += n_smallSpeed_b1
	}
	elapsed = time.Since(start).Minutes()
	fmt.Printf("---\tElapsed time: %v minutes ---\n\n", elapsed)
	fmt.Printf("Total observations, including FcwValid==0 and low-speed segments: %v\n", ivb.NumObs())
	fmt.Printf("Number segments FcwValid==0: %v\n", total_invalid)
	fmt.Printf("Number segments FcwValid==1, Speed[0] >= %v and Brake==1: %v\n", minSpeed, total_bigSpeed_b1)
	fmt.Printf("Number segments FcwValid==1, Speed[0] >= %v and Brake==0: %v\n", minSpeed, total_bigSpeed_b0)
	fmt.Printf("Number segments FcwValid==1, Speed[0] < %v and Brake==1: %v\n", minSpeed, total_smallSpeed_b1)
	fmt.Printf("Number segments FcwValid==1, Speed[0] < %v and Brake==0: %v\n", minSpeed, total_smallSpeed_b0)
	fmt.Printf("Driver\ttotal_trips_include_norange\n")
	for d, nt := range drivertrips_all {
	    fmt.Printf(strconv.FormatUint(d, 10) + "\t%v\n", nt)
	}
	fmt.Printf("\n------\n")
	fmt.Printf("Driver\ttotal_trips_hasvalid\n")
	for d, nt := range drivertrips_hasvalid {
	    fmt.Printf(strconv.FormatUint(d, 10) + "\t%v\n", nt)
	}
	
}

