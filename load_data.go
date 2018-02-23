package main

import (
	//"bufio"
	//"encoding/csv"
	"fmt"
	"github.com/brookluers/dstream/dstream"
	"github.com/gonum/matrix/mat64"
	"math"
	"os"
)

const (
	minCount  int     = 100
	nrow      int     = 100
	ncol      int     = 100
	minSpeed  float64 = 7.0
	chunkSize int     = 10000
	respvar   string  = "Brake_1sec"
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

func prodCols(c1, c2 string) dstream.ApplyFunc {
     f := func(v map[string]interface{}, z interface{}) {
       v1 := v[c1].([]float64)
       v2 := v[c2].([]float64)
       res := z.([]float64)
       for i := 0; i < len(res); i++ {
       	   res[i] = v1[i] * v2[i]
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

func imax(a, b int) int {
	if a > b {
		return a
	}
	return b
}

func lagbrakelabel(nlags int) dstream.ApplyFunc {
	f := func(v map[string]interface{}, z interface{}) {
		b := v["Brake"].([]float64) // braking indicator
		time := v["Time"].([]float64)
		var tdiff, tprev float64
		y := z.([]float64) // destination
		for i := 0; i < len(y); i++ {
			y[i] = 0
			if b[i] == 1 {
				y[i] = 1
				tprev = time[i]
				for j := i - 1; j >= imax(0, i-nlags); j-- {
					tdiff = tprev - time[j]
					if tdiff != 10 {
						break
					}
					y[j] = 1
					tprev = time[j]
				}
			}
		}
	}
	return f
}

func diagMat(v []float64) *mat64.Dense {
	n := len(v)
	m := mat64.NewDense(n, n, nil)
	for i := 0; i < n; i++ {
		m.Set(i, i, v[i])
	}
	return m
}

func loadSummary() dstream.Dstream {
	srdr, err := os.Open("data/summary.txt")
	// ("/nfs/turbo/ivbss/LvFot/summary.txt")
	//
	if err != nil {
		panic(err)
	}

	summary := dstream.FromCSV(srdr).SetFloatVars([]string{"Driver", "Trip", "StartTime", "EndTime", "Distance", "TODTripStart"}).SetChunkSize(5000).HasHeader().Done()
	summary = dstream.Apply(summary, "DriverTrip", driverTripId, "float64")
	summary = dstream.Convert(summary, "DriverTrip", "uint64")
	summary = dstream.Segment(summary, []string{"DriverTrip"})
	summary = dstream.Apply(summary, "SummaryDistance", applyIdent("Distance"), "float64")
	return summary
}

func loadDriverDat(id int, floatvars1 []string) dstream.Dstream {
	fname := fmt.Sprintf("data/small_%03d.txt", id)
	//("/nfs/turbo/ivbss/LvFot/data_%03d.txt", id)
	//

	fmt.Printf("Loading data from %v\n", fname)

	rdr, err := os.Open(fname)
	if err != nil {
		panic(err)
	}

	ivb := dstream.FromCSV(rdr).SetFloatVars(floatvars1).HasHeader().Done()

	// ID column for driver-trip
	ivb = dstream.Apply(ivb, "DriverTrip", driverTripId, "float64")
	ivb = dstream.Convert(ivb, "DriverTrip", "uint64")
	ivb = dstream.Apply(ivb, "DriverTripTime", driverTripTimeId, "float64")
	ivb = dstream.Convert(ivb, "Driver", "uint64")
	ivb = dstream.Regroup(ivb, "DriverTrip", false)

	summary := loadSummary()

	ivb = dstream.LeftJoin(ivb, summary, []string{"DriverTrip", "DriverTrip"}, []string{"SummaryDistance"})
	return ivb
}

func loadPoolDat(maxID int, floatvars1 []string) dstream.Dstream {
	fnames := make([]string, maxID)
	fmt.Printf("Loading data from:\n")
	for i := 1; i <= maxID; i++ {
		fnames[i-1] = fmt.Sprintf("data/small_%03d.txt", i)
		//  ("/nfs/turbo/ivbss/LvFot/data_%03d.txt", i)
		//

		fmt.Println(fnames[i-1])
	}

	// reader for main data files
	mrdr := dstream.NewMultiReadSeek0(fnames, true)

	ivb := dstream.FromCSV(mrdr).SetFloatVars(floatvars1).HasHeader().Done()

	// ID column for driver-trip
	ivb = dstream.Apply(ivb, "DriverTrip", driverTripId, "float64")
	ivb = dstream.Convert(ivb, "DriverTrip", "uint64")
	ivb = dstream.Apply(ivb, "DriverTripTime", driverTripTimeId, "float64")
	ivb = dstream.Convert(ivb, "Driver", "uint64")
	ivb = dstream.Regroup(ivb, "DriverTrip", false)
	summary := loadSummary()

	ivb = dstream.LeftJoin(ivb, summary, []string{"DriverTrip", "DriverTrip"}, []string{"SummaryDistance"})
	return ivb
}

func doTransforms(ivb dstream.Dstream, lagmap map[string]int) dstream.Dstream {
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
	ivb = dstream.LagChunk(ivb, lagmap)

	// Drop consecutive brake points after the first,
	// require FCW to be active and minimum speed of 7 meters per second
	// total distance for trip > 0
	ivb = dstream.Apply(ivb, "brake2", fbrake, "float64")
	ivb = dstream.Filter(ivb, map[string]dstream.FilterFunc{"brake2": selectEq(0),
		"FcwValidTarget": selectEq(1), "Speed[0]": selectGt(minSpeed),
		"SummaryDistance": selectGtNotNaN(0)})
	ivb = dstream.Apply(ivb, respvar, lagbrakelabel(10), "float64")
	// keep driver, trip, time
	ivb = dstream.DropCols(ivb, []string{"DriverTrip", "DriverTripTime", "Time$d1",
		"FcwValidTarget", "brake2", "SummaryDistance"})
	return ivb
}

func laggedInteraction(ivb dstream.Dstream, vn1 string, vn2 string, nlag int) dstream.Dstream {
     
     for i := 0; i <= nlag; i++{
     	 cn1 := fmt.Sprintf("%v[%d]", vn1, -i)
	 fmt.Printf("cn1 = %v\n", cn1)
	 cn2 := fmt.Sprintf("%v[%d]", vn2, -i)
	 fmt.Printf("cn2 = %v\n", cn2)
      	 ivb = dstream.Apply(ivb, fmt.Sprintf("%v%v[%d]", vn1, vn2, -i), prodCols(cn1, cn2), "float64")
     }

     return ivb
}
