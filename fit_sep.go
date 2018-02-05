package main

import (
	"encoding/csv"
	"fmt"
	"bufio"
	"github.com/brookluers/dimred"
	"github.com/brookluers/dstream/dstream"
	"github.com/gonum/matrix/mat64"
	"math"
	"os"
	"sort"
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
	  for i := 0; i < len(y); i++{
	      y[i] = 0
	      if (b[i] == 1) {
	      	 y[i] = 1
	      	 tprev = time[i]
	      	 for j := i - 1; j >= imax(0, i - nlags); j--{
		     tdiff = tprev - time[j]
		     if (tdiff != 10) { break }
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

func main() {
     	   
	maxDriverID := 3
	fnames := make([]string, maxDriverID)

	// fetch summary data for each driver-trip
	srdr, err := os.Open("data/summary.txt")
	      //  ("/nfs/turbo/ivbss/LvFot/summary.txt")
	if err != nil {
		panic(err)
	}

	summary := dstream.FromCSV(srdr).SetFloatVars([]string{"Driver", "Trip", "StartTime", "EndTime", "Distance", "TODTripStart"}).SetChunkSize(5000).HasHeader().Done()

	summary = dstream.Apply(summary, "DriverTrip", driverTripId, "float64")
	summary = dstream.Convert(summary, "DriverTrip", "uint64")
	summary = dstream.Segment(summary, []string{"DriverTrip"})
	summary = dstream.Apply(summary, "SummaryDistance", applyIdent("Distance"), "float64")

	var regxnames []string
	for j := maxSpeedLag; j >= 0; j-- {
		regxnames = append(regxnames, fmt.Sprintf("Speed[%d]", -j))
	}
	for j := maxRangeLag; j >= 0; j-- {
		regxnames = append(regxnames, fmt.Sprintf("FcwRange[%d]", -j))
	}

	pdim := len(regxnames)
	npc := 8
	ndir := 10

	for i := 1; i <= maxDriverID; i++ {
		fnames[i-1] = fmt.Sprintf("data/small_%03d.txt", i)
			    //("/nfs/turbo/ivbss/LvFot/data_%03d.txt", i)
		fmt.Println(fnames[i-1])

		rdr, err := os.Open(fnames[i-1])
		if err != nil {
			panic(err)
		}

		ivb := dstream.FromCSV(rdr).SetFloatVars([]string{"Driver", "Trip",
			"Time", "Speed", "Brake",
			"FcwValidTarget", "FcwRange"}).HasHeader().Done()

		// ID column for driver-trip
		ivb = dstream.Apply(ivb, "DriverTrip", driverTripId, "float64")
		ivb = dstream.Convert(ivb, "DriverTrip", "uint64")
		ivb = dstream.Apply(ivb, "DriverTripTime", driverTripTimeId, "float64")
		ivb = dstream.Convert(ivb, "Driver", "uint64")
		ivb = dstream.LeftJoin(ivb, summary, []string{"DriverTrip", "DriverTrip"}, []string{"SummaryDistance"})

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
			"FcwValidTarget": selectEq(1), "Speed[0]": selectGt(7),
			"SummaryDistance": selectGtNotNaN(0)})
		ivb = dstream.Apply(ivb, "Brake_1sec", lagbrakelabel(10), "float64")

		// keep driver, trip, time
		ivb = dstream.DropCols(ivb, []string{"DriverTrip", "DriverTripTime", "Time$d1",
			"FcwValidTarget", "brake2", "SummaryDistance"})
		ivb.Reset()

		
		// ---------- Fitting DOC -------------

		ivr0 :=  dstream.DropCols(ivb, []string{"Driver", "Brake", "Trip", "Time"})

		fmt.Printf("Variable names for DOC: %v\n", ivr0.Names())
		     //dstream.NewReg(ivb, "Brake_1sec", regxnames, "", "")
		fmt.Printf("Variable names for original dstream: %v\n", ivb.Names())
		doc0 := dimred.NewDOC(ivr0, "Brake_1sec")
		doc0.SetLogFile(fmt.Sprintf("log%03d.txt", i))
		doc0.Done()
		doc0.Fit(ndir)

		covfile, err := os.Create(fmt.Sprintf("data/sep_cov_small_%03d.txt", i))
			     // ("/scratch/stats_flux/luers/sep_cov_%03d.txt", i))
		if err != nil {
		   panic(err)
		}
		wCov := bufio.NewWriter(covfile)
		
		margcov := mat64.NewSymDense(pdim, doc0.GetMargCov())
		marg_sd_inv := make([]float64, pdim)
		for k := 0; k < pdim; k++ {
			marg_sd_inv[k] = math.Pow(margcov.At(k, k), -0.5)
		}
		sd_inv_Diag := diagMat(marg_sd_inv)
		t1 := mat64.NewDense(pdim, pdim, nil)
		t1.Mul(sd_inv_Diag, margcov)
		t2 := mat64.NewDense(pdim, pdim, nil)
		t2.Mul(t1, sd_inv_Diag)
		var t3 []float64
		for k := 0; k < pdim; k++ {
			t3 = append(t3, mat64.Row(nil, k, t2)...)
		}

		//  fmt.Printf("-----Marginal Covariance-----\n%v\n", margcov)
		margcorr := mat64.NewSymDense(pdim, t3)
		//  fmt.Printf("-----Marginal Correlation-----\n%v\n", margcorr)
		es := new(mat64.EigenSym)
		ok := es.Factorize(margcorr, true) //margcov, true)
		if !ok {
			panic("can't factorize marginal correlation matrix\n")
		}
		marg_evals := es.Values(nil)
		sort.Float64s(marg_evals)
		_, err = fmt.Fprintf(wCov, "%v\n", regxnames)
		if err != nil {
		   panic(err)
		}		
		_, err = fmt.Fprintf(wCov, "marginal eigenvalues\n%v\n", marg_evals)
		if err != nil {
		   panic(err)
		}
		fmt_margcov := mat64.Formatted(margcov, mat64.Prefix("\n     "), mat64.Squeeze())
		_, err = fmt.Fprintf(wCov, "marginal covariance\n%v\n", fmt_margcov)
		if err != nil {
		   panic(err)
		}
		_, err = fmt.Fprintf(wCov, "DOC eigenvalues\n%v\n", doc0.Eig())
		if err != nil {
		   panic(err)
		}
		cov1 := mat64.NewSymDense(pdim, doc0.GetCov(1))
		cov0 := mat64.NewSymDense(pdim, doc0.GetCov(0))
		fmt_cov1 := mat64.Formatted(cov1, mat64.Prefix("     "), mat64.Squeeze())
		fmt_cov0 := mat64.Formatted(cov0, mat64.Prefix("     "), mat64.Squeeze())
		_, err = fmt.Fprintf(wCov, "covariance, y=1\n%v\n", fmt_cov1)
		if err != nil {
		   panic(err)
		}
		_, err = fmt.Fprintf(wCov, "covariance, y=0\n%v\n", fmt_cov0)
		if err != nil {
		   panic(err)
		}
		wCov.Flush()		
		
		evec := new(mat64.Dense)
		evec.EigenvectorsSym(es)
		dirFile, err := os.Create(fmt.Sprintf("data/dir_sep_small_%03d.txt", i))
			     //  ("/scratch/stats_flux/luers/directions_sep_%03d.txt", i))
		if err != nil {
			panic(err)
		}
		wDir := csv.NewWriter(dirFile)
		var temp []string
		temp = append(temp, "varname", "meandir")

		dirs0 := make([][]float64, 1+ndir+npc)
		dirs0[0] = doc0.MeanDir()
		for j := 0; j < ndir; j++ {
			//dirs0[1 + j] = make([]float64, pdim)
			dirs0[1+j] = doc0.CovDir(j)
			temp = append(temp, fmt.Sprintf("cd%d", j+1))
		}

		pcMat := evec.View(0, pdim-npc, pdim, npc) // PCs can be applied to standardized x
		pcMatDense := mat64.DenseCopyOf(pcMat)
		pcMatDense.Mul(sd_inv_Diag, pcMatDense) // these directions can be applied to raw x
		// eigenvalues sorted in increasing order
		for j := 0; j < npc; j++ {
			dirs0[1+ndir+j] = mat64.Col(nil, j, pcMat)
			temp = append(temp, fmt.Sprintf("pc%d", j+1))
		}

		// temp = {"varname", "meandir", "cd1"..."pc1"...}
		// len(temp): 2 + ndir + npc
		// dirs0: mean direction, covariance directions, pc directions
		if err := wDir.Write(temp); err != nil {
			panic(err)
		}

		// save coefficient vectors for each dimension reduction direction
		drec := make([]string, len(temp))
		for k, na := range regxnames {
			drec[0] = na
			for j := 0; j < len(dirs0); j++ {
				drec[1+j] = fmt.Sprintf("%v", dirs0[j][k])
			}
			if err := wDir.Write(drec); err != nil {
				panic(err)
			}
			wDir.Flush()
		}
		dirFile.Close()
	}

}
