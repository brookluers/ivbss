package main

import (
	"encoding/csv"
	"fmt"
	"github.com/brookluers/dimred"
	"github.com/brookluers/dstream/dstream"
	"github.com/gonum/matrix/mat64"
	"math"
	"os"
	"sort"
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
		"FcwValidTarget": selectEq(1), "Speed[0]": selectGt(7),
		"SummaryDistance": selectGtNotNaN(0)})

	// keep driver, trip, time
	ivb = dstream.DropCols(ivb, []string{"DriverTrip", "DriverTripTime", "Time$d1", "FcwValidTarget", "brake2", "SummaryDistance"})

	ivb.Reset()

	fmt.Printf("Variable names after transformations: %v\n", ivb.Names())

	// ---------- Fitting DOC -------------
	var regxnames []string
	for j := maxSpeedLag; j >= 0; j-- {
		regxnames = append(regxnames, fmt.Sprintf("Speed[%d]", -j))
	}
	for j := maxRangeLag; j >= 0; j-- {
		regxnames = append(regxnames, fmt.Sprintf("FcwRange[%d]", -j))
	}

	ivr0 := dstream.NewReg(ivb, "Brake", regxnames, "", "") // no projection onto PCs
	doc0 := dimred.NewDOC(ivr0)
	doc0.SetLogFile("log_noPC.txt")
	doc0.Init()
	ndir := 3
	npc := 8

	doc0.Fit(ndir)

	elapsed = time.Since(start).Minutes()
	fmt.Printf("--- Finished transforming/filtering/reading data and computed DOC directions for %d drivers\n\tElapsed time: %v minutes ---\n", maxDriverID, elapsed)

	// ---- Save the directions from DOC without any PC projections ---
	var dirnames []string
	dirnames = append(dirnames, "varname", "meandir")
	dirs0 := make([][]float64, 1+ndir+npc)
	dirs0[0] = doc0.MeanDir()
	for j := 0; j < ndir; j++ {
		dirs0[1+j] = doc0.CovDir(j)
		dirnames = append(dirnames, fmt.Sprintf("cd%d", j+1))
	}

	pdim := len(regxnames)
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

	fmt.Printf("-----Marginal Covariance-----\n%v\n", margcov)
	margcorr := mat64.NewSymDense(pdim, t3)
	fmt.Printf("-----Marginal Correlation-----\n%v\n", margcorr)
	fmt.Printf("--- mean of X, Y = 1 --- \n%v\n", doc0.GetMean(1))
	fmt.Printf("--- mean of X, Y = 0 --- \n%v\n", doc0.GetMean(0))

	es := new(mat64.EigenSym)
	ok := es.Factorize(margcorr, true)
	if !ok {
		panic("can't factorize marginal correlation matrix\n")
	}
	marg_evals := es.Values(nil)
	sort.Float64s(marg_evals)
	evec := new(mat64.Dense)
	evec.EigenvectorsSym(es)

	var esum float64 = 0.0
	eprop := make([]float64, npc)
	for _, v := range marg_evals {
		esum += v
	}

	pcMat := evec.View(0, pdim-npc, pdim, npc) // PCs can be applied to standardized x
	pcCoefs := make([][]float64, npc)
	pcnames := make([]string, npc)
	pcMatDense := mat64.DenseCopyOf(pcMat)
	pcMatDense.Mul(sd_inv_Diag, pcMatDense) // these directions can be applied to raw x

	// eigenvalues sorted in increasing order
	for j := 0; j < npc; j++ {
		pcCoefs[j] = mat64.Col(nil, j, pcMatDense)
		dirs0[1+ndir+j] = mat64.Col(nil, j, pcMatDense)
		eprop[j] = marg_evals[len(marg_evals)-1-j] / esum
		fmt.Printf("Proportion of variance, PC %d: %v\n", j, eprop[j])
		dirnames = append(dirnames, fmt.Sprintf("pc%d", j))
		pcnames[j] = fmt.Sprintf("%s%d", "pc", j)
	}

	dirFile, err := os.Create("/scratch/stats_flux/luers/directions_0pc.txt")
	if err != nil {
		panic(err)
	}
	wDir := csv.NewWriter(dirFile)
	if err := wDir.Write(dirnames); err != nil {
		panic(err)
	}
	drec := make([]string, len(dirnames))
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

	vm_pc := make(map[string]int)

	pcCoefs_expand := make([][]float64, npc)
	for k, a := range ivb.Names() {
		vm_pc[a] = k
	}
	for j := 0; j < npc; j++ {
		x := make([]float64, len(vm_pc)) // matches num. columns of ivb
		for k, na := range regxnames {
			x[vm_pc[na]] = pcCoefs[j][k]
		}
		pcCoefs_expand[j] = x
	}

	ivb = dstream.Linapply(ivb, pcCoefs_expand, "pc")
	//ivb = dstream.DropCols(ivb, regxnames)

	for cpc := npc; cpc >= ndir+1; cpc-- {
		fmt.Printf("---Projecting against %d PCs, then fitting DOC---\n", cpc)
		ivb.Reset()
		start = time.Now()
		curnames := make([]string, cpc)
		for j := 0; j < cpc; j++ {
			curnames[j] = pcnames[j]
		}
		fmt.Printf("Current X names: %v\n", curnames)
		ivr := dstream.NewReg(ivb, "Brake", curnames, "", "")
		doc := dimred.NewDOC(ivr)
		doc.SetLogFile(fmt.Sprintf("log_pc%d.txt", cpc))
		doc.Init()
		doc.Fit(ndir)
		elapsed = time.Since(start).Minutes()
		fmt.Printf("----Finished computing DOC directions using %d PCs for %d drivers. Elapsed time %v minutes---\n", cpc, maxDriverID, elapsed)

		//---------Convert back to raw X coordinates------------
		cDirs := make([][]float64, 1+ndir)
		var cDirsFlat []float64
		cDirs[0] = doc.MeanDir()
		cDirsFlat = append(cDirsFlat, cDirs[0]...)
		for j := 0; j < ndir; j++ {
			cDirs[1+j] = doc.CovDir(j)
			cDirsFlat = append(cDirsFlat, cDirs[1+j]...)
		}
		cMat := mat64.NewDense(ndir+1, cpc, cDirsFlat)           // matrix of coefficients in standardized + PC-projected X coordinates
		cCoefMat := mat64.NewDense(pdim, ndir+1, nil)            // matrix of coefficients in raw X coordinates
		cCoefMat.Mul(pcMatDense.View(0, 0, pdim, cpc), cMat.T()) // each column is a vector of coefficients for the original X variables
		//-------------

		//---------------Save current coefficients in raw X coordinates ---------------
		cFile, err := os.Create(fmt.Sprintf("/scratch/stats_flux/luers/directions_%dpc.txt", cpc))
		if err != nil {
			panic(err)
		}
		cWr := csv.NewWriter(cFile)
		if err := cWr.Write(dirnames[0:(ndir + 2)]); err != nil { // first line contains names of directions
			panic(err)
		}
		cRec := make([]string, 2+ndir)

		for k, na := range regxnames {
			cRow := mat64.Row(nil, k, cCoefMat)
			cRec[0] = na
			for j := 0; j < len(cDirs); j++ {
				cRec[1+j] = fmt.Sprintf("%v", cRow[j])
			}
			if err := cWr.Write(cRec); err != nil {
				panic(err)
			}
			cWr.Flush()
		}
		cFile.Close()
		//-----------------------------------------------------------

		// ----------------- Save projected data using npc PCs---------
		if cpc == npc {
			vm := make(map[string]int)                   // map variable names to column positions
			dirs_expand := make([][]float64, len(cDirs)) // mean direction, doc directions, PC directions

			for k, a := range ivb.Names() {
				vm[a] = k
			}

			//******Project against mean direction*************
			for j := 0; j < 1; j++ {
				x := make([]float64, len(vm)) // length will be longer than current number of PCs
				for k, na := range ivr.XNames() {
					x[vm[na]] = cDirs[j][k] // coefficients
				}
				dirs_expand[j] = x // x has zeroes in position of non-X variables
			}
			ivb.Reset()
			ivb = dstream.Linapply(ivb, dirs_expand[0:1], "meandir")

			//********* Project against covariance directions***********
			for k, a := range ivb.Names() {
				vm[a] = k
			}

			for j := 1; j < ndir+1; j++ {
				x := make([]float64, len(vm))
				for k, na := range ivr.XNames() {
					x[vm[na]] = cDirs[j][k]
				}
				dirs_expand[j] = x
			}
			ivb = dstream.Linapply(ivb, dirs_expand[1:(ndir+1)], "cd")

			ivb.Reset()
			ivb = dstream.Regroup(ivb, "Driver", false)

			pfilename := fmt.Sprintf("smproj_%dpc_", cpc)
			err := dstream.ToCSV(ivb).DoneByChunk("Driver", "%03d", pfilename, ".txt")
			if err != nil {
				panic(err)
			}

			//**********Compute speed-range averages on projected axes
			ivb.Reset()
			var meandir_window1, meandir_window2, cd0_w1, cd0_w2 dstream.Dstream
			meandir_window1 = dstream.Filter(ivb, map[string]dstream.FilterFunc{"meandir0": selectBt(1.0, 1.25)})
			meandir_window2 = dstream.Filter(ivb, map[string]dstream.FilterFunc{"meandir0": selectBt(2.75, 3.0)})

			cd0_w1 = dstream.Filter(ivb, map[string]dstream.FilterFunc{"cd0": selectBt(-0.75,-0.25)})
			cd0_w2 = dstream.Filter(ivb, map[string]dstream.FilterFunc{"cd0": selectBt(0.25, 0.75)})

			sums_window1 := make([]float64, len(regxnames))
			means_window1 := make([]float64, len(regxnames))
			n_window1 := 0
			for meandir_window1.Next() {
				cn := len(meandir_window1.Get("Brake").([]float64))
				n_window1 += cn
				for k := 0; k < len(regxnames); k++ {
					x := meandir_window1.Get(regxnames[k]).([]float64)
					for i := 0; i < cn; i++ {
						sums_window1[k] += x[i]
					}
				}
			}
			for k := 0; k < len(regxnames); k++ {
				means_window1[k] = sums_window1[k] / float64(n_window1)
			}
			meandir_window1.Reset()
			ivb.Reset()

			sums_window2 := make([]float64, len(regxnames))
			means_window2 := make([]float64, len(regxnames))
			n_window2 := 0
			for meandir_window2.Next() {
				cn := len(meandir_window2.Get("Brake").([]float64))
				n_window2 += cn
				for k := 0; k < len(regxnames); k++ {
					x := meandir_window2.Get(regxnames[k]).([]float64)
					for i := 0; i < cn; i++ {
						sums_window2[k] += x[i]
					}
				}
			}
			for k := 0; k < len(regxnames); k++ {
				means_window2[k] = sums_window2[k] / float64(n_window2)
			}
			meandir_window2.Reset()
			ivb.Reset()

			sums_cd0_w1 := make([]float64, len(regxnames))
			means_cd0_w1 := make([]float64, len(regxnames))
			n_cd0_w1 := 0
			for cd0_w1.Next() {
				cn := len(cd0_w1.Get("Brake").([]float64))
				n_cd0_w1 += cn
				for k := 0; k < len(regxnames); k++ {
					x := cd0_w1.Get(regxnames[k]).([]float64)
					for i := 0; i < cn; i++ {
						sums_cd0_w1[k] += x[i]
					}
				}
			}
			for k := 0; k < len(regxnames); k++ {
				means_cd0_w1[k] = sums_cd0_w1[k] / float64(n_cd0_w1)
			}
			cd0_w1.Reset()
			ivb.Reset()

			sums_cd0_w2 := make([]float64, len(regxnames))
			means_cd0_w2 := make([]float64, len(regxnames))
			n_cd0_w2 := 0
			for cd0_w2.Next() {
				cn := len(cd0_w2.Get("Brake").([]float64))
				n_cd0_w2 += cn
				for k := 0; k < len(regxnames); k++ {
					x := cd0_w2.Get(regxnames[k]).([]float64)
					for i := 0; i < cn; i++ {
						sums_cd0_w2[k] += x[i]
					}
				}
			}
			for k := 0; k < len(regxnames); k++ {
				means_cd0_w2[k] = sums_cd0_w2[k] / float64(n_cd0_w2)
			}
			cd0_w2.Reset()
			ivb.Reset()

			fmt.Printf("Mean of X for meandir between -0.75 and -0.25: \n%v\n", means_window1)
			fmt.Printf("Mean of X for meandir between 0.25 and 0.75: \n%v\n", means_window2)
			fmt.Printf("Mean of X for cd0 between -4 and -2: \n%v\n", means_cd0_w1)
			fmt.Printf("Mean of X for cd0 between 2 and 4: \n%v\n", means_cd0_w2)
			fmt.Printf("Num. observations for meandir between -0.75 and -0.25: %v\n", n_window1)
			fmt.Printf("Num. observations for meandir between 0.25 and 0.75: %v\n", n_window2)
			fmt.Printf("Num. observations for cd0 between -4 and -2: %v\n", n_cd0_w1)
			fmt.Printf("Num. observations for cd0 between 2 and 4: %v\n", n_cd0_w2)

		}
		//-----------------------------------------------------------------

	}

}
