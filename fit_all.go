package main

import (
	"encoding/csv"
	"fmt"
	"github.com/brookluers/dimred"
	"github.com/brookluers/dstream/dstream"
	"github.com/gonum/floats"
	"math"
	"os"
	//	"sort"
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

//driverTripTimeId is an ApplyFunc that creates
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

type matxyz struct {
	data []float64
	r    int
	c    int
}

func (m *matxyz) Dims() (int, int) {
	return m.r, m.c
}

func (m *matxyz) Z(c, r int) float64 {
	return m.data[r*m.c+c]
}

func (m *matxyz) X(c int) float64 {
	return float64(c)
}

func (m *matxyz) Y(r int) float64 {
	return float64(r)
}

func (m *matxyz) Min() float64 {
	return 0
}

func (m *matxyz) Max() float64 {
	return 1
}

func standardize(vec, mat []float64) {
	v := 0.0
	p := len(vec)
	for i := 0; i < p; i++ {
		for j := 0; j <= i; j++ {
			u := vec[i] * vec[j] * mat[i*p+j]
			if j != i {
				v += 2 * u
			} else {
				v += u
			}
		}
	}
	v = math.Sqrt(v)

	floats.Scale(1/v, vec)
}

func getPos(data dstream.Dstream, name string) int {
	for k, v := range data.Names() {
		if v == name {
			return k
		}
	}
	panic("cannot find " + name)
}

func main() {

	maxDriverID := 15 //108
	fnames := make([]string, maxDriverID)
	fnames2 := make([]string, maxDriverID)
	fmt.Println("File names:")
	for i := 1; i <= maxDriverID; i++ {
		fnames[i-1] = fnames[i-1] = fmt.Sprintf("/nfs/turbo/ivbss/LvFot/data_%03d.txt", i)
		fnames2[i-1] = fmt.Sprintf("/nfs/turbo/ivbss/LvFot/data2_%03d.txt", i) //
		fmt.Println(fnames[i-1])
		fmt.Println(fnames2[i-1])
	}

	// reader for main data files
	mrdr := dstream.NewMultiReadSeek0(fnames, true)
	mrdr2 := dstream.NewMultiReadSeek0(fnames2, true)

	ivb2 := dstream.FromCSV(mrdr2).SetFloatVars([]string{"Driver", "Trip", "Time", "OutsideTemperature"}).SetChunkSize(chunkSize).HasHeader().Done()
	ivb2 = dstream.Apply(ivb2, "DriverTripTime", driverTripTimeId, "float64")
	ivb2 = dstream.Segment(ivb2, []string{"DriverTripTime"})
	ivb2 = dstream.Convert(ivb2, "DriverTripTime", "uint64")
	ivb2 = dstream.DropCols(ivb2, []string{"Driver", "Trip", "Time"})

	// fetch summary data for each driver-trip
	srdr, err := os.Open("/nfs/turbo/ivbss/LvFot/summary.txt")
	if err != nil {
		panic(err)
	}
	srdr1, err := os.Open("/scratch/stats_flux/luers/summary_trip1_starttime.txt")
	if err != nil {
		panic(err)
	}
	summary1 := dstream.FromCSV(srdr1).SetFloatVars([]string{"Driver", "trip1starttime"}).SetChunkSize(117).HasHeader().Done()
	summary1 = dstream.Segment(summary1, []string{"Driver"})
	summary1 = dstream.Convert(summary1, "Driver", "uint64")
	summary := dstream.FromCSV(srdr).SetFloatVars([]string{"Driver", "Trip", "StartTime", "EndTime", "IvbssEnable", "Distance", "TODTripStart"}).SetChunkSize(5000).HasHeader().Done()

	summary = dstream.Apply(summary, "DriverTrip", driverTripId, "float64")
	summary = dstream.Convert(summary, "DriverTrip", "uint64")
	summary = dstream.Segment(summary, []string{"DriverTrip"})
	summary = dstream.Apply(summary, "SummaryDistance", applyIdent("Distance"), "float64")
	summary.Reset()

	// TODO: add "Wipers" variable, 0-1 indicator of wiper activity
	ivb := dstream.FromCSV(mrdr).SetFloatVars([]string{"Driver", "Trip", "Time", "Speed", "Brake", "FcwValidTarget", "FcwRange"}).HasHeader().Done()
	
	// ID column for driver-trip
	ivb = dstream.Apply(ivb, "DriverTrip", driverTripId, "float64")
	ivb = dstream.Convert(ivb, "DriverTrip", "uint64")
	ivb = dstream.Apply(ivb, "DriverTripTime", driverTripTimeId, "float64")
	// ID column for 10-hz measurement within each driver-trip
	ivb = dstream.Convert(ivb, "Driver", "uint64")
	ivb = dstream.Regroup(ivb, "Driver", false)
	
	// Fetch the start time of the first trip for each driver
	ivb = dstream.LeftJoin(ivb, summary1, []string{"Driver", "Driver"}, []string{"trip1starttime"})
	ivb.Reset()
	ivb = dstream.Regroup(ivb, "DriverTrip", true)
	// Fetch trip-level summary information: start time, whether the sensors are turned on
	// (IvbssEnable), total trip distance
	ivb = dstream.LeftJoin(ivb, summary, []string{"DriverTrip", "DriverTrip"}, []string{"StartTime", "IvbssEnable", "TODTripStart", "SummaryDistance"})
	ivb.Reset()
	ivb = dstream.Convert(ivb, "DriverTripTime", "uint64")
	ivb = dstream.Segment(ivb, []string{"DriverTripTime"})

	// Fetch 10-hz temperature 
	ivb = dstream.LeftJoin(ivb, ivb2, []string{"DriverTripTime", "DriverTripTime"}, []string{"OutsideTemperature"})

	// hundredths of a second since this trip started
	ivb = dstream.Apply(ivb, "TripElapsed", diffCols("Time", "StartTime"), "float64")

	hund2days := func(v map[string]interface{}, z interface{}) {
		el := v["TripElapsed"].([]float64)
		res := z.([]float64)
		var day100ths float64 = 8640000 // 100ths of a second in a day
		for i := 0; i < len(res); i++ {
			res[i] = el[i] / day100ths
		}
	}

	// days (fractional) since this trip started
	ivb = dstream.Apply(ivb, "TripElapsedDays", hund2days, "float64")

	// number of days (can be fractional) between the start of Trip 1
	// and the start of the current trip
	// subtracting two columns with units: days since December 30, 1899
	ivb = dstream.Apply(ivb, "TripStartStudyElapsed", diffCols("TODTripStart", "trip1starttime"), "float64")

	// elapsed time since start of study, i.e.
	// the number of days (can be fractional) between the start of Trip 1
	// and the current measurement
	ivb = dstream.Apply(ivb, "OnStudyElapsed", sumCols("TripStartStudyElapsed", "TripElapsedDays"), "float64")
	ivb = dstream.DropCols(ivb, []string{"StartTime", "TripElapsed", "TripElapsedDays", "TripStartStudyElapsed", "trip1starttime"})

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
	ivb = dstream.LagChunk(ivb, map[string]int{"Speed": maxSpeedLag, "FcwRange": maxRangeLag})

	// Drop consecutive brake points after the first,
	// require FCW to be active and minimum speed of 7 meters per second
	// total distance for trip > 0
	ivb = dstream.Apply(ivb, "brake2", fbrake, "float64")
	ivb = dstream.Filter(ivb, map[string]dstream.FilterFunc{"brake2": selectEq(0),
		"FcwValidTarget": selectEq(1), "Speed[0]": selectGt(7), "SummaryDistance": selectGt(0)})

	// keep driver, trip, time
	ivb = dstream.DropCols(ivb, []string{"DriverTrip", "DriverTripTime", "Time$d1", "FcwValidTarget", "brake2", "TODTripStart", "SummaryDistance", "IvbssEnable"})
	//ivb = dstream.Convert(ivb, "Driver", "float64") // added uint64 support to linapply
	fmt.Printf("\n-----ivb.NumObs() after transformation: %v---\n", ivb.NumObs())

	ivb.Reset()

	fmt.Printf("Variable names after transformations: %v\n", ivb.Names())

	// ---------- Fitting DOC -------------
	fmt.Printf("\nnumber of variables before fit: %v\n", len(ivb.Names()))
	var regxnames []string
	regxnames = append(regxnames, "OutsideTemperature", "OnStudyElapsed")
	for j := maxSpeedLag; j >= 0; j-- {
		regxnames = append(regxnames, fmt.Sprintf("Speed[%d]", -j))
	}
	for j := maxRangeLag; j >= 0; j-- {
		regxnames = append(regxnames, fmt.Sprintf("FcwRange[%d]", -j))
	}
	fmt.Printf("\nregxnames = %v\n", regxnames)

	ivr := dstream.NewReg(ivb, "Brake", regxnames, "", "")

	doc := dimred.NewDOC(ivr)
	doc.SetLogFile("log_multi.txt")
	doc.Init()

	ndir := 2
	doc.Fit(ndir)

	fmt.Printf("nobs after fit=%d\n", ivb.NumObs())

	dirs0 := [][]float64{doc.MeanDir(), doc.CovDir(0), doc.CovDir(1)}

	// Expand to match the data set
	vm := make(map[string]int) // map variable names to column positions
	dirs := make([][]float64, 3)
	for k, a := range ivb.Names() {
		vm[a] = k
	}
	for j := 0; j < 3; j++ {
		x := make([]float64, len(vm)) // length will be longer than len(ivr.XNames())
		for k, na := range ivr.XNames() {
			x[vm[na]] = dirs0[j][k] // coefficients
		}
		dirs[j] = x // x has zeroes in position of non-X variables
	}
	ivb.Reset()
	ivb = dstream.Linapply(ivb, dirs, "dr")
	//ivb.Reset()
	//md := dstream.GetCol(ivb, "dr0").([]float64)
	//ivb.Reset()
	//cd1 := dstream.GetCol(ivb, "dr1").([]float64)
	//ivb.Reset()
	//cd2 := dstream.GetCol(ivb, "dr2").([]float64)
	//ivb.Reset()
	//uy := dstream.GetCol(ivb, "Brake").([]float64)

	ivb.Reset()
	ivb = dstream.Regroup(ivb, "Driver", false)
	var curDriver uint64
	for ivb.Next() {

		curDriver = ivb.Get("Driver").([]uint64)[0]
		md := ivb.Get("dr0").([]float64)
		cd1 := ivb.Get("dr1").([]float64)
		cd2 := ivb.Get("dr2").([]float64)
		uy := ivb.Get("Brake").([]float64)
		fmt.Printf("\nWriting projected data for driver %d to disk\n", int(curDriver))
		pFile, err := os.Create(fmt.Sprintf("/scratch/stats_flux/luers/smproj_multi_%03d.txt", int(curDriver)))
		if err != nil {
			panic(err)
		}
		wCsvProj := csv.NewWriter(pFile)
		prec := make([]string, 5)
		if err := wCsvProj.Write([]string{"Driver", "meandir", "cd1", "cd2", "y"}); err != nil {
			panic(err)
		}

		for i := 0; i < len(md); i++ {
			prec[0] = fmt.Sprintf("%d", int(curDriver))
			prec[1] = fmt.Sprintf("%v", md[i])
			prec[2] = fmt.Sprintf("%v", cd1[i])
			prec[3] = fmt.Sprintf("%v", cd2[i])
			prec[4] = fmt.Sprintf("%v", uy[i])
			if err := wCsvProj.Write(prec); err != nil {
				panic(err)
			}
			wCsvProj.Flush()
		}
		pFile.Close()
	}

}
