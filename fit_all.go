package main

import (
	"encoding/csv"
	"fmt"
	"github.com/brookluers/dimred"
	"github.com/brookluers/dstream/dstream"
	"github.com/gonum/floats"
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
	chunkSize   int	    = 10000
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
	    res[i] = dr[i] * 1000000000 + tr[i] * 1000000 + ti[i]
	}
}

//driverTripId is an ApplyFunc that creates 
func driverTripId(v map[string]interface{}, z interface{}) {
	dr := v["Driver"].([]float64)
	tr := v["Trip"].([]float64)
	res := z.([]float64)
	for i := range dr {
	    res[i] = dr[i] * 1000000000 + tr[i] * 1000000
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

func cormat(x []float64, p int) []float64 {

	y := make([]float64, len(x))
	s := make([]float64, p)
	for i := 0; i < p; i++ {
		s[i] = math.Sqrt(x[i*p+i])
	}
	for i, v := range x {
		j1 := i / p
		j2 := i % p
		y[i] = v / (s[j1] * s[j2])
	}

	return y
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

func getScores(x, br []float64, w int) ([]float64, []float64) {

	ii := make([]int, len(x))
	floats.Argsort(x, ii)
	sort.Sort(sort.Float64Slice(x))

	var b []float64
	for _, i := range ii {
		b = append(b, br[i])
	}

	z := make([]float64, len(x))
	for i := w; i < len(b)-w; i++ {
		if b[i] == 1 {
			for j := i - w; j < i+w; j++ {
				z[j]++ // count the  number of Brake==1 in the window
			}
		}
	}

	for i, _ := range z {
		z[i] /= float64(2 * w) // each window has width 2 * w
	}

	return x[w : len(x)-w], z[w : len(z)-w]
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

	maxDriverID := 2 //108
	fnames := make([]string, maxDriverID)
	fnames2 := make([]string, maxDriverID)
	fmt.Println("File names:")
	for i := 1; i <= maxDriverID; i++ {
		fnames[i-1] = fmt.Sprintf("small_%03d.txt", i)//fnames[i-1] = fmt.Sprintf("/nfs/turbo/ivbss/LvFot/data_%03d.txt", i) 
		fnames2[i-1]= fmt.Sprintf("small2_%03d.txt", i) //
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
	ivb2 = dstream.DropCols(ivb2, []string{"Driver","Trip","Time"})

	dstream.ToCSV(ivb2).Filename("ivb2_multi.txt").Done()
	fmt.Printf("\n----wrote ivb2 to ivb2.txt---\n")
	fmt.Printf("\nivb2 number observations = %v\n", ivb2.NumObs())
	// fetch summary data for each driver-trip
	srdr, err := os.Open("summary.txt")
	if err != nil{
	   panic(err)
	}
	srdr1, err := os.Open("summary_trip1_starttime.txt")
	if err != nil {
	   panic(err)
	}
	summary1 := dstream.FromCSV(srdr1).SetFloatVars([]string{"Driver","trip1starttime"}).SetChunkSize(117).HasHeader().Done()
	summary1 = dstream.Segment(summary1, []string{"Driver"})
	summary1 = dstream.Convert(summary1, "Driver", "uint64")
 	summary := dstream.FromCSV(srdr).SetFloatVars([]string{"Driver","Trip","StartTime","EndTime","IvbssEnable","Distance","TODTripStart"}).SetChunkSize(5000).HasHeader().Done()

	summary = dstream.Apply(summary, "DriverTrip", driverTripId, "float64")
	summary = dstream.Convert(summary, "DriverTrip", "uint64")
	summary = dstream.Segment(summary, []string{"DriverTrip"})
	summary = dstream.Apply(summary, "SummaryDistance", applyIdent("Distance"), "float64")
	summary.Reset()

	ivb := dstream.FromCSV(mrdr).SetFloatVars([]string{"Driver", "Trip", "Time", "Speed", "Brake", "FcwValidTarget", "FcwRange","Wipers"}).HasHeader().Done()
	ivb = dstream.Apply(ivb, "DriverTrip", driverTripId, "float64")
	ivb = dstream.Convert(ivb, "DriverTrip", "uint64")
	ivb = dstream.Apply(ivb, "DriverTripTime", driverTripTimeId, "float64")
	ivb = dstream.Convert(ivb, "Driver", "uint64")
	ivb = dstream.Regroup(ivb, "Driver", false)
	ivb = dstream.LeftJoin(ivb, summary1, []string{"Driver","Driver"}, []string{"trip1starttime"})
	ivb.Reset()
	ivb = dstream.Regroup(ivb, "DriverTrip", true)
	ivb = dstream.LeftJoin(ivb, summary, []string{"DriverTrip","DriverTrip"}, []string{"StartTime","IvbssEnable","TODTripStart", "SummaryDistance"})
	ivb.Reset()
	ivb = dstream.Convert(ivb, "DriverTripTime", "uint64")
	ivb = dstream.Segment(ivb, []string{"DriverTripTime"})
	ivb = dstream.LeftJoin(ivb, ivb2, []string{"DriverTripTime","DriverTripTime"}, []string{"OutsideTemperature"})
	
	// hundredths of a second since this trip started
	ivb = dstream.Apply(ivb, "TripElapsed", diffCols("Time","StartTime"), "float64")

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
	dstream.ToCSV(ivb).Filename("small_multi_after_joins.txt").Done()

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
	//	"IvbssEnable": selectEq(1)})

	ivb = dstream.DropCols(ivb, []string{"Driver","DriverTrip","DriverTripTime","Trip", "Time", "Time$d1", "FcwValidTarget", "brake2", "TODTripStart", "SummaryDistance","IvbssEnable"})

	fmt.Printf("\n-----ivb.NumObs() after transformation: %v---\n", ivb.NumObs())
	dstream.ToCSV(ivb).Filename("small_multi_trans.txt").Done()
	
	ivb.Reset()
	//ivb = dstream.MemCopy(ivb)

	fmt.Printf("Variable names after transformations: %v\n", ivb.Names())
	
	// ---------- Fitting DOC -------------
	fmt.Printf("\nVariable names before fitting DOC: %v\n", ivb.Names())
	fmt.Printf("\nnumber of variables before fit: %v\n", len(ivb.Names()))
	ivr := dstream.NewReg(ivb, "Brake", nil, "", "")

	doc := dimred.NewDOC(ivr)
	doc.SetLogFile("log.txt")
	doc.Init()

	ndir := 2
	doc.Fit(ndir)

	fmt.Printf("nobs after fit=%d\n", ivb.NumObs())
	fmt.Printf("%v\n", ivr.XNames()[0:31])
	fmt.Printf("%v\n", ivr.XNames()[31:62])
	//fmt.Printf("%d %d %d\n", len(doc.YMean(0)), len(doc.MeanDir()), len(doc.CovDir(0)))

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
		dirs[j] = x // x has zeroes in position of Brake variable
	}
	ivb = dstream.Linapply(ivb, dirs, "dr")
	ivb.Reset()
	md := dstream.GetCol(ivb, "dr0").([]float64)
	ivb.Reset()
	cd1 := dstream.GetCol(ivb, "dr1").([]float64)
	ivb.Reset()
	cd2 := dstream.GetCol(ivb, "dr2").([]float64)
	ivb.Reset()
	uy := dstream.GetCol(ivb, "Brake").([]float64)
	bins1, bins2, counts := hist2d(md, cd1, uy, 50)
	histFile, err := os.Create("smhist_018.txt") //os.Create("/scratch/stats_flux/luers/hist_multi.txt")
	if err != nil {
		panic(err)
	}

	pFile, err := os.Create("smproj_multi.txt")
	if err != nil {
	   panic(err)
	}

	wCsvProj := csv.NewWriter(pFile)
	defer pFile.Close()

	prec := make([]string, 4)
	if err := wCsvProj.Write([]string{"meandir", "cd1", "cd2", "y"}); err != nil {
	   panic(err)
	}
	for i := 0; i < len(md); i++{
	    prec[0] = fmt.Sprintf("%v", md[i])
	    prec[1] = fmt.Sprintf("%v", cd1[i])
	    prec[2] = fmt.Sprintf("%v", cd2[i])
    	    prec[3] = fmt.Sprintf("%v", uy[i])
	    if err := wCsvProj.Write(prec); err != nil {
	       panic(err)
	    }
	}
	wCsvProj.Flush()
	
	wCsv := csv.NewWriter(histFile)
	defer histFile.Close()
	rec := make([]string, 4)
	if err := wCsv.Write([]string{"bin_meandir", "bin_covdir", "num0", "num1"}); err != nil {
		panic(err)
	}
	for i := 0; i < len(bins1); i++ {
		rec[0] = fmt.Sprintf("%v", bins1[i])
		rec[1] = fmt.Sprintf("%v", bins2[i])
		rec[2] = fmt.Sprintf("%v", counts[i][0])
		rec[3] = fmt.Sprintf("%v", counts[i][1])
		if err := wCsv.Write(rec); err != nil {
			panic(err)
		}
	}
	wCsv.Flush()

	
}
