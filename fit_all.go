package main

import (
	"fmt"
	"log"
	"math"
	"sort"

	"github.com/brookluers/dimred"
	"github.com/brookluers/dstream/dstream"
	"github.com/gonum/floats"
	"github.com/gonum/plot"
	"github.com/gonum/plot/palette"
	"github.com/gonum/plot/plotter"
	"github.com/gonum/plot/plotutil"
	"github.com/gonum/plot/vg"
)

const (
	maxSpeedLag int     = 30
	maxRangeLag int     = 30
	minCount    int     = 100
	nrow        int     = 100
	ncol        int     = 100
	minSpeed    float64 = 10.0
)

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

//pmaxFunc returns an ApplyFunc that
// computes the maximum "horizontally" across
// a set of variables defined by vnames (variable names)
func pmaxFunc(vnames []string) dstream.ApplyFunc {
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
				if data[j][i] > y[i] {
					y[i] = data[j][i]
				}
			}
		}

	}
	return f
}

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
				z[j]++
			}
		}
	}

	for i, _ := range z {
		z[i] /= float64(2 * w)
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

func cellMeans(data dstream.Dstream, row, col []int) ([][]float64, []int) {

	var il []int
	for k := 0; k <= maxSpeedLag; k++ {
		il = append(il, getPos(data, fmt.Sprintf("Speed[%d]", -k)))
	}
	for k := 0; k <= maxRangeLag; k++ {
		il = append(il, getPos(data, fmt.Sprintf("FcwRange[%d]", -k)))
	}

	cmn := make([][]float64, nrow*ncol)
	cmc := make([]int, nrow*ncol)
	for i := 0; i < nrow*ncol; i++ {
		cmn[i] = make([]float64, len(il))
	}

	data.Reset()
	ii := 0
	for data.Next() {
		var n int
		for j, k := range il {
			v := data.GetPos(k).([]float64)
			n = len(v)
			for i := 0; i < len(v); i++ {
				jj := ii + i
				if row[jj] >= 0 && row[jj] < nrow && col[jj] >= 0 && col[jj] < ncol {
					q := row[jj]*ncol + col[jj]
					cmn[q][j] += v[i]
					if j == 0 {
						cmc[q]++
					}
				}
			}
		}
		ii += n
	}

	for q, v := range cmn {
		floats.Scale(1/float64(cmc[q]), v)
	}

	return cmn, cmc
}

func standardizeCellMeans(cmn [][]float64, cmc []int) {

	v := make([]float64, len(cmn[0]))

	// Center
	w := 0
	for i, u := range cmn {
		if cmc[i] > 0 {
			floats.AddScaled(v, float64(cmc[i]), u)
			w += cmc[i]
		}
	}
	floats.Scale(1/float64(w), v)
	for _, u := range cmn {
		floats.Sub(u, v)
	}

	// Scale
	for j, _ := range v {
		v[j] = 0
	}
	w = 0
	for i, u := range cmn {
		if cmc[i] > 0 {
			for j, _ := range u {
				v[j] += float64(cmc[i]) * u[j] * u[j]
			}
			w += cmc[i]
		}
	}
	floats.Scale(1/float64(w), v)
	for j, x := range v {
		v[j] = math.Sqrt(x)
	}
	for i, u := range cmn {
		for j, x := range u {
			cmn[i][j] = x / v[j]
		}
	}
}

func heatMap(row, col []int, y, x0, x1 []float64) ([]float64, []int) {
	missed := 0
	hit := 0
	num := make([]float64, nrow*ncol)
	denom := make([]int, nrow*ncol)
	for i, _ := range x0 {
		if row[i] >= 0 && col[i] >= 0 && row[i] < nrow && col[i] < ncol {
			denom[ncol*row[i]+col[i]]++
			if y[i] == 1 {
				num[ncol*row[i]+col[i]]++
			}
			hit++
		} else {
			missed++
		}
	}
	fmt.Printf("Missed %d\n", missed)
	fmt.Printf("Hit %d\n", hit)

	rat := make([]float64, nrow*ncol)
	for i, _ := range num {
		if denom[i] > minCount {
			rat[i] = math.Pow(num[i]/float64(denom[i]), 0.1)
		} else {
			rat[i] = -1
		}
	}

	return rat, denom
}

func getCells(x0, x1 []float64) ([]int, []int) {
	row := make([]int, len(x0))
	col := make([]int, len(x1))
	for i, _ := range x0 {
		row[i] = int(math.Floor(70*x0[i] + 50))
		col[i] = int(math.Floor(15*x1[i] + 50))
	}
	return row, col
}

// Generate a heatmap of a m x m covariance matrix, converting it to a
// correlation matrix if scale is true.
func plotcov(cov []float64, scale bool, m int, pc plotconfig, fname string) {

	pal := palette.Heat(100, 1)
	da := &covheat{data: cov, m: m, scale: scale}
	h := plotter.NewHeatMap(da, pal)

	p, err := plot.New()
	if err != nil {
		log.Panic(err)
	}
	p.Title.Text = pc.title
	p.X.Label.Text = pc.xlabel
	p.Y.Label.Text = pc.ylabel
	p.Add(h)

	if err := p.Save(4*vg.Inch, 4*vg.Inch, fname); err != nil {
		panic(err)
	}
}

// Configuration parameters for a plot.
type plotconfig struct {
	title  string
	xlabel string
	ylabel string
}

// Plot one or more lines as functions on a plot.  If scale is true,
// scale the data for each line to have unit L2 norm.
func plotlines(x [][]float64, scale bool, names []string, pc plotconfig, fname string) {

	gxy := func(x []float64) plotter.XYs {

		s := float64(0)
		if scale {
			for _, v := range x {
				s += v * v
			}
			s = math.Sqrt(s)
		} else {
			s = 1
		}

		z := make(plotter.XYs, len(x))
		for i, y := range x {
			z[i].X = float64(i)
			z[i].Y = y / s
		}
		return z
	}

	p, err := plot.New()
	if err != nil {
		panic(err)
	}

	p.Title.Text = pc.title
	p.X.Label.Text = pc.xlabel
	p.Y.Label.Text = pc.ylabel

	var ag []interface{}
	for i, z := range x {
		ag = append(ag, names[i])
		ag = append(ag, gxy(z))
	}
	err = plotutil.AddLinePoints(p, ag...)
	if err != nil {
		panic(err)
	}

	if err := p.Save(4*vg.Inch, 4*vg.Inch, fname); err != nil {
		panic(err)
	}
}

func plotscatter(x []float64, y []float64, pc plotconfig, fname string) {

	z := make(plotter.XYs, len(x))
	for i := range x {
		z[i].X = x[i]
		z[i].Y = y[i]
	}

	p, err := plot.New()
	if err != nil {
		panic(err)
	}

	s, err := plotter.NewScatter(z)
	if err != nil {
		panic(err)
	}

	p.Title.Text = pc.title
	p.X.Label.Text = pc.xlabel
	p.Y.Label.Text = pc.ylabel
	p.Add(s)

	if err := p.Save(4*vg.Inch, 4*vg.Inch, fname); err != nil {
		panic(err)
	}
}

// Implement the XYZGrid interface for making a heatmap.
type covheat struct {
	m     int       // The data are an m x m array
	data  []float64 // The data
	scale bool      // If true, scale as a correlation matrix
}

func (h *covheat) Dims() (int, int) {
	return h.m, h.m
}

func (h *covheat) Z(r, c int) float64 {
	if h.scale {
		return h.data[r*h.m+c] / math.Sqrt(h.data[r*h.m+r]*h.data[c*h.m+c])
	} else {
		return h.data[r*h.m+c]
	}
}

func (h *covheat) X(c int) float64 {
	return float64(c)
}

func (h *covheat) Y(r int) float64 {
	return float64(h.m) - float64(r)
}

func main() {

	f0 := selectEq(0)
	f1 := selectEq(1)
	f10 := selectEq(10)

	maxDriverID := 108
	fnames := make([]string, maxDriverID)
	fmt.Println("File names:")
	for i := 1; i <= maxDriverID; i++ {
		fnames[i-1] = fmt.Sprintf("/nfs/turbo/ivbss/LvFot/data_%03d.txt", i)
		fmt.Println(fnames[i-1])
	}

	mrdr := NewMultiReadSeek0(fnames, true) // every file has a header row

	ivx := dstream.FromCSV(mrdr).SetFloatVars([]string{"Trip", "Time", "Speed", "Brake", "FcwValidTarget", "FcwRange"}).SetChunkSize(10000).HasHeader()

	ivb := dstream.Apply(ivx, "brake2", fbrake, "float64")
	ivb = dstream.Filter(ivb, map[string]dstream.FilterFunc{"brake2": f0, "FcwValidTarget": f1})
	ivb = dstream.DiffChunk(ivb, map[string]int{"Time": 1})
	ivb = dstream.LagChunk(ivb, map[string]int{"Speed": maxSpeedLag, "FcwRange": maxRangeLag})
	ivb = dstream.Segment(ivb, []string{"Trip", "Time$d1"})
	ivb = dstream.Filter(ivb, map[string]dstream.FilterFunc{"Time$d1": f10})
	ivb = dstream.DropCols(ivb, []string{"Trip", "Time", "Time$d1", "FcwValidTarget", "brake2"})

	//	// Plot the distribution of block sizes
	//	var bsize []float64
	//	ivb.Reset()
	//	for ivb.Next() {
	//		n := len(ivb.GetPos(0).([]float64))
	//		if n > 0 {
	//			bsize = append(bsize, math.Log(float64(n))/math.Log(10))
	//		}
	//	}
	//	sort.Sort(sort.Float64Slice(bsize))
	//	plotlines([][]float64{bsize}, false, []string{"Size"}, plotconfig{ylabel: "Log10 block size"}, "blocksizes.pdf")

	//compute the minimum speed within the 3-second window
	spNames := make([]string, maxSpeedLag+1)
	for i := 0; i <= maxSpeedLag; i++ {
		spNames[i] = fmt.Sprintf("%s[%d]", "Speed", -i)
	}
	pmin := pminFunc(spNames)
	//pmax := pmaxFunc(spNames)
	ivb = dstream.Apply(ivb, "pminSpeed", pmin, "float64")
//	ivb = dstream.Apply(ivb, "pmaxSpeed", pmax, "float64")

	//keep windows where the minimum speed is at least minSpeed
	fmt.Printf("\nOnly retain windows where the minimum speed is at least %v\n", minSpeed)

	spOrBrake := func(v map[string]interface{}, ret interface{}) {
		//		br := v["Brake"].([]float64)
				msp := v["pminSpeed"].([]float64)
		//		masp := v["pmaxSpeed"].([]float64)

		tf := ret.([]float64)
		for ix := 0; ix < len(tf); ix++ {
			if msp[ix] > minSpeed {
				tf[ix] = 1.0
			} else {
				tf[ix] = 0.0
			}
		}
	}
	ivb = dstream.Apply(ivb, "speedThreshKeep", spOrBrake, "float64")
	ivb = dstream.Filter(ivb, map[string]dstream.FilterFunc{"speedThreshKeep": f1})

	fmt.Printf("Variable names after transformations: %v\n", ivb.Names())
	ivb = dstream.DropCols(ivb, []string{"speedThreshKeep", "pminSpeed"})

	// ---------- Fitting DOC -------------
	ivr := dstream.NewReg(ivb, "Brake", nil, "", "")

	doc := dimred.NewDOC(ivr)
	doc.SetLogFile("log.txt")
	doc.Init()

	ndir := 2
	doc.Fit(ndir)

	fmt.Printf("nobs after fit=%d\n", ivb.NumObs())
	fmt.Printf("%v\n", ivr.XNames()[0:31])
	fmt.Printf("%v\n", ivr.XNames()[31:62])

	fmt.Printf("%d %d %d\n", len(doc.YMean(0)), len(doc.MeanDir()), len(doc.CovDir(0)))

	z := [][]float64{doc.YMean(0)[0:31], doc.YMean(1)[0:31]}
	plotlines(z, false, []string{"0", "1"}, plotconfig{title: "Mean speed", xlabel: "Time lag", ylabel: "Speed"}, "meanspeed.pdf")
	z = [][]float64{doc.YMean(0)[31:62], doc.YMean(1)[31:62]}
	plotlines(z, false, []string{"0", "1"}, plotconfig{title: "Mean range", xlabel: "Time lag", ylabel: "Range"}, "meanrange_alldrivers.pdf")

	plotcov(doc.YCov(0), true, 62, plotconfig{title: "Non-braking correlation", xlabel: "Time lag", ylabel: "Time lag"}, "cov0.pdf")
	plotcov(doc.YCov(1), true, 62, plotconfig{title: "Braking correlation", xlabel: "Time lag", ylabel: "Time lag"}, "cov1_alldrivers.pdf")

	covdiff := make([]float64, 62*62)
	floats.SubTo(covdiff, doc.YCov(1), doc.YCov(0))
	plotcov(covdiff, false, 62, plotconfig{title: "Covariance difference", xlabel: "Time lag", ylabel: "Time lag"}, "covdiff_alldrivers.pdf")

	z = [][]float64{doc.MeanDir()[0:31], doc.MeanDir()[31:62]}
	plotlines(z, true, []string{"Speed", "Range"}, plotconfig{xlabel: "Time lag", ylabel: "Coefficient"}, "mean_dir_alldrivers.pdf")

	z = [][]float64{doc.CovDir(0)[0:31], doc.CovDir(1)[0:31]}
	plotlines(z, true, []string{"Cov1", "Cov2"}, plotconfig{xlabel: "Time lag", ylabel: "Speed"}, "speed_dir_alldrivers.pdf")
	z = [][]float64{doc.CovDir(0)[31:62], doc.CovDir(1)[31:62]}
	plotlines(z, true, []string{"Cov1", "Cov2"}, plotconfig{xlabel: "Time lag", ylabel: "Range"}, "range_dir_alldrivers.pdf")

	dirs0 := [][]float64{doc.MeanDir(), doc.CovDir(0), doc.CovDir(1)}

	// Expand to match the data set
	vm := make(map[string]int)
	dirs := make([][]float64, 3)
	for k, a := range ivb.Names() {
		vm[a] = k
	}
	for j := 0; j < 3; j++ {
		x := make([]float64, len(vm))
		for k, na := range ivr.XNames() {
			x[vm[na]] = dirs0[j][k]
		}
		dirs[j] = x
	}
	ivb = dstream.Linapply(ivb, dirs, "dr")

	ww := 3000
	ux := dstream.GetCol(ivb, "dr0").([]float64)
	uy := dstream.GetCol(ivb, "Brake").([]float64)
	x0, b0 := getScores(ux, uy, ww)
	plotscatter(x0, b0, plotconfig{xlabel: "Mean direction", ylabel: "P(Brake)"}, "dr0_alldrivers.png")

	ux = dstream.GetCol(ivb, "dr1").([]float64)
	uy = dstream.GetCol(ivb, "Brake").([]float64)
	x0, b0 = getScores(ux, uy, ww)
	plotscatter(x0, b0, plotconfig{xlabel: "Covariance direction 1", ylabel: "P(Brake)"}, "dr1_alldrivers.png")

	ux = dstream.GetCol(ivb, "dr2").([]float64)
	uy = dstream.GetCol(ivb, "Brake").([]float64)
	x0, b0 = getScores(ux, uy, ww)
	plotscatter(x0, b0, plotconfig{xlabel: "Covariance direction 1", ylabel: "P(Brake)"}, "dr2_alldrivers.png")
}
