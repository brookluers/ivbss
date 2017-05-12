package main

import (
	"compress/gzip"
	"encoding/json"
	"fmt"
	"math"
	"os"

	"github.com/gonum/floats"
	"github.com/kshedden/dimred"
	"github.com/kshedden/statmodel/dataprovider"
)

const (
	maxSpeedLag int = 30
	maxRangeLag int = 30
	minCount    int = 100
	nrow        int = 100
	ncol        int = 100
)

func selectEq(w float64) dataprovider.FilterColFunc {
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

func getScores(data dataprovider.Data, j0, j1, jb int) ([]float64, []float64, []float64) {

	data.Reset()
	var x0, x1, y []float64
	for data.Next() {
		z0 := data.Get(j0).([]float64)
		z1 := data.Get(j1).([]float64)
		yy := data.Get(jb).([]float64)
		x0 = append(x0, z0...)
		x1 = append(x1, z1...)
		y = append(y, yy...)
	}

	mn0 := floats.Sum(x0) / float64(len(x0))
	mn1 := floats.Sum(x1) / float64(len(x1))

	floats.AddConst(-mn0, x0)
	floats.AddConst(-mn1, x1)

	var s0, s1 float64
	for i, _ := range x0 {
		s0 += x0[i] * x0[i]
		s1 += x1[i] * x1[i]
	}
	s0 = math.Sqrt(s0 / float64(len(x0)))
	s1 = math.Sqrt(s1 / float64(len(x1)))

	floats.Scale(1/s0, x0)
	floats.Scale(1/s1, x1)

	return y, x0, x1
}

func getPos(data dataprovider.Data, name string) int {
	for k, v := range data.Names() {
		if v == name {
			return k
		}
	}
	panic("cannot find " + name)
}

func cellMeans(data dataprovider.Data, row, col []int) ([][]float64, []int) {

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
			v := data.Get(k).([]float64)
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

func main() {

	f0 := selectEq(0)
	f1 := selectEq(1)
	f10 := selectEq(10)

	rdr, err := os.Open("/nfs/turbo/ivbss/LvFot/data_001.txt")
	if err != nil {
		panic(err)
	}
	ivx := dataprovider.NewStreamFromCSV(rdr)
	ivx.SetFloatVars([]string{"Trip", "Time", "Speed", "Brake", "FcwValidTarget", "FcwRange"})
	ivx.SetChunkSize(5000)

	ivb := dataprovider.Segment(ivx, []string{"Trip"})

	ivb = dataprovider.Apply(ivb, "brake2", fbrake, "float64")
	ivb = dataprovider.FilterCol(ivb, map[string]dataprovider.FilterColFunc{
		"brake2": f0, "FcwValidTarget": f1})
	ivb = dataprovider.DiffChunk(ivb, map[string]int{"Time": 1})
	ivb = dataprovider.LagChunk(ivb, map[string]int{"Speed": maxSpeedLag, "FcwRange": maxRangeLag})

	ivb = dataprovider.Segment(ivb, []string{"Time$d1"})
	ivb = dataprovider.FilterCol(ivb, map[string]dataprovider.FilterColFunc{
		"Time$d1": f10})

	ivb = dataprovider.ChunksizeFilter(ivb, 1, 1e10)
	ivb = dataprovider.Drop(ivb, []string{"Trip", "Time", "Time$d1", "FcwValidTarget", "brake2"})

	var ivr dataprovider.Reg
	ivr = dataprovider.NewReg(ivb, "Brake", nil, "", "")

	doc := dimred.NewDOC(ivr)
	doc.SetLogFile("log.txt")
	doc.Init()
	doc.Fit(2)

	fmt.Printf("nobs after fit=%d\n", ivb.Nobs())

	dirs0 := [][]float64{doc.CovDir(0), doc.MeanDir()}
	//dirs0 := [][]float64{doc.CovDir(0), doc.CovDir(1)}

	margcov := doc.MargCov()
	standardize(dirs0[0], margcov)
	standardize(dirs0[1], margcov)

	// Expand to match the data set
	vm := make(map[string]int)
	dirs := make([][]float64, 2)
	for k, a := range ivb.Names() {
		vm[a] = k
	}
	for j := 0; j < 2; j++ {
		x := make([]float64, len(vm))
		for k, na := range ivr.XNames() {
			x[vm[na]] = dirs0[j][k]
		}
		dirs[j] = x
	}
	ivb = dataprovider.Linapply(ivb, dirs, "dr")

	j0 := getPos(ivb, "dr0")
	j1 := getPos(ivb, "dr1")
	jb := getPos(ivb, "Brake")

	y, x0, x1 := getScores(ivb, j0, j1, jb)
	row, col := getCells(x0, x1)
	cmn, cmc := cellMeans(ivb, row, col)
	standardizeCellMeans(cmn, cmc)
	rat, _ := heatMap(row, col, y, x0, x1)

	// Export the heatmap data
	fid, err := os.Create("heat1.json.gz")
	if err != nil {
		panic(err)
	}
	gzw := gzip.NewWriter(fid)
	enc := json.NewEncoder(gzw)
	for i, _ := range rat {
		if rat[i] != -1 {
			x := &struct {
				Row int
				Col int
				Z   float64
			}{
				Row: i / ncol,
				Col: i % ncol,
				Z:   rat[i],
			}
			enc.Encode(x)
		}
	}
	gzw.Close()
	fid.Close()

	// Export the trajectory data
	fid, err = os.Create("trajectory1.json.gz")
	if err != nil {
		panic(err)
	}
	gzw = gzip.NewWriter(fid)
	enc = json.NewEncoder(gzw)
	for i, v := range cmn {
		if rat[i] != -1 {
			x := &struct {
				Row int
				Col int
				Tr  []float64
			}{
				Row: i / ncol,
				Col: i % ncol,
				Tr:  v,
			}
			enc.Encode(x)
		}
	}
	gzw.Close()
	fid.Close()

	fmt.Printf("CovDir(0) = %v\n", doc.CovDir(0))
}
