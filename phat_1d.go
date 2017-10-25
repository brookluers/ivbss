package main

import (
	"encoding/csv"
	"fmt"
	"github.com/gonum/floats"
	//"github.com/brookluers/dimred"
	"github.com/brookluers/dstream/dstream"
	//"github.com/gonum/matrix/mat64"
	"math"
	"os"
	"sort"
)

const (
	minSpeed  float64 = 7.0
	chunkSize int     = 25000
	ndir_phat int     = 2
)

// getBrakeProb estimates the probability of breaking at each value of
// a numeric score.  The breaking probabilities are estimated based on
// a local mean (+/- w positions from the score value being
// conditioned on).
func getBrakeProb(sc, br []float64, w int) ([]float64, []float64) {

	ii := make([]int, len(sc))
	floats.Argsort(sc, ii)

	// Reorder br to be compatible with sc.
	var b []float64
	for _, i := range ii {
		b = append(b, br[i])
	}

	z := make([]float64, len(sc))
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

	return sc[w : len(sc)-w], z[w : len(z)-w]
}

// getBrakeProb2 estimates the probability of breaking at each value of
// a numeric score.  The breaking probabilities are estimated based on
// a local mean (+/- w positions from the score value being
// conditioned on).
// compute the probability at np evenly spaced points
// in the range of sc
func getBrakeProb2(sc, br []float64, w, np int) ([]float64, []float64, []float64) {

	ii := make([]int, len(sc))
	floats.Argsort(sc, ii)
	
	// Reorder br to be compatible with sc.
	var b []float64
	for _, i := range ii {
		b = append(b, br[i])
	}
	
	evpoints := make([]float64, np)
	p_hat := make([]float64, np)
	ny := make([]float64, np)
	stderr := make([]float64, np)
	// sc is sorted
	var minev, maxev float64
	minev = sc[w]
	maxev = sc[len(sc) - w]
	var stepsize float64 = (maxev - minev) / float64(np)
	
	for i := 0; i < np; i++{
	    evpoints[i] = minev + float64(i) * stepsize
	    cix := sort.SearchFloat64s(sc, evpoints[i]) // location of evpoints[i] within sc
	    for j := cix - w; j < cix + w; j++{
	    	if j < 0 {
		   fmt.Printf("negative index = %v\n", j)
		   j = 0
		}
	    	ny[i] += b[j]
	    }
	    p_hat[i] = ny[i] / float64(2 * w)
	    stderr[i] = math.Sqrt((p_hat[i] * (1.0 - p_hat[i])) / float64(2 * w))
	}
	
	return evpoints, p_hat, stderr
}

func main() {
    fmt.Printf("\n\t--- Computing 1-dimensional conditional probability estimates---\t\n")
	maxDriverID := 108
	
	fvars := []string{"Driver", "Brake"}
	dirvars := []string{"meandir0"}
	for j := 0; j < ndir_phat; j++ {
		dirvars = append(dirvars, fmt.Sprintf("cd%d", j))
		dirvars = append(dirvars, fmt.Sprintf("pc%d", j))
	}
	fvars = append(fvars, dirvars...)
	fmt.Printf("fvars = %v\n", fvars)

	ww := 3000
	n_evals := 250
	fmt.Printf("window size = %v\n", ww)
	for f_ix := 1; f_ix <= maxDriverID; f_ix++ {
		fmt.Printf("Processing driver %d\n", f_ix)
		rdr, err := os.Open(fmt.Sprintf("/scratch/stats_flux/luers/smproj_8pc_%03d.txt", f_ix))
		if err != nil {
			panic(err)
		}

		projected := dstream.FromCSV(rdr).SetFloatVars(fvars).SetChunkSize(chunkSize).HasHeader().Done()
		phat := make([][]float64, len(dirvars))
		scores := make([][]float64, len(dirvars))
		stderr := make([][]float64, len(dirvars))
		dy := dstream.GetCol(projected, "Brake").([]float64)
		projected.Reset()
		for k := 0; k < len(dirvars); k++ {
			dx := dstream.GetCol(projected, dirvars[k]).([]float64)
			projected.Reset()
			scores[k], phat[k], stderr[k] = getBrakeProb2(dx, dy, ww, n_evals)

		}
		resfile, err := os.Create(fmt.Sprintf("/scratch/stats_flux/luers/phat_%03d.txt", f_ix))
		if err != nil {
			panic(err)
		}
		wcsv := csv.NewWriter(resfile)
		rec := []string{"Driver"}
		for _, v := range dirvars {
			rec = append(rec, v)
			rec = append(rec, v+"_phat")
			rec = append(rec, v+"_se")
		}
		//fmt.Printf("--Header---\nrec = %v\nlen(rec)=%v\n",rec,len(rec))
		if err := wcsv.Write(rec); err != nil {
			panic(err)
		}

		for i := 0; i < len(scores[0]); i++ {
			rec[0] = fmt.Sprintf("%v", f_ix)
			for k := 0; k < len(dirvars); k++ {
				rec[1+3*k] = fmt.Sprintf("%v", scores[k][i])
				rec[1+3*k+1] = fmt.Sprintf("%v", phat[k][i])
				rec[1+3*k+2] = fmt.Sprintf("%v", stderr[k][i])
			}
			if err := wcsv.Write(rec); err != nil {
				panic(err)
			}
		}
		wcsv.Flush()
		resfile.Close()
		rdr.Close()
	}
}
