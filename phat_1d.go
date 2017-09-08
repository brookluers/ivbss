package main

import (
	"encoding/csv"
	"fmt"
	"github.com/gonum/floats"
	//"github.com/brookluers/dimred"
	"github.com/brookluers/dstream/dstream"
	//"github.com/gonum/matrix/mat64"
	//"math"
	"os"
	//"sort"
)

const (
	maxSpeedLag int     = 30 //30 samples = 30 * 100 milliseconds = 3 seconds
	maxRangeLag int     = 30
	minCount    int     = 100
	nrow        int     = 100
	ncol        int     = 100
	minSpeed    float64 = 7.0
	chunkSize   int     = 10000
	ndir_phat   int	    = 5
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



func main() {
     maxDriverID := 2
     fnames := make([]string, maxDriverID)
     for i := 1; i <= maxDriverID; i++ {
		fnames[i-1] = fmt.Sprintf("smproj_multi_%03d.txt", i) 
		//fmt.Sprintf("/scratch/stats_flux/luers/smproj_multi_%03d.txt", i)
		fmt.Println(fnames[i-1])
     }
     fvars := []string{"Driver", "Brake"}
     dirvars := []string{"meandir0"}
     for j := 0; j < ndir_phat; j++{
     	 dirvars = append(dirvars, fmt.Sprintf("cd%d", j))
	 dirvars = append(dirvars, fmt.Sprintf("pc%d", j))
     }
     fvars = append(fvars, dirvars...)
     fmt.Printf("fvars = %v\n", fvars)
     
     ww := 3000
     for f_ix := 1; f_ix <= maxDriverID; f_ix++{
     	 fmt.Printf("Processing driver %d\n", f_ix)
     	 rdr, err := os.Open(fmt.Sprintf("smproj_multi_%03d.txt", f_ix))
	 if err != nil {
	    panic(err)
	 }

     	 projected := dstream.FromCSV(rdr).SetFloatVars(fvars).SetChunkSize(25000).HasHeader().Done()
	 phat := make([][]float64, len(dirvars))
	 scores := make([][]float64, len(dirvars))
	 dy := dstream.GetCol(projected, "Brake").([]float64)
	 projected.Reset()
	 for k := 0; k < len(dirvars); k++{
	     dx := dstream.GetCol(projected, dirvars[k]).([]float64)
	     projected.Reset()
	     fmt.Printf("len(dx) = %v\n", len(dx))
	     fmt.Printf("len(dy) = %v\n", len(dy))
	     scores[k], phat[k] = getBrakeProb(dx, dy, ww)
	 }
	 resfile, err := os.Create(fmt.Sprintf("phat_%03d.txt", f_ix))
	 if err != nil {
	    panic(err)
	 }
	 wcsv := csv.NewWriter(resfile)
	 rec := []string{"Driver"}
	 for _, v := range dirvars {
	     rec = append(rec, v)
	     rec = append(rec, v + "_phat")
	 }
	 fmt.Printf("--Header---\nrec = %v\nlen(rec)=%v\n",rec,len(rec))
	 if err := wcsv.Write(rec); err != nil {
	    panic(err)
	 }
	 
	 fmt.Printf("len(dirvars) = %v\n",len(dirvars))
	 for i := 0; i < len(scores[0]); i++{
	     rec[0] = fmt.Sprintf("%v", f_ix)
	     for k := 0; k < len(dirvars); k++ {
	     	 rec[1 + 2*k] = fmt.Sprintf("%v", scores[k][i])
		 rec[1 + 2*k + 1] = fmt.Sprintf("%v", phat[k][i])
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