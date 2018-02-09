package main

import (
	"encoding/csv"
	"fmt"
	"github.com/brookluers/dimred"
	"github.com/brookluers/dstream/dstream"
	"github.com/gonum/matrix/mat64"
	"math"
	"os"
	//"sort"
	"time"
)


func main() {
	maxDriverID := 3
	
	start := time.Now()
	ivb := loadPoolDat(maxDriverID)	
	ivb = doTransforms(ivb)
	fmt.Printf("Finished data transformations\n")

	ivb.Reset()
	fmt.Printf("---Writing transformed data to disk  ---\n")	
	ivb = dstream.Regroup(ivb, "Driver", false)
	lagdat_fname :=  "data/lagdat_small_"
		     //"/scratch/stats_flux/luers/lagdat_"
	werr := dstream.ToCSV(ivb).DoneByChunk("Driver", "%03d", lagdat_fname, ".txt")
	if werr != nil {
		panic(werr)
	}
	fmt.Printf("Wrote transformed data to disk. Elapsed time: %v minutes\n\n", time.Since(start).Minutes())
	fmt.Printf("Variable names after transformations: %v\n", ivb.Names())

	// ---------- Fitting DOC -------------
	var regxnames []string
	for j := maxSpeedLag; j >= 0; j-- {
		regxnames = append(regxnames, fmt.Sprintf("Speed[%d]", -j))
	}
	for j := maxRangeLag; j >= 0; j-- {
		regxnames = append(regxnames, fmt.Sprintf("FcwRange[%d]", -j))
	}

	ivr0 := dstream.DropCols(ivb, []string{"Driver", "Brake", "Trip", "Time"})

	npc := 10
	doc0 := dimred.NewDOC(ivr0, "Brake_1sec").SetProjection(npc)
	doc0.SetLogFile("log_pooled.txt")
	doc0.Done()
	ndir := 10
	doc0.Fit(ndir) // fit DOC without any PC projections

	// ---- Save the directions from DOC without any PC projections ---
	var dirnames []string
	dirnames = append(dirnames, "varname", "meandir")
	dirs0 := make([][]float64, 1+ndir)
	dirs0[0] = doc0.MeanDir()
	for j := 0; j < ndir; j++ {
		dirs0[1+j] = doc0.CovDir(j)
		dirnames = append(dirnames, fmt.Sprintf("cd%d", j+1))
	}

	pdim := doc0.Dim() //len(regxnames)
	fmt.Printf("pdim = %v\n", pdim)
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

	dirFile, err := os.Create(fmt.Sprintf("data/directions_small_%dpc.txt", npc))
		     //"/scratch/stats_flux/luers/directions_%dpc.txt", npc)
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

}


