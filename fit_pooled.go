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

const (
	maxLagLarge int     = 30 //30 samples = 30 * 100 milliseconds = 3 seconds
	maxLagSmall int     = 10
)

func main() {
	maxDriverID := 2
	lagmap := map[string]int{"Speed": maxLagLarge, "FcwRange": maxLagLarge, 
	       	  	         "FcwRangeRate": maxLagSmall,
				 "Steer": maxLagSmall}
	floatvars1 := []string{"Driver", "Trip", "Time", "Speed", 
		      		"Brake","FcwValidTarget", "FcwRange",
				"FcwRangeRate", "Steer"}

	start := time.Now()
	ivb := loadPoolDat(maxDriverID, floatvars1)	
	ivb = doTransforms(ivb, lagmap)
	fmt.Printf("names after doTransforms: %v\n", ivb.Names())
	ivb = laggedInteraction(ivb, "Speed", "FcwRange", maxLagLarge)
	fmt.Printf("Finished data transformations\n")
	
	ivb.Reset()
	ivb = dstream.Regroup(ivb, "Driver", false)

	// fmt.Printf("Variable names after transformations: %v\n", ivb.Names())

	// ---------- Fitting DOC -------------

	ivr0 := dstream.DropCols(ivb, []string{"Driver", "Brake", "Trip", "Time"})

	npc := 0
	doc0 := dimred.NewDOC(ivr0, respvar) // .SetProjection(npc)
	doc0.SetLogFile(fmt.Sprintf("log_pooled_%dpc.txt", npc))
	doc0.Done()
	ndir := 10
	doc0.Fit(ndir) // fit DOC without any PC projections

	// ---- Save the directions from DOC  ----
	var dirnames []string
	dirnames = append(dirnames, "varname", "doc0")
	dirs0 := make([][]float64, 1+ndir)
	dirs0[0] = doc0.MeanDir()
	for j := 0; j < ndir; j++ {
		dirs0[1+j] = doc0.CovDir(j)
		dirnames = append(dirnames, fmt.Sprintf("doc%d", j+1))
	}

	pdim := doc0.Dim() 
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
		 //  ("/scratch/stats_flux/luers/directions_%dpc.txt", npc))
		     //

	if err != nil {
		panic(err)
	}
	wDir := csv.NewWriter(dirFile)
	if err := wDir.Write(dirnames); err != nil {
		panic(err)
	}

	xnames := make([]string, ivr0.NumVar() - 1) // names of the "x" variables

	vm := make(map[string]int) // map from variable names to column positions

	// create map of variable names using original ivb dstream,
	// which contains Driver ids
	for k, na := range ivb.Names() { 
	    vm[na] = k
	}

	// store the names of the non-response variables
	ct := 0
	for _, na := range ivr0.Names() {
	    if na != respvar {
	       xnames[ct] = na
	       ct++
	    }
	}

	drec := make([]string, len(dirnames))
	for k, na := range xnames {
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
	
	// project against DOC directions
	dirs_expand := make([][]float64, len(dirs0))
	for j := 0; j < len(dirs0); j++ {
	    x := make([]float64, len(vm)) // number variables in the original dstream
	    for k, na := range xnames {
	    	x[vm[na]] = dirs0[j][k] // zeroes in the non-x-variable positions
	    }
	    dirs_expand[j] = x
	}

	ivb.Reset()
	ivb = dstream.Linapply(ivb, dirs_expand, "doc")
	fmt.Printf("\nFinished Linapply on ivb\n")
	fmt.Printf("ivr0.Names() = %v\n", ivr0.Names())
	fmt.Printf("ivb.Names() = %v\n", ivb.Names())
	start = time.Now()
	lagdat_fname := "data/lagdat_small_"
		     // "/scratch/stats_flux/luers/lagdat_" 
	ivbproj := dstream.DropCols(ivb, xnames)
	proj_fname := fmt.Sprintf("data/projdat_small_%dpc_", npc)
		   // ("/scratch/stats_flux/luers/projdat_%dpc_", npc)
	werr := dstream.ToCSV(ivb).DoneByChunk("Driver", "%03d", lagdat_fname, ".txt")
	if werr != nil {
		fmt.Printf("failed to write transformed data to disk\n")
	}
	werr = dstream.ToCSV(ivbproj).DoneByChunk("Driver", "%03d", proj_fname, ".txt")
	if werr != nil {
	   	fmt.Printf("failed to write projected scores to disk\n")
	}
	fmt.Printf("\nWrote transformed and doc-projected data to disk.\nElapsed time: %v minutes\n\n", time.Since(start).Minutes())

}

