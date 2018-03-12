package main

import (
	"encoding/csv"
	"fmt"
	"github.com/brookluers/dimred"
	"github.com/brookluers/dstream/dstream"
	"os"
)

func main() {
	maxDriverID := 2
	
	ndoc := 8
	projcols := make([]string, ndoc + 1)
	for j :=0; j <= ndoc; j++{
	    projcols[j] = fmt.Sprintf("doc%d", j)
	}
 	npc := 0
	ivproj := loadProjDat(maxDriverID, "data/projdat_0pc_%03d.txt", append(projcols, "Driver", "Brake_1sec"))
	respvar_lda := "Driver_minus1"
	ivproj = dstream.Apply(ivproj, respvar_lda, idm1, "float64")
	ivproj_nobr := dstream.DropCols(ivproj, []string{"Driver", "Brake_1sec"})
	fmt.Printf("ivproj.Names() = %v\n", ivproj.Names())
	fmt.Printf("ivproj_nobr.Names() = %v\n", ivproj_nobr.Names())
	lfit := dimred.NewLDA(ivproj_nobr, respvar_lda, maxDriverID)
	lfit.SetLogFile(fmt.Sprintf("ldalog_%d_drivers.txt", maxDriverID))
	lfit.Done()
	lfit.Fit()
		
	ivproj.Reset()
	ivproj = dstream.Convert(ivproj, "Driver", "uint64")
	ivproj = dstream.Regroup(ivproj, "Driver", false)

	// ---- Save the directions from LDA  ----
	ndir := lfit.GetNDir()
	dirnames := make([]string, ndir + 1)
	dirnames[0] = "varname"

	dirs0 := make([][]float64, ndir)
	for j := 0; j < ndir; j++ {
		dirs0[j] = lfit.LDir(j)
		dirnames[j+1] = fmt.Sprintf("lda%d", j)
	}

	dirFile, err := os.Create(fmt.Sprintf("data/lda_coord_%d_drivers_%dpc.txt", maxDriverID, npc))
		     //"/scratch/stats_flux/luers/directions_%dpc.txt", npc)
	if err != nil {
		panic(err)
	}
	wDir := csv.NewWriter(dirFile)
	if err := wDir.Write(dirnames); err != nil {
		panic(err)
	}

	xnames := make([]string, ivproj_nobr.NumVar() - 1) // names of the "x" variables

	vm := make(map[string]int) // map from variable names to column positions

	// create map of variable names using original ivb dstream,
	// which contains Driver ids
	for k, na := range ivproj.Names() { 
	    vm[na] = k
	}

	// store the names of the non-response variables
	ct := 0
	for _, na := range ivproj_nobr.Names() {
	    if na != respvar_lda {
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
	
	// project against LDA directions
	dirs_expand := make([][]float64, len(dirs0))
	for j := 0; j < len(dirs0); j++ {
	    x := make([]float64, len(vm)) // number variables in the original dstream
	    for k, na := range xnames {
	    	x[vm[na]] = dirs0[j][k] // zeroes in the non-x-variable positions
	    }
	    dirs_expand[j] = x
	}

	ivproj.Reset()
	ivproj = dstream.Linapply(ivproj, dirs_expand, "lda")
	fmt.Printf("\nFinished Linapply on ivproj\n")
	fmt.Printf("ivproj_nobr.Names() = %v\n", ivproj_nobr.Names())
	fmt.Printf("ivproj.Names() = %v\n", ivproj.Names())

	ldproj_fname :=  fmt.Sprintf("data/ldproj_%dpc_", npc)
		     //"/scratch/stats_flux/luers/lagdat_"
	werr := dstream.ToCSV(ivproj).DoneByChunk("Driver", "%03d", ldproj_fname, ".txt")
	if werr != nil {
		panic(werr)
	}

}

