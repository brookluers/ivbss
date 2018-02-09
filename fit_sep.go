package main

import (
	"encoding/csv"
	"fmt"
	"bufio"
	"github.com/brookluers/dimred"
	"github.com/brookluers/dstream/dstream"
	"github.com/gonum/matrix/mat64"
	"math"
	"os"
	"sort"
)

func main() {
     	   
	maxDriverID := 3

	var regxnames []string
	for j := maxSpeedLag; j >= 0; j-- {
		regxnames = append(regxnames, fmt.Sprintf("Speed[%d]", -j))
	}
	for j := maxRangeLag; j >= 0; j-- {
		regxnames = append(regxnames, fmt.Sprintf("FcwRange[%d]", -j))
	}

	pdim := len(regxnames)
	npc := 8
	ndir := 10

	for i := 1; i <= maxDriverID; i++ {
		ivb := loadDriverDat(i)
		ivb = doTransforms(ivb)

		// ---------- Fitting DOC -------------

		ivr0 :=  dstream.DropCols(ivb, []string{"Driver", "Brake", "Trip", "Time"})

		fmt.Printf("Variable names for DOC: %v\n", ivr0.Names())
		fmt.Printf("Variable names for original dstream: %v\n", ivb.Names())
		doc0 := dimred.NewDOC(ivr0, "Brake_1sec")
		doc0.SetLogFile(fmt.Sprintf("log%03d.txt", i))
		doc0.Done()
		doc0.Fit(ndir)

		covfile, err := os.Create(fmt.Sprintf("data/sep_cov_small_%03d.txt", i))
			 //os.Create(fmt.Sprintf("/scratch/stats_flux/luers/sep_cov_%03d.txt", i))
		if err != nil {
		   panic(err)
		}
		wCov := bufio.NewWriter(covfile)
		
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

		//  fmt.Printf("-----Marginal Covariance-----\n%v\n", margcov)
		margcorr := mat64.NewSymDense(pdim, t3)
		//  fmt.Printf("-----Marginal Correlation-----\n%v\n", margcorr)
		es := new(mat64.EigenSym)
		ok := es.Factorize(margcorr, true) //margcov, true)
		if !ok {
			panic("can't factorize marginal correlation matrix\n")
		}
		marg_evals := es.Values(nil)
		sort.Float64s(marg_evals)
		_, err = fmt.Fprintf(wCov, "%v\n", regxnames)
		if err != nil {
		   panic(err)
		}		
		_, err = fmt.Fprintf(wCov, "marginal eigenvalues\n%v\n", marg_evals)
		if err != nil {
		   panic(err)
		}
		fmt_margcov := mat64.Formatted(margcov, mat64.Prefix("\n     "), mat64.Squeeze())
		_, err = fmt.Fprintf(wCov, "marginal covariance\n%v\n", fmt_margcov)
		if err != nil {
		   panic(err)
		}
		_, err = fmt.Fprintf(wCov, "DOC eigenvalues\n%v\n", doc0.Eig())
		if err != nil {
		   panic(err)
		}
		cov1 := mat64.NewSymDense(pdim, doc0.GetCov(1))
		cov0 := mat64.NewSymDense(pdim, doc0.GetCov(0))
		fmt_cov1 := mat64.Formatted(cov1, mat64.Prefix("     "), mat64.Squeeze())
		fmt_cov0 := mat64.Formatted(cov0, mat64.Prefix("     "), mat64.Squeeze())
		_, err = fmt.Fprintf(wCov, "covariance, y=1\n%v\n", fmt_cov1)
		if err != nil {
		   panic(err)
		}
		_, err = fmt.Fprintf(wCov, "covariance, y=0\n%v\n", fmt_cov0)
		if err != nil {
		   panic(err)
		}
		wCov.Flush()		
		
		evec := new(mat64.Dense)
		evec.EigenvectorsSym(es)
		dirFile, err := os.Create(fmt.Sprintf("data/dir_small_sep_%03d.txt", i))
			 //os.Create(fmt.Sprintf("/scratch/stats_flux/luers/dir_sep_%03d.txt", i))
		if err != nil {
			panic(err)
		}
		wDir := csv.NewWriter(dirFile)
		var temp []string
		temp = append(temp, "varname", "meandir")

		dirs0 := make([][]float64, 1+ndir+npc)
		dirs0[0] = doc0.MeanDir()
		for j := 0; j < ndir; j++ {
			//dirs0[1 + j] = make([]float64, pdim)
			dirs0[1+j] = doc0.CovDir(j)
			temp = append(temp, fmt.Sprintf("cd%d", j+1))
		}

		pcMat := evec.View(0, pdim-npc, pdim, npc) // PCs can be applied to standardized x
		pcMatDense := mat64.DenseCopyOf(pcMat)
		pcMatDense.Mul(sd_inv_Diag, pcMatDense) // these directions can be applied to raw x
		// eigenvalues sorted in increasing order
		for j := 0; j < npc; j++ {
			dirs0[1+ndir+j] = mat64.Col(nil, npc - 1 - j, pcMatDense)
			temp = append(temp, fmt.Sprintf("pc%d", j))
		}

		// temp = {"varname", "meandir", "cd1"..."pc1"...}
		// len(temp): 2 + ndir + npc
		// dirs0: mean direction, covariance directions, pc directions
		if err := wDir.Write(temp); err != nil {
			panic(err)
		}

		// save coefficient vectors for each dimension reduction direction
		drec := make([]string, len(temp))
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

}
