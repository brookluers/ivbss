package main

import (
	"bufio"
	//"encoding/csv"
	"fmt"
	"github.com/brookluers/dimred"
	"github.com/brookluers/dstream/dstream"
	"github.com/gonum/matrix/mat64"
	"math"
	"os"
	"sort"
)

const (
	maxLagLarge int = 30 //30 samples = 30 * 100 milliseconds = 3 seconds
	maxLagSmall int = 10
)

func main() {

	maxDriverID := 35

	lagmap := map[string]int{"Speed": maxLagLarge, "FcwRange": maxLagLarge,
		"Steer":        maxLagSmall}
	floatvars1 := []string{"Driver", "Trip", "Time", "Speed",
		"Brake", "FcwValidTarget", "FcwRange", "Steer"}

	//npc := 0
	ndir := 10

	for i := 1; i <= maxDriverID; i++ {

		ivb := loadDriverDat(i, floatvars1)
		ivb = doTransforms(ivb, lagmap)
		fmt.Printf("names after doTransforms: %v\n", ivb.Names())
		ivb = laggedInteraction(ivb, "Speed", "FcwRange", maxLagLarge)
		ivb.Reset()
		ivb = dstream.Regroup(ivb, "Driver", false)

		// ---------- Fitting DOC -------------

		ivr0 := dstream.DropCols(ivb, []string{"Driver", "Brake", "Trip", "Time"})

		fmt.Printf("Variable names for DOC: %v\n", ivr0.Names())
		fmt.Printf("Variable names for original dstream: %v\n", ivb.Names())
		doc0 := dimred.NewDOC(ivr0, respvar)
		doc0.SetLogFile(fmt.Sprintf("log%03d.txt", i))
		doc0.Done()
		doc0.Fit(ndir)
		pdim := doc0.Dim() 

		xnames := make([]string, ivr0.NumVar()-1) // names of the "x" variables

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

		covfile, err := os.Create(fmt.Sprintf("/scratch/stats_flux/luers/sep_cov_%03d.txt", i))
		//(fmt.Sprintf("data/sep_cov_small_%03d.txt", i))
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
		_, err = fmt.Fprintf(wCov, "%v\n", xnames)
		if err != nil {
			panic(err)
		}
		_, err = fmt.Fprintf(wCov, "marginal eigenvalues\n%v\n", marg_evals)
		if err != nil {
			panic(err)
		}
		_, err = fmt.Fprintf(wCov, "marginal covariance\n%v\n", margcov)
		if err != nil {
			panic(err)
		}
		_, err = fmt.Fprintf(wCov, "DOC eigenvalues\n%v\n", doc0.Eig())
		if err != nil {
			panic(err)
		}
		fmt.Fprintf(wCov, "--- mean of X, Y = 1 --- \n%v\n", doc0.GetMean(1))
		fmt.Fprintf(wCov, "--- mean of X, Y = 0 --- \n%v\n", doc0.GetMean(0))
		_, err = fmt.Fprintf(wCov, "covariance, y=1\n%v\n", doc0.GetCov(1))
		if err != nil {
			panic(err)
		}
		_, err = fmt.Fprintf(wCov, "covariance, y=0\n%v\n", doc0.GetCov(0))
		if err != nil {
			panic(err)
		}
		wCov.Flush()
	}
}
