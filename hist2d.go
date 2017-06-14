package main

import (
	"fmt"
	"github.com/gonum/floats"
)

func getBins(x1sorted, x2sorted []float64, nbins int) ([]float64, []float64){
     n := len(x1sorted)
     x1width := x1sorted[n-1] - x1sorted[0]
     x2width := x2sorted[n-1] - x2sorted[0]
     x1step := x1width / float64(nbins)
     x2step := x2width / float64(nbins)

     b1 := make([]float64, nbins)
     b2 := make([]float64, nbins)
     for i := 0; i < nbins; i++{
     	 b1[i] = x1sorted[0] + float64(i) * x1step
	 b2[i] = x2sorted[0] + float64(i) * x2step
     }
     return b1, b2
}

//firstLeft returns the index i of the first element of dat
// such that dat[i] >= v
// dat is assumed to be sorted
func firstLeft(dat []float64, v float64) int{
     for i, x := range dat{
     	 if x >= v {
	    return i
	 }
     }
     return -1
}

//hist2d computes a 2-dimensional histogram with
// nbins evenly spaced bins along the dimensions x1 and x2
// br contains {1, 0} indicators
// returns the lower-left corners of each rectangle 
// and separately computed counts for the indicators contained in
// br (i.e. how many 1s and how many zeroes are in each rectangle)
func hist2d(x1, x2, br []float64, nbins int) ([]float64, []float64, [][]int){
     n := len(x1)
     if n != len(x2) {
     	panic(fmt.Errorf("unequally sized data slices"))
     }
     
     ix1 := make([]int, n)
     ix2 := make([]int, n)
     floats.Argsort(x1, ix1)
     floats.Argsort(x2, ix2)
//     fmt.Printf("sorted ix1: %v\n", ix1)
//     fmt.Printf("sorted ix2: %v\n", ix2)
     b1, b2 := getBins(x1, x2, nbins) // lower-left corner bin boundaries
     
     datIx := 0
     assign1 := make([]int, n)
     assign2 := make([]int, n)
     for bix := 1; bix < nbins; bix++{
     	 datIx += firstLeft(x1[datIx:], b1[bix])
	 for j := datIx; j < n; j++{
	     assign1[ix1[j]]++
	 }
     }
     
     datIx = 0
     for bix :=1; bix < nbins; bix++{
     	 datIx += firstLeft(x2[datIx:], b2[bix])
	 for j := datIx; j < n; j++{
	     assign2[ix2[j]]++
	 }
     }
//     fmt.Printf("assign1: %v\n", assign1)
//     fmt.Printf("assign2: %v\n", assign2)

     counts := make([][]int, nbins*nbins)
     for i := 0; i < nbins*nbins; i++{
     	 counts[i] = make([]int, 2)
     }

     for datIx = 0; datIx < n; datIx++{
     	 if br[datIx] > 0{
	    counts[assign1[datIx] + nbins * assign2[datIx]][1]++
	 } else {
	    counts[assign1[datIx] + nbins * assign2[datIx]][0]++
	 }
     }
     b1_coord := make([]float64, nbins*nbins)
     b2_coord := make([]float64, nbins*nbins)
     for i := 0; i < nbins; i++{
     	 for j := 0; j < nbins; j++{
	     b1_coord[j + i * nbins] = b1[i]
	     b2_coord[j + i * nbins] = b2[j]
	 }
     }
//     fmt.Printf("counts: %v\n", counts)
     return b1_coord, b2_coord, counts
}


