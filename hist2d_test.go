package main

import (
       "fmt"
       "testing"
)


func TestHist2d(t *testing.T){
      data1 := []float64{-4.0, -4.0, -2.0, 0.0, -2.0, -4.0}
      data2 := []float64{0.0, 0.0, 2.0, 4.0, 0.0, 4.0}
      fmt.Printf("unsorted data1: %v\n", data1)
      fmt.Printf("unsorted data2: %v\n", data2)
      brtest := []float64{1.0, 1.0, 0.0, 0.0, 1.0, 1.0}
      fmt.Printf("corresponding indicators: %v\n", brtest)
      fmt.Println()
     b1, b2, counts := hist2d(data1, data2, brtest, 2)
     fmt.Printf("data1: %v\n",data1)
     fmt.Printf("data2: %v\n", data2)

     fmt.Printf("b1: %v\n", b1)
     fmt.Printf("b2: %v\n", b2)
     fmt.Printf("counts: %v\n", counts)
     
}