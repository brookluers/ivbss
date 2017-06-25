package main

import (
       "fmt"
       "testing"
       "os"
       "encoding/csv"
       "github.com/brookluers/dstream/dstream"
)


func TestHist2d(t *testing.T){
     f, err := os.Open("bivariate_doc.txt")
     if err != nil {
     	panic(err)
     }
     defer f.Close()
     ds := dstream.FromCSV(f).SetChunkSize(20000).SetFloatVars([]string{"x1","x2","y"}).HasHeader()
     dss := dstream.MemCopy(ds)
     
      data1 := dstream.GetCol(dss, "x1").([]float64)
      data2 := dstream.GetCol(dss, "x2").([]float64)
      y := dstream.GetCol(dss, "y").([]float64)

      fmt.Printf("100 unsorted data1 records: %v\n", data1[0:100])
      fmt.Printf("100 unsorted data2 records: %v\n", data2[0:100])

      fmt.Printf("100 y records: %v\n", y[0:100])
      fmt.Println()
     b1, b2, counts := hist2d(data1, data2, y, 50) // 50 bins

     fmt.Printf("b1: %v\n", b1)
     fmt.Printf("b2: %v\n", b2)
     fmt.Printf("counts: %v\n", counts)
     
     fw, err := os.Create("biv_doc_hist.txt")
     if err != nil{
     	panic(err)
     }
     
     fwc := csv.NewWriter(fw)
     rec := make([]string, 4)
     fwc.Write([]string{"bin1", "bin2", "num0", "num1"})
     for i := 0; i < len(b1); i++{
     	 rec[0] = fmt.Sprintf("%v", b1[i])
	 rec[1] = fmt.Sprintf("%v", b2[i])
	 rec[2] = fmt.Sprintf("%v", counts[i][0])
	 rec[3] = fmt.Sprintf("%v", counts[i][1])
     	 fwc.Write(rec)
     }
     fwc.Flush()
     
}