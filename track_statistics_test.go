/* Copyright (C) 2016 Philipp Benner
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package gonetics

/* -------------------------------------------------------------------------- */

//import   "fmt"
import   "math"
import   "testing"

/* -------------------------------------------------------------------------- */

func TestTrackCrosscorrelation(t *testing.T) {

  // Reproduce cross-correlation diagram from GNU-R csaw package:
  //
  // n <- 20
  // bamFile <- system.file("exdata", "rep1.bam", package="csaw")
  //
  // x <- correlateReads(bamFile, max.dist=n)
  // plot(1:(n+1), x, xlab="delay (bp)", ylab="ccf")
  //
  // output:
  
  y1 := []float64{
     0.022738,  0.037315,  0.039222,  0.034716,  0.013813, -0.011113,
    -0.011442, -0.011026, -0.007629,  0.003967,  0.011386, -0.004300,
    -0.041619, -0.050234, -0.026522, -0.011391, -0.004422,  0.002634,
    -0.017630, -0.037979, -0.026556 }

  // printVector := func(x []float64) {
  //   for i := 0; i < len(x); i++ {
  //     if i != 0 {
  //       fmt.Printf(", ")
  //     }
  //     fmt.Printf("%f", x[i])
  //   }
  //   fmt.Println()
  // }

  genome := Genome{}
  genome.Import("track_statistics_test.genome")
  reads  := NewEmptyGRanges(0)
  reads.ImportBed6("track_statistics_test.bed")

  _, y2, _, err := CrosscorrelateReads(reads.AsReadChannel(), genome, 21, 1)

  if err != nil {
    t.Error(err); return
  }
  if len(y1) != len(y2) {
    t.Error("cross-correlation test failed")
  }
  for i := 0; i < len(y1); i++ {
    if math.Abs(y1[i]-y2[i]) > 0.006 {
      t.Error("cross-correlation test failed")
    }
  }
}
