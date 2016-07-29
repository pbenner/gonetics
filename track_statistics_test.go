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
     0.0005482399,
     0.0211633071,  0.0129903831,  0.0138011575,  0.0106894531,
    -0.0274756566,  0.0017533539, -0.0181861776, -0.0067751107,
    -0.0038123595,  0.0020096182,  0.0022340036, -0.0357639711,
    -0.0350814350, -0.0163243527, -0.0074565095, -0.0130368542,
     0.0046453726, -0.0101288881, -0.0408354593, -0.0066084254 }

  // printVector := func(x []float64) {
  //   for i := 0; i < len(x); i++ {
  //     if i != 0 {
  //       fmt.Printf(", ")
  //     }
  //     fmt.Printf("%f", x[i])
  //   }
  //   fmt.Println()
  // }

  genome := ReadGenome("track_statistics_test.genome")
  reads  := NewEmptyGRanges(0)
  reads.ReadBed6("track_statistics_test.bed")

  _, y2, err := CrosscorrelateReads(reads, genome, 21, 1)

  if err != nil {
    t.Error(err)
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
