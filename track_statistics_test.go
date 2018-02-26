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

  y1 := []float64{
     0.005469647585207884,  0.026416616452943666,  0.0172013186746760,  0.01854281422250843,  0.0156587590929066,
    -0.022210691999063996,  0.006276843530171762, -0.0134952475741970, -0.00159695870026277,  0.0012526501798705,
     0.007118485724604633,  0.007250657259603093, -0.0303384641978693, -0.02966409307884197, -0.0108739022872460,
    -0.001991840077911761, -0.008190911191579012,  0.0096702852216801, -0.00501680661679200, -0.0350469017277961,
    -0.0017546965671620916 }

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
