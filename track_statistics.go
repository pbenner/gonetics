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

import "fmt"

/* -------------------------------------------------------------------------- */

type AcResult struct {
  X      []int
  Y      []float64
  Mean     float64
  Variance float64
}

/* -------------------------------------------------------------------------- */

// Compute the sample autocorrelation normalized by mean and variance. The
// arguments [from] and [to] specify the range of the delay in basepairs.
func (track Track) Autocorrelation(from, to int) (AcResult, error) {
  if from < 0 || to < from {
    return AcResult{}, fmt.Errorf("Autocorrelation(): invalid parameters")
  }
  n := (to-from)/track.Binsize // number of points in the resulting  autocorrelation
  m := 0.0                     // number of data points
  // allocate result
  acr  := AcResult{}
  acr.X = make([]int,     n)
  acr.Y = make([]float64, n)
  // compute mean and covariance
  for _, sequence := range track.Data {
    s := 0.0
    t := 0.0
    // loop over sequence
    for i := 0; i < len(sequence); i++ {
      s += sequence[i]
      t += sequence[i]*sequence[i]
    }
    k := float64(len(sequence))
    acr.Mean     = m/(m+k)*acr.Mean     + 1/(m+k)*s
    acr.Variance = m/(m+k)*acr.Variance + 1/(m+k)*t
    m += k
  }
  acr.Variance -= acr.Mean*acr.Mean
  // compute delays used for indexing (i.e. normalized by binsize)
  for j, l := 0, from; l < to; j, l = j+1, l+track.Binsize {
    acr.X[j] = l/track.Binsize
  }
  // compute autocorrelation
  m = 0
  for _, sequence := range track.Data {
    s := make([]float64, n)
    // loop over sequence
    for i := 0; i < len(sequence); i++ {
      for j := 0; j < n && i+acr.X[j] < len(sequence); j++ {
        s[j] += (sequence[i]-acr.Mean)*(sequence[i+acr.X[j]]-acr.Mean)
      }
    }
    k := float64(len(sequence))
    for j := 0; j < n ; j++ {
      acr.Y[j] = m/(m+k)*acr.Y[j] + 1/(m+k)*s[j]
    }
    m += k
  }
  // normalize result and convert delays
  for j := 0; j < n ; j++ {
    acr.X[j] *= track.Binsize
    acr.Y[j] /= acr.Variance
  }
  return acr, nil
}
