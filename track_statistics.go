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

// Compute the sample autocorrelation. If [normalize[ is true the result is
// normalized by mean and variance. The arguments [from] and [to] specify the
// range of the delay in basepairs.
func (track Track) Autocorrelation(from, to int, normalize bool) (x []int, y []float64, err error) {
  if from < 0 || to < from {
    err = fmt.Errorf("Autocorrelation(): invalid parameters")
    return
  }
  n := (to-from)/track.Binsize // number of points in the resulting  autocorrelation
  m := 0.0                     // number of data points
  // sample mean and covariance
  mean     := 0.0
  variance := 1.0
  // allocate result
  x = make([]int,     n)
  y = make([]float64, n)
  if normalize {
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
      mean     = m/(m+k)*mean     + 1/(m+k)*s
      variance = m/(m+k)*variance + 1/(m+k)*t
      m += k
    }
    variance -= mean*mean
  }
  // compute delays used for indexing (i.e. normalized by binsize)
  for j, l := 0, from; l < to; j, l = j+1, l+track.Binsize {
    x[j] = l/track.Binsize
  }
  // compute autocorrelation
  m = 0
  for _, sequence := range track.Data {
    s := make([]float64, n)
    // loop over sequence
    for i := 0; i < len(sequence); i++ {
      for j := 0; j < n && i+x[j] < len(sequence); j++ {
        s[j] += (sequence[i]-mean)*(sequence[i+x[j]]-mean)
      }
    }
    k := float64(len(sequence))
    for j := 0; j < n ; j++ {
      y[j] = m/(m+k)*y[j] + 1/(m+k)*s[j]
    }
    m += k
  }
  // normalize result and convert delays
  for j := 0; j < n ; j++ {
    x[j] *= track.Binsize
    y[j] /= variance
  }
  return
}

/* -------------------------------------------------------------------------- */

// Compute the sample cross-correlation between track1 and track2. The
// arguments [from] and [to] specify the range of the delay in basepairs.
func (track1 Track) Crosscorrelation(track2 Track, from, to int) (x []int, y []float64, err error) {
  if from < 0 || to < from {
    err = fmt.Errorf("Crosscorrelation(): invalid parameters")
    return
  }
  if track1.Binsize != track2.Binsize {
    err = fmt.Errorf("Crosscorrelation(): track binsizes do not match")
  }
  for name, sequence1 := range track1.Data {
    sequence2, ok := track2.Data[name]
    // skip sequence if it is not present in the second track
    if !ok {
      continue
    }
    if len(sequence1) != len(sequence2) {
      err = fmt.Errorf("Crosscorrelation(): track sequence lengths do not match")
    }
  }
  b := track1.Binsize
  n := (to-from)/b  // number of points in the resulting  autocorrelation
  m := 0.0          // number of data points
  // allocate result
  x = make([]int,     n)
  y = make([]float64, n)
  // compute delays used for indexing (i.e. normalized by binsize)
  for j, l := 0, from; l < to; j, l = j+1, l+b {
    x[j] = l/b
  }
  // compute autocorrelation
  m = 0
  for name, sequence1 := range track1.Data {
    sequence2, ok := track2.Data[name]
    // skip sequence if it is not present in the second track
    if !ok {
      continue
    }
    s := make([]float64, n)
    // loop over sequence
    for i := 0; i < len(sequence1); i++ {
      for j := 0; j < n && i+x[j] < len(sequence1); j++ {
        s[j] += sequence1[i]*sequence2[i+x[j]]
      }
    }
    k := float64(len(sequence1))
    for j := 0; j < n ; j++ {
      y[j] = m/(m+k)*y[j] + 1/(m+k)*s[j]
    }
    m += k
  }
  // convert delays
  for j := 0; j < n ; j++ {
    x[j] *= b
  }
  return
}
