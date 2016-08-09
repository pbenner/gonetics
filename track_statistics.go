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
import "math"

/* -------------------------------------------------------------------------- */

// Compute the sample cross-correlation between track1 and track2. If
// [normalize] is true the result is normalized by mean and variance. The
// arguments [from] and [to] specify the range of the delay in basepairs.
func (track1 Track) Crosscorrelation(track2 Track, from, to int, normalize bool) (x []int, y []float64, err error) {
  if from < 0 || to < from {
    err = fmt.Errorf("Crosscorrelation(): invalid parameters")
    return
  }
  if track1.Binsize != track2.Binsize {
    err = fmt.Errorf("Crosscorrelation(): track binsizes do not match")
    return
  }
  for name, sequence1 := range track1.Data {
    sequence2, ok := track2.Data[name]
    // skip sequence if it is not present in the second track
    if !ok {
      continue
    }
    if len(sequence1) != len(sequence2) {
      err = fmt.Errorf("Crosscorrelation(): track sequence lengths do not match")
      return
    }
  }
  b := track1.Binsize
  n := (to-from)/b  // number of points in the resulting  autocorrelation
  m := 0.0          // number of data points
  // sample mean and covariance
  mean1     := 0.0
  mean2     := 0.0
  variance1 := 1.0
  variance2 := 1.0
  // allocate result
  x = make([]int,     n)
  y = make([]float64, n)
  // compute delays used for indexing (i.e. normalized by binsize)
  for j, l := 0, from; l < to; j, l = j+1, l+b {
    x[j] = l/b
  }
  if normalize {
    // compute mean and covariance
    for name, sequence1 := range track1.Data {
      sequence2, ok := track2.Data[name]
      // skip sequence if it is not present in the second track
      if !ok {
        continue
      }
      s1 := 0.0
      s2 := 0.0
      t1 := 0.0
      t2 := 0.0
      // loop over sequence
      for i := 0; i < len(sequence1); i++ {
        s1 += sequence1[i]
        s2 += sequence2[i]
        t1 += sequence1[i]*sequence1[i]
        t2 += sequence2[i]*sequence2[i]
      }
      k := float64(len(sequence1))
      mean1     = m/(m+k)*mean1     + 1/(m+k)*s1
      mean2     = m/(m+k)*mean2     + 1/(m+k)*s2
      variance1 = m/(m+k)*variance1 + 1/(m+k)*t1
      variance2 = m/(m+k)*variance2 + 1/(m+k)*t2
      m += k
    }
    variance1 -= mean1*mean1
    variance2 -= mean2*mean2
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
        s[j] += (sequence1[i]-mean1)*(sequence2[i+x[j]]-mean2)
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
    y[j] /= math.Sqrt(variance1*variance2)
  }
  return
}

/* -------------------------------------------------------------------------- */

// Compute crosscorrelation between reads on the forward and reverse strand.
func CrosscorrelateReads(reads GRanges, genome Genome, maxDelay, binsize int) ([]int, []float64, error) {
  track1 := AllocTrack("forward", genome, binsize)
  track2 := AllocTrack("reverse", genome, binsize)

  s := reads.SetLengths(1)
  // move reads on the reverse strand one bp to the right
  for i := 0; i < s.Length(); i++ {
    if s.Strand[i] == '-' {
      s.Ranges[i].From++
      s.Ranges[i].To  ++
    }
  }
  // split reads into forward and reverse strand
  forward := s.FilterStrand('+')
  reverse := s.FilterStrand('-')
  track1.AddReads(forward, 0)
  track2.AddReads(reverse, 0)

  return track1.Crosscorrelation(track2, 0, maxDelay, true)
}

/* -------------------------------------------------------------------------- */

// Compute the sample autocorrelation. If [normalize] is true the result is
// normalized by mean and variance. The arguments [from] and [to] specify the
// range of the delay in basepairs.
func (track Track) Autocorrelation(from, to int, normalize bool) (x []int, y []float64, err error) {
  return track.Crosscorrelation(track, from, to, normalize)
}