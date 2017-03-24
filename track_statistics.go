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

type BinSummaryStatistics func(sum, min, max, n float64) float64

func BinMean(sum, min, max, n float64) float64 {
  return sum/n
}
func BinMax (sum, min, max, n float64) float64 {
  return max
}
func BinMin (sum, min, max, n float64) float64 {
  return min
}

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

/* estimate mean fragment length
 * -------------------------------------------------------------------------- */

func EstimateFragmentLength(reads GRanges, genome Genome, maxDelay, binsize int, fraglenRange [2]int) (int, []int, []float64, error) {

  if reads.Length() == 0 {
    return -1, nil, nil, fmt.Errorf("no reads given")
  }

  // compute mean read length
  readLength := uint64(0)
  for i := 0; i < reads.Length(); i++ {
    readLength += uint64(reads.Ranges[i].To - reads.Ranges[i].From)
  }
  readLength /= uint64(reads.Length())

  x, y, err := CrosscorrelateReads(reads, genome, maxDelay, binsize)

  if err != nil {
    return -1, nil, nil, err
  }
  // set initial feasible range
  from := int(readLength + readLength/2)
  to   := maxDelay
  // check if a feasible range is given
  if fraglenRange[0] != -1 {
    from = fraglenRange[0]
  }
  if fraglenRange[1] != -1 {
    to   = fraglenRange[1]
  }
  // find peak
  i_max := -1
  v_max := 0.0
  for i := 1; i < len(x)-1; i++ {
    // skip everything close to the read length (i.e. fantom peaks)
    if x[i] < from {
      continue
    }
    if x[i] >= to {
      break
    }
    // test for local maximum
    if y[i-1] < y[i] && y[i] > y[i+1] {
      // update global maximum if necessary
      if v_max < y[i] {
        i_max = i
        v_max = y[i]
      }
    }
  }
  if i_max == -1 {
    return -1, nil, nil, fmt.Errorf("no crosscorrelation peak found")
  }
  if v_max < y[len(y)-1] {
    return -1, nil, nil, fmt.Errorf("it seems that maxDelay is too small")
  }
  return x[i_max], x, y, nil
}

/* -------------------------------------------------------------------------- */

// Compute the sample autocorrelation. If [normalize] is true the result is
// normalized by mean and variance. The arguments [from] and [to] specify the
// range of the delay in basepairs.
func (track Track) Autocorrelation(from, to int, normalize bool) (x []int, y []float64, err error) {
  return track.Crosscorrelation(track, from, to, normalize)
}

/* -------------------------------------------------------------------------- */

type TrackSummaryStatistics struct {
  Name string
  Mean float64
  Max  float64
  Min  float64
}

func (statistics TrackSummaryStatistics) String() string {
  s := "Track `%s' summary statistics\n"
  s += "- Maximum: %f\n"
  s += "- Minimum: %f\n"
  s += "- Mean   : %f"
  return fmt.Sprintf(s, statistics.Name, statistics.Max, statistics.Min, statistics.Mean)
}

func (track Track) SummaryStatistics() TrackSummaryStatistics {
  statistics := TrackSummaryStatistics{}
  statistics.Name = track.Name
  statistics.Max  = math.Inf(-1)
  statistics.Min  = math.Inf(+1)
  n := 0

  // compute number of data points
  for _, sequence := range track.Data {
    for i := 0; i < len(sequence); i++ {
      if !math.IsNaN(sequence[i]) {
        n++
      }
    }
  }
  // compute statistics
  for _, sequence := range track.Data {
    for i := 0; i < len(sequence); i++ {
      // skip NaN values
      if math.IsNaN(sequence[i]) {
        continue
      }
      if sequence[i] < statistics.Min {
        statistics.Min = sequence[i]
      }
      if sequence[i] > statistics.Max {
        statistics.Max = sequence[i]
      }
      // update mean
      statistics.Mean += 1.0/float64(n)*sequence[i]
    }
  }
  return statistics
}

/* -------------------------------------------------------------------------- */

type TrackHistogram struct {
  Name string
  X    []float64
  Y    []float64
}

func (histogram TrackHistogram) String() string {
  s := "Track `%s' histogram\n"
  s += "- X: %v\n"
  s += "- Y: %v"
  return fmt.Sprintf(s, histogram.Name, histogram.X, histogram.Y)
}

func (track Track) Histogram(from, to float64, bins int) TrackHistogram {

  histogram := TrackHistogram{}

  if from >= to || bins <= 0 {
    return histogram
  }
  // allocate memory
  histogram.Name = track.Name
  histogram.X = make([]float64, bins)
  histogram.Y = make([]float64, bins)

  c := float64(bins)/(to-from)

  // compute x values
  for i := 0; i < bins; i++ {
    histogram.X[i] = from + float64(i)/c
  }
  // compute y values
  for _, sequence := range track.Data {
    for i := 0; i < len(sequence); i++ {
      // skip NaN values
      if math.IsNaN(sequence[i]) {
        continue
      }
      j := int(math.Floor((sequence[i] - from)*c))

      if j >= 0 && j < bins {
        histogram.Y[j] += 1.0
      }
    }
  }
  return histogram
}

func (track Track) CumulativeHistogram(from, to float64, bins int) TrackHistogram {

  histogram := track.Histogram(from, to, bins)
  // cumulative sum
  sum := 0.0

  for i := 0; i < len(histogram.Y); i++ {
    sum += histogram.Y[i]
    histogram.Y[i] = sum
  }
  return histogram
}
