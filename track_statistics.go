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
import "sort"

/* -------------------------------------------------------------------------- */

type BinSummaryStatistics func(sum, sumSquares, min, max, n float64) float64

func BinMean(sum, sumSquares, min, max, n float64) float64 {
  return sum/n
}
func BinMax (sum, sumSquares, min, max, n float64) float64 {
  return max
}
func BinMin (sum, sumSquares, min, max, n float64) float64 {
  return min
}
func BinDiscreteMean(sum, sumSquares, min, max, n float64) float64 {
  return math.Floor(sum/n + 0.5)
}
func BinDiscreteMax (sum, sumSquares, min, max, n float64) float64 {
  return math.Floor(max)
}
func BinDiscreteMin (sum, sumSquares, min, max, n float64) float64 {
  return math.Floor(min)
}
func BinVariance(sum, sumSquares, min, max, n float64) float64 {
  return sumSquares/n - sum/n*sum/n
}

func BinSummaryStatisticsFromString(str string) BinSummaryStatistics {
  switch str {
  case "mean":
    return BinMean
  case "max":
    return BinMax
  case "min":
    return BinMin
  case "discrete mean":
    return BinDiscreteMean
  case "discrete max":
    return BinDiscreteMax
  case "discrete min":
    return BinDiscreteMin
  case "variance":
    return BinVariance
  }
  return nil
}

/* -------------------------------------------------------------------------- */

// Compute the sample cross-correlation between track1 and track2. If
// [normalize] is true the result is normalized by mean and variance. The
// arguments [from] and [to] specify the range of the delay in basepairs.
func TrackCrosscorrelation(track1, track2 Track, from, to int, normalize bool) (x []int, y []float64, err error) {
  var sequence1 TrackSequence
  var sequence2 TrackSequence
  if from < 0 || to < from {
    err = fmt.Errorf("Crosscorrelation(): invalid parameters")
    return
  }
  if track1.GetBinSize() != track2.GetBinSize() {
    err = fmt.Errorf("Crosscorrelation(): track binSizes do not match")
    return
  }
  for _, name := range track1.GetSeqNames() {
    sequence1, err = track1.GetSequence(name); if err != nil {
      return
    }
    // skip sequence if it is not present in the second track
    sequence2, err = track2.GetSequence(name); if err != nil {
      continue
    }
    if sequence1.NBins() != sequence2.NBins() {
      err = fmt.Errorf("Crosscorrelation(): track sequence lengths do not match")
      return
    }
  }
  b := track1.GetBinSize()
  n := divIntUp(to-from, b) // number of points in the resulting autocorrelation
  m := 0.0                  // number of data points
  // sample mean and covariance
  mean1     := 0.0
  mean2     := 0.0
  variance1 := 1.0
  variance2 := 1.0
  // allocate result
  x = make([]int,     n)
  y = make([]float64, n)
  // compute delays used for indexing (i.e. normalized by binSize)
  for j, l := 0, from; l < to; j, l = j+1, l+b {
    x[j] = l/b
  }
  if normalize {
    // compute mean and covariance
    for _, name := range track1.GetSeqNames() {
      sequence1, err = track1.GetSequence(name); if err != nil {
        return
      }
      // skip sequence if it is not present in the second track
      sequence2, err = track2.GetSequence(name); if err != nil {
        continue
      }
      s1 := 0.0
      s2 := 0.0
      t1 := 0.0
      t2 := 0.0
      // loop over sequence
      for i := 0; i < sequence1.NBins(); i++ {
        s1 += sequence1.AtBin(i)
        s2 += sequence2.AtBin(i)
        t1 += sequence1.AtBin(i)*sequence1.AtBin(i)
        t2 += sequence2.AtBin(i)*sequence2.AtBin(i)
      }
      k := float64(sequence1.NBins())
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
  for _, name := range track1.GetSeqNames() {
    sequence1, err = track1.GetSequence(name); if err != nil {
      return
    }
    // skip sequence if it is not present in the second track
    sequence2, err = track2.GetSequence(name); if err != nil {
      continue
    }
    s := make([]float64, n)
    // loop over sequence
    for i := 0; i < sequence1.NBins(); i++ {
      for j := 0; j < n && i+x[j] < sequence1.NBins(); j++ {
        s[j] += (sequence1.AtBin(i)-mean1)*(sequence2.AtBin(i+x[j])-mean2)
      }
    }
    k := float64(sequence1.NBins())
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
func CrosscorrelateReads(reads ReadChannel, genome Genome, maxDelay, binSize int) ([]int, []float64, int, uint64, error) {
  track1 := AllocSimpleTrack("forward", genome, binSize)
  track2 := AllocSimpleTrack("reverse", genome, binSize)

  // number of reads
  n := uint64(0)
  // mean read length
  readLength := uint64(0)

  for read := range reads {
    if read.Strand == '+' {
      // set length to one
      read.Range.To = read.Range.From+1
      // add read
      if err := (GenericMutableTrack{track1}).AddRead(read, 0); err == nil {
        readLength += uint64(read.Range.To - read.Range.From); n++
      }
    } else
    if read.Strand == '-' {
      // set length to one
      read.Range.From = read.Range.To-1
      // move reads on the reverse strand one bp to the right
      read.Range.From++
      read.Range.To  ++
      // add read
      if err := (GenericMutableTrack{track2}).AddRead(read, 0); err == nil {
        readLength += uint64(read.Range.To - read.Range.From); n++
      }
    }
  }
  if n == 0 {
    return nil, nil, 0, n, fmt.Errorf("computing cross-correlation failed: no reads available")
  }
  readLength /= uint64(n)

  x, y, err := TrackCrosscorrelation(track1, track2, 0, maxDelay, true)

  return x, y, int(readLength), n, err
}

/* estimate mean fragment length
 * -------------------------------------------------------------------------- */

var ErrFraglenEstimate = fmt.Errorf("estimating fragment length failed")

func EstimateFragmentLength(reads ReadChannel, genome Genome, maxDelay, binSize int, fraglenRange [2]int) (int, []int, []float64, uint64, error) {

  x, y, readLength, n, err := CrosscorrelateReads(reads, genome, maxDelay, binSize)

  if err != nil {
    return -1, nil, nil, n, err
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
    return -1, x, y, n, fmt.Errorf("%w: %s", ErrFraglenEstimate, "no crosscorrelation peak found")
  }
  if v_max < y[len(y)-1] {
    return -1, x, y, n, fmt.Errorf("%w: %s", ErrFraglenEstimate, "it seems that maxDelay is too small")
  }
  return x[i_max], x, y, n, nil
}

/* -------------------------------------------------------------------------- */

// Compute the sample autocorrelation. If [normalize] is true the result is
// normalized by mean and variance. The arguments [from] and [to] specify the
// range of the delay in basepairs.
func TrackAutocorrelation(track Track, from, to int, normalize bool) (x []int, y []float64, err error) {
  return TrackCrosscorrelation(track, track, from, to, normalize)
}

/* -------------------------------------------------------------------------- */

type TrackSummaryStatistics struct {
  Name string
  Mean float64
  Max  float64
  Min  float64
  N    int
}

func (statistics TrackSummaryStatistics) String() string {
  var s string
  if statistics.Name == "" {
    s = fmt.Sprintf("Track summary statistics\n")
  } else {
    s = fmt.Sprintf("Track `%s' summary statistics\n", statistics.Name)
  }
  s += "- N      : %d\n"
  s += "- Maximum: %f\n"
  s += "- Minimum: %f\n"
  s += "- Mean   : %f"
  return fmt.Sprintf(s, statistics.N, statistics.Max, statistics.Min, statistics.Mean)
}

func (track GenericTrack) SummaryStatistics() TrackSummaryStatistics {
  statistics := TrackSummaryStatistics{}
  statistics.Name = track.GetName()
  statistics.N    = 0
  statistics.Max  = math.Inf(-1)
  statistics.Min  = math.Inf(+1)

  // compute number of data points
  for _, name := range track.GetSeqNames() {
    sequence, err := track.GetSequence(name); if err != nil {
      continue
    }
    for i := 0; i < sequence.NBins(); i++ {
      if !math.IsNaN(sequence.AtBin(i)) {
        statistics.N++
      }
    }
  }
  // compute statistics
  for _, name := range track.GetSeqNames() {
    sequence, err := track.GetSequence(name); if err != nil {
      continue
    }
    for i := 0; i < sequence.NBins(); i++ {
      // skip NaN values
      if math.IsNaN(sequence.AtBin(i)) {
        continue
      }
      if sequence.AtBin(i) < statistics.Min {
        statistics.Min = sequence.AtBin(i)
      }
      if sequence.AtBin(i) > statistics.Max {
        statistics.Max = sequence.AtBin(i)
      }
      // update mean
      statistics.Mean += 1.0/float64(statistics.N)*sequence.AtBin(i)
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

func (track GenericTrack) histogramNoBinning() TrackHistogram {
  histogram := TrackHistogram{}
  m := make(map[float64]int)
  for _, name := range track.GetSeqNames() {
    sequence, err := track.GetSequence(name); if err != nil {
      continue
    }
    for i := 0; i < sequence.NBins(); i++ {
      m[sequence.AtBin(i)]++
    }
  }
  histogram.X = make([]float64, len(m))
  histogram.Y = make([]float64, len(m))
  { i := 0
    for x, _ := range m {
      histogram.X[i] = x; i++
    }
  }
  sort.Float64s(histogram.X)
  for i := 0; i < len(histogram.X); i++ {
    histogram.Y[i] = float64(m[histogram.X[i]])
  }
  return histogram
}

func (track GenericTrack) Histogram(from, to float64, bins int) TrackHistogram {
  if bins == 0 {
    return track.histogramNoBinning()
  }
  histogram := TrackHistogram{}

  if from >= to {
    return histogram
  }
  // allocate memory
  histogram.Name = track.GetName()
  histogram.X = make([]float64, bins)
  histogram.Y = make([]float64, bins)

  c := float64(bins)/(to-from)

  // compute x values
  for i := 0; i < bins; i++ {
    histogram.X[i] = from + float64(i)/c
  }
  // compute y values
  for _, name := range track.GetSeqNames() {
    sequence, err := track.GetSequence(name); if err != nil {
      continue
    }
    for i := 0; i < sequence.NBins(); i++ {
      // skip NaN values
      if math.IsNaN(sequence.AtBin(i)) {
        continue
      }
      j := int(math.Floor((sequence.AtBin(i) - from)*c))

      if j >= 0 && j < bins {
        histogram.Y[j] += 1.0
      }
    }
  }
  return histogram
}

func (track GenericTrack) CumulativeHistogram(from, to float64, bins int) TrackHistogram {

  histogram := GenericTrack{track}.Histogram(from, to, bins)
  // cumulative sum
  sum := 0.0

  for i := 0; i < len(histogram.Y); i++ {
    sum += histogram.Y[i]
    histogram.Y[i] = sum
  }
  return histogram
}
