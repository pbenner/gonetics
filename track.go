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

import "errors"
import "fmt"
import "math"
import "os"
import "sort"

/* -------------------------------------------------------------------------- */

type TMapType map[string][]float64

// A track is a container for experimental data mapped to genomic
// locations. The data is binned in order to reduce memory usage.
// The first position in a sequence is numbered 0.
type Track struct {
  Name    string
  Data    TMapType
  Binsize int
}

/* constructor
 * -------------------------------------------------------------------------- */

func NewTrack(name string, seqnames []string, sequences [][]float64, binsize int) Track {
  if len(seqnames) != len(sequences) {
    panic("invalid arguments")
  }
  data := make(TMapType)

  for i := 0; i < len(seqnames); i++ {
    data[seqnames[i]] = sequences[i]
  }
  return Track{name, data, binsize}
}

func AllocTrack(name string, genome Genome, binsize int) Track {
  data := make(TMapType)

  for i := 0; i < genome.Length(); i++ {
    // by convention drop the last positions if they do not fully
    // cover the last bin (i.e. round down), this is required by
    // wig related tools
    data[genome.Seqnames[i]] = make([]float64,
      divIntDown(genome.Lengths[i], binsize))
  }
  return Track{name, data, binsize}
}

func EmptyTrack(name string) Track {
  data := make(TMapType)
  return Track{name, data, 0}
}

/* access methods
 * -------------------------------------------------------------------------- */

func (track Track) Clone() Track {
  name    := track.Name
  binsize := track.Binsize
  data    := make(TMapType)

  for name, sequence := range track.Data {
    t := make([]float64, len(sequence))
    copy(t, sequence)
    data[name] = t
  }
  return Track{name, data, binsize}
}

func (track Track) Index(position int) int {
  if position < 0 {
    panic("negative position")
  }
  return position/track.Binsize
}

func (track Track) Length(seqname string) (int, error) {
  seq, ok := track.Data[seqname]
  if !ok {
    return 0, errors.New("invalid seqname")
  }
  return len(seq)*track.Binsize, nil
}

func (track Track) At(seqname string, position int) (float64, error) {
  seq, ok := track.Data[seqname]
  if !ok {
    return 0, errors.New("invalid seqname")
  }
  return seq[track.Index(position)], nil
}

func (track Track) Set(seqname string, position int, value float64) error {
  seq, ok := track.Data[seqname]
  if !ok {
    return errors.New("invalid seqname")
  }
  seq[track.Index(position)] = value

  return nil
}

func (track Track) Add(seqname string, position int, value float64) error {
  seq, ok := track.Data[seqname]
  if !ok {
    return errors.New("invalid seqname")
  }
  seq[track.Index(position)] += value

  return nil
}

func (track Track) Sub(seqname string, position int, value float64) error {
  seq, ok := track.Data[seqname]
  if !ok {
    return errors.New("invalid seqname")
  }
  seq[track.Index(position)] -= value

  return nil
}

func (track Track) GetSlice(r GRangesRow) ([]float64, error) {
  seq, ok := track.Data[r.Seqnames[r.i]]
  if !ok {
    return nil, errors.New("invalid seqname")
  }
  from := r.Ranges[r.i].From/track.Binsize
  to   := r.Ranges[r.i].To  /track.Binsize
  if from >= len(seq) {
    return nil, nil
  }
  if to < 0 {
    return nil, nil
  }
  if from < 0 {
    from = 0
  }
  if to > len(seq) {
    to = len(seq)
  }
  return seq[from:to], nil
}

/* add read counts to the track
 * -------------------------------------------------------------------------- */

// Add reads to track. All reads are extended in 3' direction to have
// a length of [d]. This is the same as the macs2 `extsize' parameter.
// Reads are not extended if [d] is zero.
func (track Track) AddReads(reads GRanges, d int) {
  sum_reads_outside := 0
  for i := 0; i < reads.Length(); i++ {
    seq, ok := track.Data[reads.Seqnames[i]]
    if !ok {
      continue
    }
    from := reads.Ranges[i].From
    to   := reads.Ranges[i].To
    if d != 0 && to - from < d {
      // extend read in 3' direction
      if reads.Strand[i] == '+' {
        to = from + d - 1
      } else if reads.Strand[i] == '-' {
        from = to - d + 1
        if from < 0 { from = 0 }
      } else {
        panic("AddReads(): no strand information given!")
      }
    }
    for j := track.Index(from); j <= track.Index(to); j++ {
      jfrom := iMax(from, (j+0)*track.Binsize)
      jto   := iMin(to  , (j+1)*track.Binsize)
      if j >= len(seq) {
        sum_reads_outside++
        break
      } else {
        seq[j] += float64(jto-jfrom)/float64(track.Binsize)
      }
    }
  }
  if sum_reads_outside > 0 {
    fmt.Fprintf(os.Stderr, "AddReads(): %d read(s) are outside the genome!\n",
      sum_reads_outside)
  }
}

// Combine treatment and control from a ChIP-seq experiment into a single track.
// At each genomic location, the number of binned reads from the treatment
// experiment is divided by the number of control reads. To avoid division by
// zero, a pseudocount is added to both treatment and control. The parameter
// d determines the extension of reads.
func NormalizedTrack(name string, treatment, control []GRanges, genome Genome, d, binsize int, c1, c2 float64, logScale bool) Track {
  track1 := AllocTrack(name, genome, binsize)
  track2 := AllocTrack("",   genome, binsize)
  for _, r := range treatment {
    track1.AddReads(r, d)
  }
  for _, r := range control {
    track2.AddReads(r, d)
  }
  track1.Normalize(track2, c1, c2, logScale)

  return track1
}

func (treatment Track) Normalize(control Track, c1, c2 float64, logScale bool) {
  for seqname, _ := range treatment.Data {
    seq1     := treatment.Data[seqname]
    seq2, ok := control  .Data[seqname]
    if !ok {
      panic(fmt.Sprintf("Control has no reads on sequence `%s'!", seqname))
    }
    for i := 0; i < len(seq1); i++ {
      if logScale {
        seq1[i] = math.Log((seq1[i]+c1)/(seq2[i]+c2)*c2/c1)
      } else {
        seq1[i] = (seq1[i]+c1)/(seq2[i]+c2)*c2/c1
      }
    }
  }
}

// Smoothen track data with an adaptive window method. For each region the smallest window
// size among windowSizes is selected which contains at least minCounts counts. If the
// minimum number of counts is not reached, the larges window size is selected.
func (track Track) Smoothen(minCounts float64, windowSizes []int) {
  if len(windowSizes) == 0 {
    return
  }
  sumSlice := func(s []float64) float64 {
    sum := 0.0
    for i := 0; i < len(s); i++ {
      sum += s[i]
    }
    return sum
  }
  offset1 := divIntUp  (windowSizes[0]-1, 2)
  offset2 := divIntDown(windowSizes[0]-1, 2)
  // sort window sizes so that the smalles window size comes first
  sort.Ints(windowSizes)
  // number of window sizes
  nw := len(windowSizes)
  // loop over sequences
  for _, seq := range track.Data {
    rst := make([]float64, len(seq))
    // loop over sequence
    for i := offset1; i < len(seq)-offset2; i++ {
      counts := math.Inf(-1)
      wsize  := -1
      for k := 0; counts < minCounts && k < nw; k++ {
        from := i - divIntUp  (windowSizes[k]-1, 2)
        to   := i + divIntDown(windowSizes[k]-1, 2)
        if from < 0 {
          to   = iMin(len(seq)-1, to-from)
          from = 0
        }
        if to >= len(seq) {
          from = iMax(0, from-(to-len(seq)+1))
          to   = len(seq)-1
        }
        counts = sumSlice(seq[from:to+1])
        wsize  = to-from+1
      }
      if wsize != -1 {
        rst[i] = counts/float64(wsize)
      }
    }
    copy(seq, rst)
  }
}

/* map/reduce
 * -------------------------------------------------------------------------- */

func (track *Track) Map(f func(float64) float64) {
  for _, sequence := range track.Data {
    for i := 0; i < len(sequence); i++ {
      sequence[i] = f(sequence[i])
    }
  }
}

func (track *Track) Reduce(f func(float64, float64) float64, x0 float64) map[string]float64 {
  result := make(map[string]float64)

  for name, sequence := range track.Data {
    if len(sequence) == 0 {
      continue
    }
    tmp := f(x0, sequence[0])

    for i := 1; i < len(sequence); i++ {
      tmp = f(tmp, sequence[i])
    }
    result[name] = tmp
  }
  return result
}
