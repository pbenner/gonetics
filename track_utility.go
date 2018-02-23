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

type cumDist struct {
  x []float64
  y []int
  n   int
}

func (obj cumDist) Len() int {
  return len(obj.x)
}

func (obj cumDist) Less(i, j int) bool {
  return obj.x[i] < obj.x[j]
}

func (obj cumDist) Swap(i, j int) {
  obj.x[i], obj.x[j] = obj.x[j], obj.x[i]
  obj.y[i], obj.y[j] = obj.y[j], obj.y[i]
}

func newCumDist(m map[float64]int) cumDist {
  x := make([]float64, len(m))
  y := make([]int,     len(m))
  i := 0
  for k, v := range m {
    x[i] = k
    y[i] = v
    i++
  }
  c := cumDist{x, y, 0}
  sort.Sort(c)

  // compute cumulative distribution
  n := 0
  for i := 0; i < len(x); i++ {
    n += y[i]; y[i] = n
  }
  c.n = n

  return c
}

/* add read counts to the track
 * -------------------------------------------------------------------------- */

// Add a single read to the track. Single end reads are extended in 3' direction
// to have a length of [d]. This is the same as the macs2 `extsize' parameter.
// Reads are not extended if [d] is zero. If [addOverlap] is true, the
// percentage of overlap between reads and bins is added. The function
// returns an error if the read's position is out of range
func (track GenericMutableTrack) AddRead(read Read, d int, addOverlap bool) error {
  seq, err := track.GetSequence(read.Seqname); if err != nil {
    return err
  }
  from := read.Range.From
  to   := read.Range.To
  if !read.PairedEnd && d > 0 {
    // extend read in 3' direction
    if read.Strand == '+' {
      to = from + d
    } else if read.Strand == '-' {
      from = to - d
      if from < 0 { from = 0 }
    } else {
      return fmt.Errorf("strand information is missing for read `%v'", read)
    }
  }
  binSize := track.GetBinSize()
  if from/binSize >= seq.NBins() {
    return fmt.Errorf("read %+v is out of range", read)
  }
  for j := from/binSize; j <= (to-1)/binSize; j++ {
    if j >= seq.NBins() {
      break
    } else {
      if addOverlap {
        jfrom := iMax(from, (j+0)*binSize)
        jto   := iMin(to  , (j+1)*binSize)
        seq.SetBin(j, seq.AtBin(j) + float64(jto-jfrom)/float64(binSize))
      } else {
        seq.SetBin(j, seq.AtBin(j) + 1.0)
      }
    }
  }
  return nil
}

// Add reads to track. All single end reads are extended in 3' direction
// to have a length of [d]. This is the same as the macs2 `extsize' parameter.
// Reads are not extended if [d] is zero. If [addOverlap] is true, the
// percentage of overlap between reads and bins is added. The function
// returns the number of reads added to the track
func (track GenericMutableTrack) AddReads(reads ReadChannel, d int, addOverlap bool) int {
  n := 0
  for read := range reads {
    if err := track.AddRead(read, d, addOverlap); err == nil {
      n++
    }
  }
  return n
}

// Combine treatment and control from a ChIP-seq experiment into a single track.
// At each genomic location, the number of binned reads from the treatment
// experiment is divided by the number of control reads. To avoid division by
// zero, a pseudocount is added to both treatment and control. The parameter
// d determines the extension of reads.
func (track GenericMutableTrack) Normalize(treatment, control Track, c1, c2 float64, logScale bool) error {
  if c1 <= 0.0 || c2 <= 0.0 {
    return fmt.Errorf("pseudocounts must be strictly positive")
  }
  for _, name := range track.GetSeqNames() {
    seq, err := track.GetMutableSequence(name); if err != nil {
      return err
    }
    seq1, err := treatment.GetSequence(name); if err != nil {
      return err
    }
    seq2, err := control  .GetSequence(name); if err != nil {
      continue
    }
    for i := 0; i < seq1.NBins(); i++ {
      if logScale {
        seq.SetBin(i, math.Log((seq1.AtBin(i)+c1)/(seq2.AtBin(i)+c2)*c2/c1))
      } else {
        seq.SetBin(i, (seq1.AtBin(i)+c1)/(seq2.AtBin(i)+c2)*c2/c1)
      }
    }
  }
  return nil
}

// Quantile normalize track to reference
func (track GenericMutableTrack) QuantileNormalize(trackRef Track) error {
  mapRef := make(map[float64]int)
  mapIn  := make(map[float64]int)
  mapTr  := make(map[float64]float64)

  if err := (GenericMutableTrack{}).Map(trackRef, func(seqname string, position int, value float64) float64 {
    if !math.IsNaN(value) {
      mapRef[value] += 1
    }
    return 0.0
  }); err != nil {
    return err
  }
  if err := (GenericMutableTrack{}).Map(track, func(seqname string, position int, value float64) float64 {
    if !math.IsNaN(value) {
      mapIn[value] += 1
    }
    return 0.0
  }); err != nil {
    return err
  }
  distRef := newCumDist(mapRef)
  distIn  := newCumDist(mapIn)

  if len(distRef.x) == 0 {
    return nil
  }
  // set first value to keep data on the same range
  mapTr[distIn.x[0]] = distRef.x[0]

  for i, j := 1, 1; i < len(distRef.x); i++ {
    pRef := float64(distRef.y[i])/float64(distRef.n)
    for ; j < len(distIn.x); j++ {
      pIn := float64(distIn.y[j])/float64(distIn.n)
      if pIn > pRef {
        break
      }
      // map input x_j to reference x_i
      mapTr[distIn.x[j]] = distRef.x[i]
    }
  }

  if err := track.Map(track, func(seqname string, position int, value float64) float64 {
    if math.IsNaN(value) {
      return value
    } else {
      return mapTr[value]
    }
  }); err != nil {
    return err
  }
  return nil
}

// Smoothen track data with an adaptive window method. For each region the smallest window
// size among windowSizes is selected which contains at least minCounts counts. If the
// minimum number of counts is not reached, the larges window size is selected.
func (track GenericMutableTrack) Smoothen(minCounts float64, windowSizes []int) error {
  if len(windowSizes) == 0 {
    return nil
  }
  offset1 := divIntUp  (windowSizes[0]-1, 2)
  offset2 := divIntDown(windowSizes[0]-1, 2)
  // sort window sizes so that the smalles window size comes first
  sort.Ints(windowSizes)
  // number of window sizes
  nw := len(windowSizes)
  // loop over sequences
  for _, name := range track.GetSeqNames() {
    seq, err := track.GetMutableSequence(name); if err != nil {
      return err
    }
    rst := make([]float64, seq.NBins())
    // loop over sequence
    for i := offset1; i < seq.NBins()-offset2; i++ {
      counts := math.Inf(-1)
      wsize  := -1
      for k := 0; counts < minCounts && k < nw; k++ {
        from := i - divIntUp  (windowSizes[k]-1, 2)
        to   := i + divIntDown(windowSizes[k]-1, 2)
        if from < 0 {
          to   = iMin(seq.NBins()-1, to-from)
          from = 0
        }
        if to >= seq.NBins() {
          from = iMax(0, from-(to-seq.NBins()+1))
          to   = seq.NBins()-1
        }
        counts = 0
        for i := from; i < to+1; i++ {
          counts += seq.AtBin(i)
        }
        wsize  = to-from+1
      }
      if wsize != -1 {
        rst[i] = counts/float64(wsize)
      }
    }
    for i := 0; i < seq.NBins(); i++ {
      seq.SetBin(i, rst[i])
    }
  }
  return nil
}

/* map/reduce
 * -------------------------------------------------------------------------- */

// Apply a function f to the sequences of track2. The function is given as
// arguments the name of the sequence, the position, and the value at that
// position. The return value is stored in track1 if it is not nil.
func (track1 GenericMutableTrack) Map(track2 Track, f func(string, int, float64) float64) error {
  if track1.MutableTrack == nil {
    binSize := track2.GetBinSize()
    for _, name := range track2.GetSeqNames() {
      seq2, err := track2.GetSequence(name); if err != nil {
        return err
      }
      for i := 0; i < seq2.NBins(); i++ {
        f(name, i*binSize, seq2.AtBin(i))
      }
    }
  } else {
    if track1.GetBinSize() != track2.GetBinSize() {
      return fmt.Errorf("binSizes do not match")
    }
    binSize := track1.GetBinSize()
    for _, name := range track1.GetSeqNames() {
      seq1, err := track1.GetMutableSequence(name); if err != nil {
        return err
      }
      seq2, err := track2.GetSequence(name); if err != nil {
        return err
      }
      if seq1.NBins() != seq2.NBins() {
        return fmt.Errorf("sequence lengths do not match for `%d'", name)
      }
      for i := 0; i < seq2.NBins(); i++ {
        seq1.SetBin(i, f(name, i*binSize, seq2.AtBin(i)))
      }
    }
  }
  return nil
}

// Apply a window function f to the sequences of track2. The function is given
// as arguments the name of the sequence, the position, and the value at that
// position. The return value is stored in track1 if it is not nil.
func (track1 GenericMutableTrack) WindowMap(track2 Track, windowSize int, f func(string, int, ...float64) float64) error {
  if windowSize <= 0 {
    return fmt.Errorf("invalid window size")
  }
  v := make([]float64, windowSize)
  if track1.MutableTrack == nil {
    binSize := track2.GetBinSize()
    for _, name := range track2.GetSeqNames() {
      seq2, err := track2.GetSequence(name); if err != nil {
        return err
      }
      for i := 0; i < seq2.NBins(); i++ {
        // fill v slice
        for j := 0; j < windowSize; j++ {
          if k := i - windowSize/2 + j; k < 0 || k >= seq2.NBins() {
            v[j] = math.NaN()
          } else {
            v[j] = seq2.AtBin(k)
          }
        }
        // call function
        f(name, i*binSize, v...)
      }
    }
  } else {
    if track1.GetBinSize() != track2.GetBinSize() {
      return fmt.Errorf("binSizes do not match")
    }
    binSize := track1.GetBinSize()
    for _, name := range track1.GetSeqNames() {
      seq1, err := track1.GetMutableSequence(name); if err != nil {
        return err
      }
      seq2, err := track2.GetSequence(name); if err != nil {
        return err
      }
      if seq1.NBins() != seq2.NBins() {
        return fmt.Errorf("sequence lengths do not match for `%d'", name)
      }
      for i := 0; i < seq2.NBins(); i++ {
        // fill v slice
        for j := 0; j < windowSize; j++ {
          if k := i - windowSize/2 + j; k < 0 || k >= seq2.NBins() {
            v[j] = math.NaN()
          } else {
            v[j] = seq2.AtBin(k)
          }
        }
        // call function
        seq1.SetBin(i, f(name, i*binSize, v...))
      }
    }
  }
  return nil
}

func (track GenericMutableTrack) MapList(tracks []Track, f func(string, int, ...float64) float64) error {
  if len(tracks) == 0 {
    return nil
  }
  // number of tracks
  n := len(tracks)
  v := make([]float64, n)
  if track.MutableTrack == nil {
    // check bin sizes
    for i := 1; i < n; i++ {
      if tracks[0].GetBinSize() != tracks[i].GetBinSize() {
        return fmt.Errorf("binSizes do not match")
      }
    }
    binSize := tracks[0].GetBinSize()
    for _, name := range tracks[0].GetSeqNames() {
      sequences := []TrackSequence{}
      nbins     := -1
      // collect source sequences
      for k, t := range tracks {
        if seq, err := t.GetSequence(name); err == nil {
          if nbins == -1 {
            nbins = seq.NBins()
          }
          if seq.NBins() != nbins {
            return fmt.Errorf("sequence `%s' in track `%d' has invalid length (`%d' instead of `%d')", name, k, seq.NBins(), nbins)
          }
          sequences = append(sequences, seq)
        }
      }
      // reduce length of v if some tracks are missing a sequence
      v := v[0:len(sequences)]
      // loop over sequence
      for i := 0; i < nbins; i++ {
        // copy values to local vector
        for j := 0; j < len(sequences); j++ {
          v[j] = sequences[j].AtBin(i)
        }
        // apply function
        f(name, i*binSize, v...)
      }
    }
  } else {
    // check bin sizes
    for _, t := range tracks {
      if track.GetBinSize() != t.GetBinSize() {
        return fmt.Errorf("binSizes do not match")
      }
    }
    binSize := track.GetBinSize()
    for _, name := range track.GetSeqNames() {
      dst, err := track.GetMutableSequence(name); if err != nil {
        return err
      }
      sequences := []TrackSequence{}
      // collect source sequences
      for k, t := range tracks {
        if seq, err := t.GetSequence(name); err == nil {
          if seq.NBins() != dst.NBins() {
            return fmt.Errorf("sequence `%s' in track `%d' has invalid length (`%d' instead of `%d')", name, k, seq.NBins(), dst.NBins())
          }
          sequences = append(sequences, seq)
        }
      }
      // reduce length of v if some tracks are missing a sequence
      v := v[0:len(sequences)]
      // loop over sequence
      for i := 0; i < dst.NBins(); i++ {
        // copy values to local vector
        for j := 0; j < len(sequences); j++ {
          v[j] = sequences[j].AtBin(i)
        }
        // apply function
        dst.SetBin(i, f(name, i*binSize, v...))
      }
    }
  }
  return nil
}

func (track GenericMutableTrack) WindowMapList(tracks []Track, windowSize int, f func(string, int, ...[]float64) float64) error {
  if len(tracks) == 0 {
    return nil
  }
  if windowSize <= 0 {
    return fmt.Errorf("invalid window size")
  }
  // number of tracks
  n := len(tracks)
  v := make([][]float64, n)
  for i := 0; i < n; i++ {
    v[i] = make([]float64, windowSize)
  }
  if track.MutableTrack == nil {
    // check bin sizes
    for i := 1; i < n; i++ {
      if tracks[0].GetBinSize() != tracks[i].GetBinSize() {
        return fmt.Errorf("binSizes do not match")
      }
    }
    binSize := tracks[0].GetBinSize()
    for _, name := range tracks[0].GetSeqNames() {
      sequences := []TrackSequence{}
      nbins     := -1
      // collect source sequences
      for k, t := range tracks {
        if seq, err := t.GetSequence(name); err == nil {
          if nbins == -1 {
            nbins = seq.NBins()
          }
          if seq.NBins() != nbins {
            return fmt.Errorf("sequence `%s' in track `%d' has invalid length (`%d' instead of `%d')", name, k, seq.NBins(), nbins)
          }
          sequences = append(sequences, seq)
        }
      }
      // reduce length of v if some tracks are missing a sequence
      v := v[0:len(sequences)]
      // loop over sequence
      for i := 0; i < nbins; i++ {
        // copy values to local vector
        for j := 0; j < len(sequences); j++ {
          for k := 0; k < windowSize; k++ {
            if t := i - windowSize/2 + k; t < 0 || t >= sequences[j].NBins() {
              v[j][k] = math.NaN()
            } else {
              v[j][k] = sequences[j].AtBin(t)
            }
          }
        }
        // apply function
        f(name, i*binSize, v...)
      }
    }
  } else {
    // check bin sizes
    for _, t := range tracks {
      if track.GetBinSize() != t.GetBinSize() {
        return fmt.Errorf("binSizes do not match")
      }
    }
    binSize := track.GetBinSize()
    for _, name := range track.GetSeqNames() {
      dst, err := track.GetMutableSequence(name); if err != nil {
        return err
      }
      sequences := []TrackSequence{}
      // collect source sequences
      for k, t := range tracks {
        if seq, err := t.GetSequence(name); err == nil {
          if seq.NBins() != dst.NBins() {
            return fmt.Errorf("sequence `%s' in track `%d' has invalid length (`%d' instead of `%d')", name, k, seq.NBins(), dst.NBins())
          }
          sequences = append(sequences, seq)
        }
      }
      // reduce length of v if some tracks are missing a sequence
      v := v[0:len(sequences)]
      // loop over sequence
      for i := 0; i < dst.NBins(); i++ {
        // copy values to local vector
        for j := 0; j < len(sequences); j++ {
          for k := 0; k < windowSize; k++ {
            if t := i - windowSize/2 + k; t < 0 || t >= sequences[j].NBins() {
              v[j][k] = math.NaN()
            } else {
              v[j][k] = sequences[j].AtBin(t)
            }
          }
        }
        // apply function
        dst.SetBin(i, f(name, i*binSize, v...))
      }
    }
  }
  return nil
}

func (track GenericTrack) Reduce(f func(string, int, float64, float64) float64, x0 float64) map[string]float64 {
  result  := make(map[string]float64)
  binSize := track.GetBinSize()

  for _, name := range track.GetSeqNames() {
    sequence, err := track.GetSequence(name); if err != nil {
      continue
    }
    if sequence.NBins() == 0 {
      continue
    }
    tmp := f(name, 0, x0, sequence.AtBin(0))

    for i := 1; i < sequence.NBins(); i++ {
      tmp = f(name, i*binSize, tmp, sequence.AtBin(i))
    }
    result[name] = tmp
  }
  return result
}
