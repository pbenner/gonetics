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

/* add read counts to the track
 * -------------------------------------------------------------------------- */

// Add reads to track. All reads are extended in 3' direction to have
// a length of [d]. This is the same as the macs2 `extsize' parameter.
// Reads are not extended if [d] is zero. If [addOverlap] is true, the
// percentage of overlap between reads and bins is added.
func (track GenericMutableTrack) AddReads(reads GRanges, d int, addOverlap bool) error {
  for i := 0; i < reads.Length(); i++ {
    seq, err := track.GetSequence(reads.Seqnames[i]); if err != nil {
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
        return fmt.Errorf("strand information is missing for read `%d'", i)
      }
    }
    binSize := track.GetBinSize()
    for j := from/binSize; j <= to/binSize; j++ {
      jfrom := iMax(from, (j+0)*binSize)
      jto   := iMin(to  , (j+1)*binSize)
      if j >= seq.NBins() {
        break
      } else {
        if addOverlap {
          seq.SetBin(j, seq.AtBin(j) + float64(jto-jfrom)/float64(binSize))
        } else {
          seq.SetBin(j, seq.AtBin(j) + 1.0)
        }
      }
    }
  }
  return nil
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
