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
// Reads are not extended if [d] is zero.
func TrackAddReads(track Track, reads GRanges, d int) error {
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
    binsize := track.GetBinsize()
    for j := from/binsize; j <= to/binsize; j++ {
      jfrom := iMax(from, (j+0)*binsize)
      jto   := iMin(to  , (j+1)*binsize)
      if j >= len(seq) {
        break
      } else {
        seq[j] += float64(jto-jfrom)/float64(binsize)
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
func TrackAddReadsAndNormalize(track1, track2 Track, treatment, control []GRanges, d int, c1, c2 float64, logScale bool) error {
  for _, r := range treatment {
    TrackAddReads(track1, r, d)
  }
  for _, r := range control {
    TrackAddReads(track2, r, d)
  }
  return TrackNormalize(track1, track2, c1, c2, logScale)
}

func TrackNormalize(treatment, control Track, c1, c2 float64, logScale bool) error {
  if c1 <= 0.0 || c2 <= 0.0 {
    return fmt.Errorf("pseudocounts must be strictly positive")
  }
  for _, name := range treatment.GetSeqNames() {
    seq1, err := treatment.GetSequence(name); if err != nil {
      return err
    }
    seq2, err := control  .GetSequence(name); if err != nil {
      continue
    }
    for i := 0; i < len(seq1); i++ {
      if logScale {
        seq1[i] = math.Log((seq1[i]+c1)/(seq2[i]+c2)*c2/c1)
      } else {
        seq1[i] = (seq1[i]+c1)/(seq2[i]+c2)*c2/c1
      }
    }
  }
  return nil
}

// Smoothen track data with an adaptive window method. For each region the smallest window
// size among windowSizes is selected which contains at least minCounts counts. If the
// minimum number of counts is not reached, the larges window size is selected.
func TrackSmoothen(track Track, minCounts float64, windowSizes []int) error {
  if len(windowSizes) == 0 {
    return nil
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
  for _, name := range track.GetSeqNames() {
    seq, err := track.GetSequence(name); if err != nil {
      return err
    }
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
  return nil
}

/* map/reduce
 * -------------------------------------------------------------------------- */

func TrackMap(track1, track2 Track, f func(string, int, float64) float64) error {
  if track1.GetBinsize() != track2.GetBinsize() {
    return fmt.Errorf("binsizes do not match")
  }
  binsize := track1.GetBinsize()
  for _, name := range track1.GetSeqNames() {
    seq1, err := track1.GetSequence(name); if err != nil {
      return err
    }
    seq2, err := track2.GetSequence(name); if err != nil {
      return err
    }
    if len(seq1) != len(seq2) {
      return fmt.Errorf("sequence lengths do not match for `%d'", name)
    }
    for i := 0; i < len(seq2); i++ {
      seq1[i] = f(name, i*binsize, seq2[i])
    }
  }
  return nil
}

func TrackMapList(track Track, tracks []Track, f func(string, int, ...float64) float64) error {
  // number of tracks
  n := len(tracks)
  v := make([]float64, n)
  // check bin sizes
  for _, t := range tracks {
    if track.GetBinsize() != t.GetBinsize() {
      return fmt.Errorf("binsizes do not match")
    }
  }
  binsize := track.GetBinsize()
  for _, name := range track.GetSeqNames() {
    dst, err := track.GetSequence(name); if err != nil {
      return err
    }
    sequences := [][]float64{}
    // collect source sequences
    for k, t := range tracks {
      if seq, err := t.GetSequence(name); err == nil {
        if len(seq) != len(dst) {
          return fmt.Errorf("sequence `%s' in track `%d' has invalid length (`%d' instead of `%d')", name, k, len(seq), len(dst))
        }
        sequences = append(sequences, seq)
      }
    }
    // reduce length of v if some tracks are missing a sequence
    v := v[0:len(sequences)]
    // loop over sequence
    for i := 0; i < len(dst); i++ {
      // copy values to local vector
      for j := 0; j < len(sequences); j++ {
        v[j] = sequences[j][i]
      }
      // apply function
      dst[i] = f(name, i*binsize, v...)
    }
  }
  return nil
}

func TrackReduce(track Track, f func(string, int, float64, float64) float64, x0 float64) map[string]float64 {
  result  := make(map[string]float64)
  binsize := track.GetBinsize()

  for _, name := range track.GetSeqNames() {
    sequence, err := track.GetSequence(name); if err != nil {
      continue
    }
    if len(sequence) == 0 {
      continue
    }
    tmp := f(name, 0, x0, sequence[0])

    for i := 1; i < len(sequence); i++ {
      tmp = f(name, i*binsize, tmp, sequence[i])
    }
    result[name] = tmp
  }
  return result
}
