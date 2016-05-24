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

import "bufio"
import "errors"
import "fmt"
import "math"
import "os"

import . "github.com/pbenner/pshape/Utility"

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

func NewTrack(name string, genome Genome, binsize int) Track {
  data := make(TMapType)

  for i := 0; i < genome.Length(); i++ {
    data[genome.Seqnames[i]] = make([]float64,
      divIntDown(genome.Lengths[i], binsize))
  }
  return Track{name, data, binsize}
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

/* add read counts to the track
 * -------------------------------------------------------------------------- */

// Add reads to track. All reads are extended in 3' direction to have
// a length of `d'. This is the same as the macs2 `extsize' parameter.
func (track Track) AddReads(reads GRanges, d int) {
  sum_reads_outside := 0
  for i := 0; i < reads.Length(); i++ {
    seq, ok := track.Data[reads.Seqnames[i]]
    if !ok {
      continue
    }
    from := reads.Ranges[i].From
    to   := reads.Ranges[i].To
    if to - from < d {
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
  track1 := NewTrack(name, genome, binsize)
  track2 := NewTrack("",   genome, binsize)
  for _, r := range treatment {
    track1.AddReads(r, d)
  }
  for _, r := range control {
    track2.AddReads(r, d)
  }
  for seqname, _ := range track1.Data {
    seq1     := track1.Data[seqname]
    seq2, ok := track2.Data[seqname]
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
  return track1
}

/* i/o
 * -------------------------------------------------------------------------- */

// Export the track to wiggle format. If fixedStep is false, a value is
// printed only if it is not zero. For sparse data this will significantly
// reduce the size of the file.
func (track Track) WriteWiggle(filename, description string, fixedStep bool) {
  f, err := os.Create(filename); Check(err)
  defer f.Close()

  w := bufio.NewWriter(f)
  defer w.Flush()

  fmt.Fprintf(w, "track type=wiggle_0 name=\"%s\" description=\"%s\"\n", track.Name, description)

  if fixedStep {
    for seqname, seq := range track.Data {
      fmt.Fprintf(w, "fixedStep chrom=%s start=1 step=%d\n", seqname, track.Binsize)

      for i := 0; i < len(seq); i++ {
        fmt.Fprintf(w, "%f\n", seq[i])
      }
    }
  } else {
    for seqname, seq := range track.Data {
      fmt.Fprintf(w, "variableStep chrom=%s span=%d\n", seqname, track.Binsize)

      for i := 0; i < len(seq); i++ {
        if seq[i] != 0.0 {
          fmt.Fprintf(w, "%d %f\n", i*track.Binsize+1, seq[i])
        }
      }
    }
  }
}
