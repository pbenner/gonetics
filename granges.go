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
import "sort"

/* -------------------------------------------------------------------------- */

type GRanges struct {
  Seqnames   []string
  Ranges     []Range
  Strand     []byte
  Meta
}

/* constructors
 * -------------------------------------------------------------------------- */

func NewGRanges(seqnames []string, from, to []int, strand []byte) GRanges {
  n := len(seqnames)
  if len(  from) != n || len(    to) != n ||
    (len(strand) != 0 && len(strand) != n) {
    panic("NewGRanges(): invalid arguments!")
  }
  if len(strand) == 0 {
    strand = make([]byte, n)
    for i := 0; i < n; i++ {
      strand[i] = '*'
    }
  }
  ranges := make([]Range, n)
  for i := 0; i < n; i++ {
    // create range
    ranges[i] = NewRange(from[i], to[i])
    // check if strand is valid
    if strand[i] != '+' && strand[i] != '-' && strand[i] != '*' {
      panic("NewGRanges(): Invalid strand!")
    }
  }
  return GRanges{seqnames, ranges, strand, Meta{}}
}

func NewEmptyGRanges(n int) GRanges {
  seqnames := make([]string, n)
  ranges   := make([]Range, n)
  strand   := make([]byte, n)
  for i := 0; i < n; i++ {
    strand[i] = '*'
  }
  return GRanges{seqnames, ranges, strand, Meta{}}
}

func (r *GRanges) Clone() GRanges {
  result := GRanges{}
  n := r.Length()
  result.Seqnames = make([]string, n)
  result.Ranges   = make([]Range, n)
  result.Strand   = make([]byte, n)
  copy(result.Seqnames, r.Seqnames)
  copy(result.Ranges,   r.Ranges)
  copy(result.Strand,   r.Strand)
  result.Meta = r.Meta.Clone()
  return result
}

/* -------------------------------------------------------------------------- */

func (r GRanges) Length() int {
  return len(r.Ranges)
}

func (r GRanges) Row(i int) GRangesRow {
  return NewGRangesRow(r, i)
}

func (r1 GRanges) Append(r2 GRanges) GRanges {
  result := GRanges{}

  result.Seqnames = append(r1.Seqnames, r2.Seqnames...)
  result.Ranges   = append(r1.Ranges,   r2.Ranges...)
  result.Strand   = append(r1.Strand,   r2.Strand...)

  result.Meta = r1.Meta.Append(r2.Meta)

  return result
}

func (r GRanges) Remove(indices []int) GRanges {
  if len(indices) == 0 {
    return r.Clone()
  }
  indices = removeDuplicatesInt(indices)
  sort.Ints(indices)

  n := r.Length()
  m := n - len(indices)
  // convert indices to subset indices
  idx := make([]int, m)
  for i, j, k := 0, 0, 0; i < r.Length(); i++ {
    for k < len(indices)-1 && i > indices[k] {
      k++
    }
    if i != indices[k] {
      idx[j] = i
      j++
    }
  }
  result := r.Subset(idx)
  result.Meta = r.Meta.Subset(idx)

  return result
}

func (r GRanges) RemoveOverlapsWith(subject GRanges) GRanges {
  queryHits, _ := FindOverlaps(r, subject)
  return r.Remove(queryHits)
}

func (r GRanges) KeepOverlapsWith(subject GRanges) GRanges {
  queryHits, _ := FindOverlaps(r, subject)
  queryHits     = removeDuplicatesInt(queryHits)
  // avoid shuffling the rows
  sort.Ints(queryHits)
  return r.Subset(queryHits)
}

func (r GRanges) Subset(indices []int) GRanges {
  n := len(indices)
  seqnames := make([]string, n)
  from     := make([]int, n)
  to       := make([]int, n)
  strand   := make([]byte, n)

  for i := 0; i < n; i++ {
    seqnames[i] = r.Seqnames[indices[i]]
    from    [i] = r.Ranges  [indices[i]].From
    to      [i] = r.Ranges  [indices[i]].To
    strand  [i] = r.Strand  [indices[i]]
  }
  result := NewGRanges(seqnames, from, to, strand)
  result.Meta = r.Meta.Subset(indices)

  return result
}

func (r GRanges) Slice(ifrom, ito int) GRanges {
  // check parameters
  if ifrom < 0 {
    ifrom = 0
  }
  if ifrom > r.Length() {
    ifrom = r.Length()
  }
  if ito < ifrom {
    ito = ifrom
  }
  if ito > r.Length() {
    ito = r.Length()
  }

  n := ito-ifrom
  seqnames := make([]string, n)
  from     := make([]int, n)
  to       := make([]int, n)
  strand   := make([]byte, n)

  for i := ifrom; i < ito; i++ {
    seqnames[i-ifrom] = r.Seqnames[i]
    from    [i-ifrom] = r.Ranges  [i].From
    to      [i-ifrom] = r.Ranges  [i].To
    strand  [i-ifrom] = r.Strand  [i]
  }
  result := NewGRanges(seqnames, from, to, strand)
  result.Meta = r.Meta.Slice(ifrom, ito)

  return result
}

func (r GRanges) Intersection(s GRanges) GRanges {
  queryHits, subjectHits := FindOverlaps(r, s)
  n        := len(queryHits)
  seqnames := make([]string, n)
  from     := make([]int, n)
  to       := make([]int, n)
  strand   := make([]byte, n)

  for i := 0; i < n; i++ {
    iQ := queryHits[i]
    iS := subjectHits[i]
    gr := r.Ranges[iQ].Intersection(s.Ranges[iS])
    seqnames[i] = r.Seqnames[iQ]
    strand  [i] = r.Strand  [iQ]
    from    [i] = gr.From
    to      [i] = gr.To
  }
  result := NewGRanges(seqnames, from, to, strand)
  result.Meta = r.Meta.Subset(queryHits)

  return result
}

func (r GRanges) Sort(name string, reverse bool) (GRanges, error) {
  j, err := r.sortedIndices(name, reverse)
  if err != nil {
    return GRanges{}, err
  }
  return r.Subset(j), nil
}

// Remove all entries that are not in the given genome.
func (r GRanges) FilterGenome(genome Genome) GRanges {
  idx      := []int{}
  seqnames := make(map[string]int)
  for i := 0; i < genome.Length(); i++ {
    seqnames[genome.Seqnames[i]] = genome.Lengths[i]
  }
  for i := 0; i < r.Length(); i++ {
    length, ok := seqnames[r.Seqnames[i]]
    if ok && r.Ranges[i].To <= length {
      continue
    }
    idx = append(idx, i)
  }
  return r.Remove(idx)
}

func (r GRanges) FilterStrand(s byte) GRanges {
  idx := []int{}
  for i := 0; i < r.Length(); i++ {
    if r.Strand[i] != s {
      idx = append(idx, i)
    }
  }
  return r.Remove(idx)
}

// Set length of each range to the given value. Ranges
// with no strand information are not changed.
func (r GRanges) SetLengths(n int) GRanges {
  s := r.Clone()
  // negative values are not allowed
  if n < 0 {
    n = 0
  }
  for i := 0; i < s.Length(); i++ {
    // forward strand
    if s.Strand[i] == '+' {
      s.Ranges[i].To = s.Ranges[i].From+n
    }
    // reverse strand
    if s.Strand[i] == '-' {
      s.Ranges[i].From = s.Ranges[i].To-n
    }
    // if strand is '*' do nothing
  }
  return s
}

// Add data from a track to the GRanges object. The data will be
// contained in a meta-data column with the same name as the track.
// It is required that each range has the same length.
func (r *GRanges) ImportTrack(track Track, revNegStrand bool) (*GRanges, error) {
  n := r.Length()
  m := -1
  data    := make([][]float64, n)
  binsize := make([]int, n)
  // fill matrix
  for i := 0; i < n; i++ {
    from := r.Ranges[i].From
    to   := r.Ranges[i].To
    seq  := r.Seqnames[i]
    if m == -1 {
      m = divIntUp(to - from, track.Binsize)
    } else if m != divIntUp(to - from, track.Binsize) {
      return nil, fmt.Errorf("varying window sizes are not allowed")
    }
    // all rows are using the same binsize
    binsize[i] = track.Binsize
    // loop over window
    data[i] = make([]float64, m)
    if r.Strand[i] == '+' || revNegStrand == false {
      for j, k := 0, from; k < to; k, j = k+track.Binsize, j+1 {
        value, err := track.At(seq, k)
        if err == nil {
          data[i][j] = value
        }
      }
    } else if r.Strand[i] == '-' {
      for j, k := 0, to-1; k >= from; k, j = k-track.Binsize, j+1 {
        value, err := track.At(seq, k)
        if err == nil {
          data[i][j] = value
        }
      }
    } else {
      return nil, fmt.Errorf("range has no strand information")
    }
  }
  r.AddMeta("binsize", binsize)
  r.AddMeta(track.Name, data)
  return r, nil
}

func (r *GRanges) ImportBigWig(reader *BigWigReader, binsize int, revNegStrand bool) (*GRanges, error) {
  counts := [][]float64{}
  for i := 0; i < r.Length(); i++ {
    seqname := r.Seqnames[i]
    from    := r.Ranges[i].From
    to      := r.Ranges[i].To
    strand  := r.Strand[i]
    n   := (to - from)/binsize
    seq := make([]float64, n)
    for record := range reader.Query(seqname, from, to, binsize) {
      if record.Error != nil {
        return nil, fmt.Errorf("%v", record.Error)
      }
      if revNegStrand == false || strand == '+' {
        if idx := (record.From - from)/binsize; idx >= 0 && idx < n {
          seq[idx] = record.Sum/float64(record.Valid)
        }
      } else {
        if idx := (record.From - from)/binsize; idx >= 0 && idx < n {
          seq[n-1-idx] = record.Sum/float64(record.Valid)
        }
      }
    }
    counts = append(counts, seq)
  }
  r.AddMeta("counts", counts)
  return r, nil
}

/* convert to gene object
 * -------------------------------------------------------------------------- */

func (granges GRanges) Genes() Genes {
  return newGenes(granges)
}

/* convert to string
 * -------------------------------------------------------------------------- */

func (granges GRanges) String() string {
  return granges.PrettyPrint(10)
}
