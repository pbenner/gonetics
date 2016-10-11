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

func mergeDuplicates(r GRanges) GRanges {
  type grline struct {
    seqname  string
    from     int
    to       int
  }
  m := make(map[grline][]int)
  for i := 0; i < r.Length(); i++ {
    t := grline{r.Seqnames[i], r.Ranges[i].From, r.Ranges[i].To}
    m[t] = append(m[t], i)
  }
  // get first index for filtering
  filterIdx := []int{}
  // map first index to all rows with the same ranges entry
  filterMap := make(map[int][]int)
  // assign the same value to all rows that should be merged
  mergeIdx  := make([]int, r.Length())
  for _, indices := range m {
    filterIdx = append(filterIdx, indices[0])
    filterMap[indices[0]] = indices
  }
  sort.Ints(filterIdx)
  for i, firstIdx := range filterIdx {
    for _, idx := range filterMap[firstIdx] {
      mergeIdx[idx] = i
    }
  }
  result := r.Subset(filterIdx)
  result.Meta = r.Meta.Merge(mergeIdx)
  return result
}

// Generate a GRanges object where each range is a window around
// the transcription start site (TSS) of a gene. The arguments offset1,
// and offset2 determine the size of the window starting from the TSS
// in 5' and 3' direction respectively. Genomic bounds are not
// checked. (This is required by the GRanges ImportTrack() method!)
func Promoters(genes Genes, offset1, offset2 int) (GRanges, error) {
  result := GRanges{}
  names  := []string{}
  // fill matrix
  for i := 0; i < genes.Length(); i++ {
    from, to := 0,0
    if genes.Strand[i] == '+' {
      from = genes.Ranges[i].From - offset1
      to   = genes.Ranges[i].From + offset2 + 1
    } else if genes.Strand[i] == '-' {
      from = genes.Ranges[i].To - 1 - offset2
      to   = genes.Ranges[i].To - 1 + offset1 + 1
    } else {
      return result, fmt.Errorf("gene `%d' has no strand information", i)
    }
    result.Seqnames = append(result.Seqnames, genes.Seqnames[i])
    result.Ranges   = append(result.Ranges,   NewRange(from, to))
    result.Strand   = append(result.Strand,   genes.Strand[i])
    // name of the region
    names = append(names, genes.Names[i])
  }
  // add names as a meta column
  result.Meta = genes.Meta.Clone()
  result.DeleteMeta("cds")
  result.AddMeta("names", names)
  return mergeDuplicates(result), nil
}
