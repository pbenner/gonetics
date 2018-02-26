/* Copyright (C) 2016-2017 Philipp Benner
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
import "os"

/* -------------------------------------------------------------------------- */

func (r *GRanges) ReadBigWig(bwr *BigWigReader, name string, f BinSummaryStatistics, binSize, binOverlap int, init float64, revNegStrand bool) error {
  rev := func(s []float64) []float64 {
    for i := 0; i < len(s)/2; i++ {
      s[i], s[len(s)-1-i] = s[len(s)-1-i], s[i]
    }
    return s
  }
  counts := [][]float64{}
  for i := 0; i < r.Length(); i++ {
    seqname := r.Seqnames[i]
    from    := r.Ranges[i].From
    to      := r.Ranges[i].To
    strand  := r.Strand[i]
    if s, bs, err := bwr.QuerySlice(seqname, from, to, f, binSize, binOverlap, init); err != nil {
      return err
    } else {
      if binSize == 0 {
        binSize = bs
      }
      if revNegStrand == false || strand == '+' {
        counts = append(counts, s)
      } else {
        counts = append(counts, rev(s))
      }
    }
  }
  r.AddMeta(name, counts)
  return nil
}

func (r *GRanges) ImportBigWig(filename string, name string, s BinSummaryStatistics, binSize, binOverlap int, init float64, revNegStrand bool) error {
  f, err := os.Open(filename)
  if err != nil {
    return err
  }
  defer f.Close()

  bwr, err := NewBigWigReader(f)
  if err != nil {
    return err
  }
  if err := r.ReadBigWig(bwr, name, s, binSize, binOverlap, init, revNegStrand); err != nil {
    return fmt.Errorf("importing bigWig file from `%s' failed: %v", filename, err)
  }
  return nil
}
