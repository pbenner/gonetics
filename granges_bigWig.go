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
  counts := [][]float64{}
  for i := 0; i < r.Length(); i++ {
    seqname := r.Seqnames[i]
    from    := r.Ranges[i].From
    to      := r.Ranges[i].To
    strand  := r.Strand[i]
    n := 0
    r := []BbiSummaryRecord{}
    s := []float64{}
    // first collect all records
    for record := range bwr.Query(seqname, from, to, binSize) {
      if binSize == 0 {
        binSize = record.To - record.From
      }
      if len(r) == 0 {
        n = (to - from)/binSize
        r = make([]BbiSummaryRecord, n)
        s = make([]float64, n)
      }
      if record.Error != nil {
        return fmt.Errorf("%v", record.Error)
      }
      if revNegStrand == false || strand == '+' {
        if idx := (record.From - from)/binSize; idx >= 0 && idx < n {
          r[idx] = record.BbiSummaryRecord
        }
      } else {
        if idx := (record.From - from)/binSize; idx >= 0 && idx < n {
          r[n-1-idx] = record.BbiSummaryRecord
        }
      }
    }
    // convert summary records to sequence
    if init != 0.0 {
      for i := 0; i < len(s); i++ {
        s[i] = init
      }
    }
    if binOverlap != 0 {
      t := BbiSummaryRecord{}
      for i := 0; i < len(s); i++ {
        t.Reset()
        for j := i-binOverlap; j <= i+binOverlap; j++ {
          if j < 0 || j >= len(s) {
            continue
          }
          t.AddRecord(r[j])
        }
        if t.Valid > 0 {
          s[i] = f(t.Sum, t.SumSquares, t.Min, t.Max, t.Valid)
        }
      }
    } else {
      for i, t := range r {
        if t.Valid > 0 {
          s[i] = f(t.Sum, t.SumSquares, t.Min, t.Max, t.Valid)
        }
      }
    }
    counts = append(counts, s)
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
