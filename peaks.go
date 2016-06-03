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

//import "fmt"
import "bufio"
import "os"
import "strconv"
import "strings"

import . "github.com/pbenner/pshape/Utility"

/* constructors
 * -------------------------------------------------------------------------- */

func NewGPeaks(seqnames []string, from, to, absSummit []int, pileup, pvalue, foldEnrichment, qvalue []float64) GRanges {
  r := NewGRanges(seqnames, from, to, []byte{})
  r.AddMeta("abs_summit",      absSummit)
  r.AddMeta("pileup",          pileup)
  r.AddMeta("-log10(pvalue)",  pvalue)
  r.AddMeta("fold_enrichment", foldEnrichment)
  r.AddMeta("-log10(qvalue)",  qvalue)
  return r
}

/* i/o
 * -------------------------------------------------------------------------- */

func ReadXlsPeaks(filename string) GRanges {

  f, err := os.Open(filename)
  Check(err)
  defer f.Close()
  // check if we already saw the header
  header := false

  seqnames       := []string{}
  from           := []int{}
  to             := []int{}
  absSummit      := []int{}
  pileup         := []float64{}
  pvalue         := []float64{}
  foldEnrichment := []float64{}
  qvalue         := []float64{}

  scanner := bufio.NewScanner(f)
  for scanner.Scan() {
    fields := strings.Fields(scanner.Text())
    if len(fields) == 0 {
      continue
    }
    if fields[0] == "#" {
      continue
    }
    if len(fields) != 10 {
      panic("Invalid peaks xls file!")
    }
    if header == false {
      // first line is the header
      if fields[0] != "chr"             || fields[1] != "start"  || fields[2] != "end" ||
        (fields[4] != "abs_summit"      || fields[5] != "pileup" || fields[6] != "-log10(pvalue)") ||
        (fields[7] != "fold_enrichment" || fields[8] != "-log10(qvalue)") {
        panic("Invalid Xls header!")
      }
      header = true
      continue
    }
    t1, e1 := strconv.ParseInt(fields[1], 10, 64) // from
    t2, e2 := strconv.ParseInt(fields[2], 10, 64) // to
    t3, e3 := strconv.ParseInt(fields[4], 10, 64) // abs_summit
    t4, e4 := strconv.ParseFloat(fields[5], 64)   // pileup
    t5, e5 := strconv.ParseFloat(fields[6], 64)   // pvalue
    t6, e6 := strconv.ParseFloat(fields[7], 64)   // fold_enrichment
    t7, e7 := strconv.ParseFloat(fields[8], 64)   // qvalue
    Check(e1); Check(e2); Check(e3); Check(e4); Check(e5); Check(e6); Check(e7)

    seqnames       = append(seqnames,       fields[0])
    from           = append(from,           int(t1))
    to             = append(to,             int(t2)+1)
    absSummit      = append(absSummit,      int(t3))
    pileup         = append(pileup,         float64(t4))
    pvalue         = append(pvalue,         float64(t5))
    foldEnrichment = append(foldEnrichment, float64(t6))
    qvalue         = append(qvalue,         float64(t7))
  }
  return NewGPeaks(seqnames, from, to, absSummit, pileup, pvalue, foldEnrichment, qvalue)
}
