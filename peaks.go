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
import "bufio"
import "os"
import "io"
import "strconv"
import "strings"

/* -------------------------------------------------------------------------- */

type GPeaks struct {
  GRanges
}

/* constructors
 * -------------------------------------------------------------------------- */

func NewGPeaks(seqnames []string, from, to, absSummit []int, pileup, pvalue, foldEnrichment, qvalue []float64) GPeaks {
  r := NewGRanges(seqnames, from, to, []byte{})
  r.AddMeta("abs_summit",      absSummit)
  r.AddMeta("pileup",          pileup)
  r.AddMeta("-log10(pvalue)",  pvalue)
  r.AddMeta("fold_enrichment", foldEnrichment)
  r.AddMeta("-log10(qvalue)",  qvalue)
  return GPeaks{r}
}

/* -------------------------------------------------------------------------- */

func (p GPeaks) AbsSummit() []int {
  return p.GetMetaInt("abs_summit")
}

func (p GPeaks) Pileup() []float64 {
  return p.GetMetaFloat("pileup")
}

func (p GPeaks) Pvalue() []float64 {
  return p.GetMetaFloat("-log10(pvalue)")
}

func (p GPeaks) FoldEnrichment() []float64 {
  return p.GetMetaFloat("fold_enrichment")
}

func (p GPeaks) Qvalue() []float64 {
  return p.GetMetaFloat("-log10(qvalue)")
}

/* i/o
 * -------------------------------------------------------------------------- */

func ReadXlsPeaks(r io.Reader) (GPeaks, error) {

  var gpeaks GPeaks

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

  scanner := bufio.NewScanner(r)
  for scanner.Scan() {
    fields := strings.Fields(scanner.Text())
    if len(fields) == 0 {
      continue
    }
    if fields[0] == "#" {
      continue
    }
    if len(fields) != 10 {
      return gpeaks, fmt.Errorf("invalid peaks xls file")
    }
    if header == false {
      // first line is the header
      if fields[0] != "chr"             || fields[1] != "start"  || fields[2] != "end" ||
        (fields[4] != "abs_summit"      || fields[5] != "pileup" || fields[6] != "-log10(pvalue)") ||
        (fields[7] != "fold_enrichment" || fields[8] != "-log10(qvalue)") {
        return gpeaks, fmt.Errorf("Invalid Xls header!")
      }
      header = true
      continue
    }
    t1, e := strconv.ParseInt(fields[1], 10, 64) // from
    if e != nil {
      return gpeaks, e
    }
    t2, e := strconv.ParseInt(fields[2], 10, 64) // to
    if e != nil {
      return gpeaks, e
    }
    t3, e := strconv.ParseInt(fields[4], 10, 64) // abs_summit
    if e != nil {
      return gpeaks, e
    }
    t4, e := strconv.ParseFloat(fields[5], 64)   // pileup
    if e != nil {
      return gpeaks, e
    }
    t5, e := strconv.ParseFloat(fields[6], 64)   // pvalue
    if e != nil {
      return gpeaks, e
    }
    t6, e := strconv.ParseFloat(fields[7], 64)   // fold_enrichment
    if e != nil {
      return gpeaks, e
    }
    t7, e := strconv.ParseFloat(fields[8], 64)   // qvalue
    if e != nil {
      return gpeaks, e
    }


    seqnames       = append(seqnames,       fields[0])
    from           = append(from,           int(t1))
    to             = append(to,             int(t2)+1)
    absSummit      = append(absSummit,      int(t3))
    pileup         = append(pileup,         float64(t4))
    pvalue         = append(pvalue,         float64(t5))
    foldEnrichment = append(foldEnrichment, float64(t6))
    qvalue         = append(qvalue,         float64(t7))
  }
  gpeaks = NewGPeaks(seqnames, from, to, absSummit, pileup, pvalue, foldEnrichment, qvalue)

  return gpeaks, nil
}

func ImportXlsPeaks(filename string) (GPeaks, error) {
  f, err := os.Open(filename)
  if err != nil {
    return GPeaks{}, err
  }
  defer f.Close()

  return ReadXlsPeaks(f)
}
