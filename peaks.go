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

package bioinf

/* -------------------------------------------------------------------------- */

import "bufio"
import "bytes"
import "fmt"
import "os"
import "strconv"
import "strings"

import . "github.com/pbenner/pshape/Utility"

/* -------------------------------------------------------------------------- */

type GPeaks struct {
  GRanges
  AbsSummit      []int
  Pvalue         []float64
  FoldEnrichment []float64
}

/* constructors
 * -------------------------------------------------------------------------- */

func NewGPeaks(seqnames []string, from, to, absSummit []int, pvalue, foldEnrichment []float64) GPeaks {
  r := NewGRanges(seqnames, from, to, []byte{})
  n := r.Length()
  if len(absSummit) != n || len(pvalue) != n || len(foldEnrichment) != n {
    panic("NewGPeaks(): invalid arguments!")
  }
  return GPeaks{r, absSummit, pvalue, foldEnrichment}
}

/* convert to string
 * -------------------------------------------------------------------------- */

func (gpeaks GPeaks) String() string {
  var buffer bytes.Buffer
  // number of lines to print
  const n int = 10

  printRow := func(i int) {
    if i != 0 {
      buffer.WriteString("\n")
    }
    buffer.WriteString(
      fmt.Sprintf("%10d %10s [%10d, %10d] | %10d %14f %15f",
        i+1,
        gpeaks.Seqnames[i],
        gpeaks.Ranges[i].From,
        gpeaks.Ranges[i].To,
        gpeaks.AbsSummit[i],
        gpeaks.Pvalue[i],
        gpeaks.FoldEnrichment[i]))
  }

  // pring header
  buffer.WriteString(
    fmt.Sprintf("%10s %10s %24s | %10s %14s %15s\n",
      "", "seqnames", "ranges",
      "abs_summit", "-log10(pvalue)", "fold_enrichment"))

  // select rows to print
  if gpeaks.Length() <= n+1 {
    // print all entries
    for i := 0; i < gpeaks.Length(); i++ {
      printRow(i)
    }
  } else {
    // print first n/2 rows
    for i := 0; i < n/2; i++ {
      printRow(i)
    }
    buffer.WriteString(
      fmt.Sprintf("\n%10s %10s %24s | %10s", "", "...", "...", "..."))
    // print last n/2 rows
    for i := gpeaks.Length() - n/2; i < gpeaks.Length(); i++ {
      printRow(i)
    }
  }

  return buffer.String()
}

/* i/o
 * -------------------------------------------------------------------------- */

func ReadXlsPeaks(filename string) GPeaks {

  f, err := os.Open(filename)
  Check(err)
  defer f.Close()
  // check if we already saw the header
  header := false

  seqnames       := []string{}
  from           := []int{}
  to             := []int{}
  absSummit      := []int{}
  pvalue         := []float64{}
  foldEnrichment := []float64{}

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
      if fields[0] != "chr"        || fields[1] != "start"          || fields[2] != "end" ||
        (fields[4] != "abs_summit" || fields[6] != "-log10(pvalue)" || fields[7] != "fold_enrichment") {
        panic("Invalid Xls header!")
      }
      header = true
      continue
    }
    t1, e1 := strconv.ParseInt(fields[1], 10, 64) // from
    t2, e2 := strconv.ParseInt(fields[2], 10, 64) // to
    t3, e3 := strconv.ParseInt(fields[4], 10, 64) // abs_summit
    t4, e4 := strconv.ParseFloat(fields[6], 64)   // pvalue
    t5, e5 := strconv.ParseFloat(fields[7], 64)   // fold_enrichment
    Check(e1); Check(e2); Check(e3); Check(e4); Check(e5)

    seqnames       = append(seqnames,       fields[0])
    from           = append(from,           int(t1))
    to             = append(to,             int(t2))
    absSummit      = append(absSummit,      int(t3))
    pvalue         = append(pvalue,         float64(t4))
    foldEnrichment = append(foldEnrichment, float64(t5))
  }
  return NewGPeaks(seqnames, from, to, absSummit, pvalue, foldEnrichment)
}
