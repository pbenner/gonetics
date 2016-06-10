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
import "bytes"
import "fmt"
import "io"
import "io/ioutil"
import "strings"

/* -------------------------------------------------------------------------- */

func (genes Genes) PrettyPrint(n int) string {
  var buffer bytes.Buffer
  writer := bufio.NewWriter(&buffer)

  // pretty print meta data and create a scanner reading
  // the resulting string
  metaStr     := genes.Meta.PrettyPrint(n)
  metaReader  := strings.NewReader(metaStr)
  metaScanner := bufio.NewScanner(metaReader)

  // compute the width of a single cell
  updateMaxWidth := func(format string, widths []int, j int, args ...interface{}) {
    width, _ := fmt.Fprintf(ioutil.Discard, format, args...)
    if width > widths[j] {
      widths[j] = width
    }
  }
  // compute widths of all cells in row i
  updateMaxWidths := func(i int, widths []int) {
    updateMaxWidth("%d", widths, 0, i+1)
    updateMaxWidth("%s", widths, 1, genes.Names[i])
    updateMaxWidth("%s", widths, 2, genes.Seqnames[i])
    updateMaxWidth("%d", widths, 3, genes.Tx[i].From)
    updateMaxWidth("%d", widths, 4, genes.Tx[i].To)
    updateMaxWidth("%d", widths, 5, genes.Cds[i].From)
    updateMaxWidth("%d", widths, 6, genes.Cds[i].To)
    updateMaxWidth("%c", widths, 7, genes.Strand[i])
  }
  printMetaRow := func(writer io.Writer) {
    if genes.MetaLength() != 0 {
      fmt.Fprintf(writer, " | ")
      metaScanner.Scan()
      fmt.Fprintf(writer, "%s", metaScanner.Text())
    }
  }
  printHeader := func(writer io.Writer, format string) {
    fmt.Fprintf(writer, format,
      "", "names", "seqnames", "transcripts", "cds", "strand")
    printMetaRow(writer)
    fmt.Fprintf(writer, "\n")
  }
  printRow := func(writer io.Writer, format string, i int) {
    if i != 0 {
      fmt.Fprintf(writer, "\n")
    }
    fmt.Fprintf(writer, format,
      i+1,
      genes.Names[i],
      genes.Seqnames[i],
      genes.Tx[i].From,
      genes.Tx[i].To,
      genes.Cds[i].From,
      genes.Cds[i].To,
      genes.Strand[i])
    printMetaRow(writer)
  }
  applyRows := func(f1 func(i int), f2 func()) {
    if genes.Length() <= n+1 {
      // apply to all entries
      for i := 0; i < genes.Length(); i++ { f1(i) }
    } else {
      // apply to first n/2 rows
      for i := 0; i < n/2; i++ { f1(i) }
      // between first and last n/2 rows
      f2()
      // apply to last n/2 rows
      for i := genes.Length() - n/2; i < genes.Length(); i++ { f1(i) }
    }
  }
  // maximum column widths
  widths := []int{1, 5, 8, 1, 1, 1, 1, 6}
  // determine column widths
  applyRows(func(i int) { updateMaxWidths(i, widths) }, func() {})
  // generate format strings
  formatRow    := fmt.Sprintf("%%%dd %%%ds %%%ds [%%%dd, %%%dd) [%%%dd, %%%dd) %%%dc",
    widths[0], widths[1], widths[2], widths[3], widths[4], widths[5], widths[6], widths[7])
  formatHeader := fmt.Sprintf("%%%ds %%%ds %%%ds %%%ds %%%ds %%%ds",
    widths[0], widths[1], widths[2], widths[3]+widths[4]+4, widths[5]+widths[6]+4, widths[7])
  // pring header
  printHeader(writer, formatHeader)
  // print rows
  applyRows(
    func(i int) {
      printRow(writer, formatRow, i)
    },
    func() {
      fmt.Fprintf(writer, "\n")
      fmt.Fprintf(writer, formatHeader, "", "...", "...", "...", "...", "")
      printMetaRow(writer)
    })
  writer.Flush()

  return buffer.String()
}
