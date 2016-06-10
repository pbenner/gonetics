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

/* -------------------------------------------------------------------------- */

func (genes Genes) WriteTable(filename string, header, compress bool) {
  var buffer bytes.Buffer

  w := bufio.NewWriter(&buffer)

  // print header
  if header {
    fmt.Fprintf(w, "%16s %10s %6s %10s %10s %10s %10s",
      "names", "seqnames", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd")
    genes.Meta.WriteTableRow(w, -1)
    fmt.Fprintf(w, "\n")
  }
  // print data
  for i := 0; i < genes.Length(); i++ {
    fmt.Fprintf(w,  "%16s", genes.Names[i])
    fmt.Fprintf(w, " %10s", genes.Seqnames[i])
    fmt.Fprintf(w, " %6c",  genes.Strand[i])
    fmt.Fprintf(w, " %10d", genes.Tx[i].From)
    fmt.Fprintf(w, " %10d", genes.Tx[i].To)
    fmt.Fprintf(w, " %10d", genes.Cds[i].From)
    fmt.Fprintf(w, " %10d", genes.Cds[i].To)
    genes.Meta.WriteTableRow(w, i)
    fmt.Fprintf(w, "\n")
  }
  w.Flush()
  writeFile(filename, &buffer, compress)
}
