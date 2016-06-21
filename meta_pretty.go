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
import "bytes"
import "bufio"
import "io"
import "io/ioutil"

/* -------------------------------------------------------------------------- */

func (meta Meta) PrettyPrint(n int) string {
  var buffer bytes.Buffer
  writer := bufio.NewWriter(&buffer)

  printCellSlice := func(writer io.Writer, widths []int, i, j int, data interface{}) int {
    var tmpBuffer bytes.Buffer
    tmpWriter := bufio.NewWriter(&tmpBuffer)
    switch v := data.(type) {
    case [][]string:
      for k := 0; k < len(v[i]); k++ {
        fmt.Fprintf(tmpWriter, " %s", v[i][k])
      }
    case [][]float64:
      for k := 0; k < len(v[i]); k++ {
        fmt.Fprintf(tmpWriter, " %f", v[i][k])
      }
    case [][]int:
      for k := 0; k < len(v[i]); k++ {
        fmt.Fprintf(tmpWriter, " %d", v[i][k])
      }
    default:
      panic("invalid meta data")
    }
    tmpWriter.Flush()
    format := fmt.Sprintf(" %%%ds", widths[j]-1)
    l, _ := fmt.Fprintf(writer, format, tmpBuffer.String())
    return l
  }
  printCell := func(writer io.Writer, widths []int, i, j int) int {
    length := 0
    switch v := meta.MetaData[j].(type) {
    case []string :
      format := fmt.Sprintf(" %%%ds", widths[j]-1)
      length, _ = fmt.Fprintf(writer, format, v[i])
    case []float64:
      format := fmt.Sprintf(" %%%df", widths[j]-1)
      length, _ = fmt.Fprintf(writer, format, v[i])
    case []int    :
      format := fmt.Sprintf(" %%%dd", widths[j]-1)
      length, _ = fmt.Fprintf(writer, format, v[i])
    default:
      length += printCellSlice(writer, widths, i, j, v)
    }
    return length
  }
  printRow := func(writer io.Writer, widths []int, i int) {
    if i != 0 {
      fmt.Fprintf(writer, "\n")
    }
    for j := 0; j < meta.MetaLength(); j++ {
      printCell(writer, widths, i, j)
    }
  }
  // compute widths of all cells in row i
  updateMaxWidths := func(i int, widths []int) {
    for j := 0; j < meta.MetaLength(); j++ {
      width := printCell(ioutil.Discard, widths, i, j)
      if width > widths[j] {
        widths[j] = width
      }
    }
  }
  printHeader := func(writer io.Writer, widths []int) {
    for j := 0; j < meta.MetaLength(); j++ {
      format := fmt.Sprintf(" %%%ds", widths[j]-1)
      fmt.Fprintf(writer, format, meta.MetaName[j])
    }
    fmt.Fprintf(writer, "\n")
  }
  applyRows := func(f1 func(i int), f2 func()) {
    if meta.Length() <= n+1 {
      // apply to all entries
      for i := 0; i < meta.Length(); i++ { f1(i) }
    } else {
      // apply to first n/2 rows
      for i := 0; i < n/2; i++ { f1(i) }
      // between first and last n/2 rows
      f2()
      // apply to last n/2 rows
      for i := meta.Length() - n/2; i < meta.Length(); i++ { f1(i) }
    }
  }
  // maximum column widths
  widths := make([]int, meta.MetaLength())
  for j := 0; j < meta.MetaLength(); j++ {
    widths[j], _ = fmt.Fprintf(ioutil.Discard, " %s", meta.MetaName[j])
  }
  // determine column widths
  applyRows(func(i int) { updateMaxWidths(i, widths) }, func() {})
  // pring header
  printHeader(writer, widths)
  // print rows
  applyRows(
    func(i int) { printRow(writer, widths, i) },
    func() {
      fmt.Fprintf(writer, "\n")
        for j := 0; j < meta.MetaLength(); j++ {
          format := fmt.Sprintf(" %%%ds", widths[j]-1)
          fmt.Fprintf(writer, format, "...")
        }
    })
  writer.Flush()

  return buffer.String()
}
