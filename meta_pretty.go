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

func (meta Meta) WritePretty(writer io.Writer, n int, args ...interface{}) error {
  useScientific := false
  for _, arg := range args {
    switch a := arg.(type) {
    case OptionPrintScientific:
      useScientific = a.Value
    default:
    }
  }
  printCellSlice := func(writer io.Writer, widths []int, i, j int, data interface{}) (int, error) {
    var tmpBuffer bytes.Buffer
    tmpWriter := bufio.NewWriter(&tmpBuffer)
    switch v := data.(type) {
    case [][]string:
      for k := 0; k < len(v[i]); k++ {
        if _, err := fmt.Fprintf(tmpWriter, " %s", v[i][k]); err != nil {
          return 0, err
        }
      }
    case [][]float64:
      if useScientific {
        for k := 0; k < len(v[i]); k++ {
          if _, err := fmt.Fprintf(tmpWriter, " %e", v[i][k]); err != nil {
            return 0, err
          }
        }
      } else {
        for k := 0; k < len(v[i]); k++ {
          if _, err := fmt.Fprintf(tmpWriter, " %f", v[i][k]); err != nil {
            return 0, err
          }
        }
      }
    case [][]int:
      for k := 0; k < len(v[i]); k++ {
        if _, err := fmt.Fprintf(tmpWriter, " %d", v[i][k]); err != nil {
          return 0, err
        }
      }
    default:
      panic("invalid meta data")
    }
    tmpWriter.Flush()
    format := fmt.Sprintf(" %%%ds", widths[j]-1)
    l, _ := fmt.Fprintf(writer, format, tmpBuffer.String())
    return l, nil
  }
  printCell := func(writer io.Writer, widths []int, i, j int) (int, error) {
    switch v := meta.MetaData[j].(type) {
    case []string:
      format := fmt.Sprintf(" %%%ds", widths[j]-1)
      return fmt.Fprintf(writer, format, v[i])
    case []float64:
      if useScientific {
        format := fmt.Sprintf(" %%%de", widths[j]-1)
        return fmt.Fprintf(writer, format, v[i])
      } else {
        format := fmt.Sprintf(" %%%df", widths[j]-1)
        return fmt.Fprintf(writer, format, v[i])
      }
    case []int:
      format := fmt.Sprintf(" %%%dd", widths[j]-1)
      return fmt.Fprintf(writer, format, v[i])
    case []Range:
      format := fmt.Sprintf(" [ %%d, %%d ]")
      return fmt.Fprintf(writer, format, v[i].From, v[i].To)
    default:
      return printCellSlice(writer, widths, i, j, v)
    }
  }
  printRow := func(writer io.Writer, widths []int, i int) error {
    if i != 0 {
      if _, err := fmt.Fprintf(writer, "\n"); err != nil {
        return err
      }
    }
    for j := 0; j < meta.MetaLength(); j++ {
      if _, err := printCell(writer, widths, i, j); err != nil {
        return err
      }
    }
    return nil
  }
  // compute widths of all cells in row i
  updateMaxWidths := func(i int, widths []int) error {
    for j := 0; j < meta.MetaLength(); j++ {
      if width, err := printCell(ioutil.Discard, widths, i, j); err != nil {
        return err
      } else {
        if width > widths[j] {
          widths[j] = width
        }
      }
    }
    return nil
  }
  printHeader := func(writer io.Writer, widths []int) error {
    for j := 0; j < meta.MetaLength(); j++ {
      format := fmt.Sprintf(" %%%ds", widths[j]-1)
      if _, err := fmt.Fprintf(writer, format, meta.MetaName[j]); err != nil {
        return err
      }
    }
    if _, err := fmt.Fprintf(writer, "\n"); err != nil {
      return err
    }
    return nil
  }
  applyRows := func(f1 func(i int) error, f2 func() error) error {
    if meta.Length() <= n+1 {
      // apply to all entries
      for i := 0; i < meta.Length(); i++ {
        if err := f1(i); err != nil {
          return err
        }
      }
    } else {
      // apply to first n/2 rows
      for i := 0; i < n/2; i++ {
        if err := f1(i); err != nil {
          return err
        }
      }
      // between first and last n/2 rows
      if err := f2(); err != nil {
        return err
      }
      // apply to last n/2 rows
      for i := meta.Length() - n/2; i < meta.Length(); i++ {
        if err := f1(i); err != nil {
          return err
        }
      }
    }
    return nil
  }
  // maximum column widths
  widths := make([]int, meta.MetaLength())
  for j := 0; j < meta.MetaLength(); j++ {
    if width, err := fmt.Fprintf(ioutil.Discard, " %s", meta.MetaName[j]); err != nil {
      return err
    } else {
      widths[j] = width
    }
  }
  // determine column widths
  if err := applyRows(func(i int) error { return updateMaxWidths(i, widths) }, func() error { return nil}); err != nil {
    return err
  }
  // pring header
  if err := printHeader(writer, widths); err != nil {
    return err
  }
  // print rows
  if err := applyRows(
    func(i int) error { return printRow(writer, widths, i) },
    func() error {
      if _, err := fmt.Fprintf(writer, "\n"); err != nil {
        return err
      }
      for j := 0; j < meta.MetaLength(); j++ {
        format := fmt.Sprintf(" %%%ds", widths[j]-1)
        if _, err := fmt.Fprintf(writer, format, "..."); err != nil {
          return err
        }
      }
      return nil
    }); err != nil {
    return err
  }
  return nil
}

func (meta *Meta) PrintPretty(n int, args ...interface{}) string {
  var buffer bytes.Buffer
  writer := bufio.NewWriter(&buffer)

  if err := meta.WritePretty(writer, n, args...); err != nil {
    return ""
  }
  writer.Flush()

  return buffer.String()
}
