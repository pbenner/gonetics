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
import "bytes"
import "io"
import "io/ioutil"
import "strings"

/* -------------------------------------------------------------------------- */

func (granges GRanges) WritePretty(writer io.Writer, n int) error {
  // pretty print meta data and create a scanner reading
  // the resulting string
  metaStr     := granges.Meta.PrintPretty(n)
  metaReader  := strings.NewReader(metaStr)
  metaScanner := bufio.NewScanner(metaReader)

  // compute the width of a single cell
  updateMaxWidth := func(format string, widths []int, j int, args ...interface{}) error {
    width, err := fmt.Fprintf(ioutil.Discard, format, args...)
    if err != nil {
      return err
    }
    if width > widths[j] {
      widths[j] = width
    }
    return nil
  }
  // compute widths of all cells in row i
  updateMaxWidths := func(i int, widths []int) error {
    if err := updateMaxWidth("%d", widths, 0, i+1); err != nil {
      return err
    }
    if err := updateMaxWidth("%s", widths, 1, granges.Seqnames[i]); err != nil {
      return err
    }
    if err := updateMaxWidth("%d", widths, 2, granges.Ranges[i].From); err != nil {
      return err
    }
    if err := updateMaxWidth("%d", widths, 3, granges.Ranges[i].To); err != nil {
      return err
    }
    if err := updateMaxWidth("%c", widths, 4, granges.Strand[i]); err != nil {
      return err
    }
    return nil
  }
  printMetaRow := func(writer io.Writer) error {
    if granges.MetaLength() != 0 {
      if _, err := fmt.Fprintf(writer, " | "); err != nil {
        return err
      }
      metaScanner.Scan()
      if _, err := fmt.Fprintf(writer, "%s", metaScanner.Text()); err != nil {
        return err
      }
    }
    return nil
  }
  printHeader := func(writer io.Writer, format string) error {
    if _, err := fmt.Fprintf(writer, format, "", "seqnames", "ranges", "strand"); err != nil {
      return err
    }
    return printMetaRow(writer)
  }
  printRow := func(writer io.Writer, format string, i int) error {
    if _, err := fmt.Fprintf(writer, "\n"); err != nil {
      return err
    }
    if _, err := fmt.Fprintf(writer, format,
      i+1,
      granges.Seqnames[i],
      granges.Ranges[i].From,
      granges.Ranges[i].To,
      granges.Strand[i]); err != nil {
      return err
    }
    return printMetaRow(writer)
  }
  applyRows := func(f1 func(i int) error, f2 func() error) error {
    if granges.Length() <= n+1 {
      // apply to all entries
      for i := 0; i < granges.Length(); i++ {
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
      for i := granges.Length() - n/2; i < granges.Length(); i++ {
        if err := f1(i); err != nil {
          return err
        }
      }
    }
    return nil
  }
  // maximum column widths
  widths := []int{1, 8, 1, 1, 6}
  // determine column widths
  if err := applyRows(func(i int) error { return updateMaxWidths(i, widths) }, func() error { return nil }); err != nil {
    return err
  }
  // generate format strings
  formatRow    := fmt.Sprintf("%%%dd %%%ds [%%%dd, %%%dd) %%%dc",
    widths[0], widths[1], widths[2], widths[3], widths[4])
  formatHeader := fmt.Sprintf("%%%ds %%%ds %%%ds %%%ds",
    widths[0], widths[1], widths[2]+widths[3]+4, widths[4])
  // pring header
  if err := printHeader(writer, formatHeader); err != nil {
    return err
  }
  // print rows
  if err := applyRows(
    func(i int) error {
      return printRow(writer, formatRow, i)
    },
    func() error {
      if _, err := fmt.Fprintf(writer, "\n"); err != nil {
        return err
      }
      if _, err := fmt.Fprintf(writer, formatHeader, "", "...", "...", ""); err != nil {
        return err
      }
      return printMetaRow(writer)
    }); err != nil {
    return err
  }
  return nil
}

func (granges GRanges) PrintPretty(n int) string {
  var buffer bytes.Buffer
  writer := bufio.NewWriter(&buffer)

  if err := granges.WritePretty(writer, n); err != nil {
    return ""
  }
  writer.Flush()

  return buffer.String()
}
