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

func (genes Genes) WriteTable(writer io.Writer, header bool, args ...interface{}) error {
  // pretty print meta data and create a scanner reading
  // the resulting string
  meta        := genes.Meta.Clone()
  meta.DeleteMeta("names")
  meta.DeleteMeta("cds")
  metaStr     := meta.PrintTable(header, args...)
  metaReader  := strings.NewReader(metaStr)
  metaScanner := bufio.NewScanner(metaReader)

  // compute the width of a single cell
  updateMaxWidth := func(format string, widths []int, j int, args ...interface{}) error {
    if width, err := fmt.Fprintf(ioutil.Discard, format, args...); err != nil {
      return err
    } else {
      if width > widths[j] {
        widths[j] = width
      }
    }
    return nil
  }
  // compute widths of all cells in row i
  updateMaxWidths := func(i int, widths []int) error {
    if err := updateMaxWidth("%s", widths, 0, genes.Names[i]); err != nil {
      return err
    }
    if err := updateMaxWidth("%s", widths, 1, genes.Seqnames[i]); err != nil {
      return err
    }
    if err := updateMaxWidth("%c", widths, 2, genes.Strand[i]); err != nil {
      return err
    }
    if err := updateMaxWidth("%d", widths, 3, genes.Ranges[i].From); err != nil {
      return err
    }
    if err := updateMaxWidth("%d", widths, 4, genes.Ranges[i].To); err != nil {
      return err
    }
    if err := updateMaxWidth("%d", widths, 5, genes.Cds[i].From); err != nil {
      return err
    }
    if err := updateMaxWidth("%d", widths, 6, genes.Cds[i].To); err != nil {
      return err
    }
    return nil
  }
  printMetaRow := func(writer io.Writer) error {
    if meta.MetaLength() != 0 {
      if _, err := fmt.Fprintf(writer, " "); err != nil {
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
    if _, err := fmt.Fprintf(writer, format,
      "names", "seqnames", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd"); err != nil {
      return err
    }
    if err := printMetaRow(writer); err != nil {
      return err
    }
    if _, err := fmt.Fprintf(writer, "\n"); err != nil {
      return err
    }
    return nil
  }
  printRow := func(writer io.Writer, format string, i int) error {
    if i != 0 {
      if _, err := fmt.Fprintf(writer, "\n"); err != nil {
        return err
      }
    }
    if _, err := fmt.Fprintf(writer, format,
      genes.Names[i],
      genes.Seqnames[i],
      genes.Strand[i],
      genes.Ranges[i].From,
      genes.Ranges[i].To,
      genes.Cds[i].From,
      genes.Cds[i].To); err != nil {
      return err
    }
    if err := printMetaRow(writer); err != nil {
      return err
    }
    return nil
  }
  applyRows := func(f1 func(i int) error) error {
    // apply to all entries
    for i := 0; i < genes.Length(); i++ {
      if err := f1(i); err != nil {
        return err
      }
    }
    return nil
  }
  // maximum column widths
  widths := []int{5, 8, 6, 7, 5, 8, 6}
  // determine column widths
  if err := applyRows(func(i int) error { return updateMaxWidths(i, widths) }); err != nil {
    return err
  }
  // generate format strings
  formatRow    := fmt.Sprintf("%%%ds %%%ds %%%dc %%%dd %%%dd %%%dd %%%dd",
    widths[0], widths[1], widths[2], widths[3], widths[4], widths[5], widths[6])
  formatHeader := fmt.Sprintf("%%%ds %%%ds %%%ds %%%ds %%%ds %%%ds %%%ds",
    widths[0], widths[1], widths[2], widths[3], widths[4], widths[5], widths[6])
  // pring header
  if header {
    if err := printHeader(writer, formatHeader); err != nil {
      return err
    }
  }
  // print rows
  if err := applyRows(
    func(i int) error {
      if err := printRow(writer, formatRow, i); err != nil {
        return err
      }
      return nil
    }); err != nil {
    return err
  }
  return nil
}

/* -------------------------------------------------------------------------- */

func (genes Genes) PrintTable(header, strand bool, args ...interface{}) string {
  var buffer bytes.Buffer
  writer := bufio.NewWriter(&buffer)

  if err := genes.WriteTable(writer, header, args...); err != nil {
    return ""
  }
  writer.Flush()

  return buffer.String()
}

/* -------------------------------------------------------------------------- */

func (genes Genes) ExportTable(filename string, header, strand, compress bool, args ...interface{}) error {
  var buffer bytes.Buffer

  w := bufio.NewWriter(&buffer)
  if err := genes.WriteTable(w, header, args...); err != nil {
    return err
  }
  w.Flush()

  return writeFile(filename, &buffer, compress)
}
