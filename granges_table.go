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
import "compress/gzip"
import "io"
import "io/ioutil"
import "os"
import "strconv"
import "strings"

/* -------------------------------------------------------------------------- */

// Export GRanges as a table. The first line contains the header
// of the table.
func (granges GRanges) WriteTable(writer io.Writer, header, strand bool, args ...interface{}) error {
  // pretty print meta data and create a scanner reading
  // the resulting string
  metaStr     := granges.Meta.PrintTable(header, args...)
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
    if err := updateMaxWidth("%s", widths, 0, granges.Seqnames[i]); err != nil {
      return err
    }
    if err := updateMaxWidth("%d", widths, 1, granges.Ranges[i].From); err != nil {
      return err
    }
    if err := updateMaxWidth("%d", widths, 2, granges.Ranges[i].To); err != nil {
      return err
    }
    if err := updateMaxWidth("%c", widths, 3, granges.Strand[i]); err != nil {
      return err
    }
    return nil
  }
  printMetaRow := func(writer io.Writer) error {
    if granges.MetaLength() != 0 {
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
    if strand {
      if _, err := fmt.Fprintf(writer, format, "seqnames", "from", "to", "strand"); err != nil {
        return err
      }
    } else {
      if _, err := fmt.Fprintf(writer, format, "seqnames", "from", "to"); err != nil {
        return err
      }
    }
    return printMetaRow(writer)
  }
  printRow := func(writer io.Writer, format string, i int) error {
    if i != 0 {
      if _, err := fmt.Fprintf(writer, "\n"); err != nil {
        return err
      }
    }
    if strand {
      if _, err := fmt.Fprintf(writer, format,
        granges.Seqnames[i],
        granges.Ranges[i].From,
        granges.Ranges[i].To,
        granges.Strand[i]); err != nil {
        return err
        }
    } else {
      if _, err := fmt.Fprintf(writer, format,
        granges.Seqnames[i],
        granges.Ranges[i].From,
        granges.Ranges[i].To); err != nil {
        return err
      }
    }
    return printMetaRow(writer)
  }
  applyRows := func(f1 func(i int) error) error {
    // apply to all entries
    for i := 0; i < granges.Length(); i++ {
      if err := f1(i); err != nil {
        return err
      }
    }
    return nil
  }
  // maximum column widths
  widths := []int{8, 4, 2, 6}
  // determine column widths
  if err := applyRows(func(i int) error { return updateMaxWidths(i, widths) }); err != nil {
    return err
  }
  // generate format strings
  var formatRow, formatHeader string
  if strand {
    formatRow    = fmt.Sprintf("%%%ds %%%dd %%%dd %%%dc", widths[0], widths[1], widths[2], widths[3])
    formatHeader = fmt.Sprintf("%%%ds %%%ds %%%ds %%%ds", widths[0], widths[1], widths[2], widths[3])
  } else {
    formatRow    = fmt.Sprintf("%%%ds %%%dd %%%dd", widths[0], widths[1], widths[2])
    formatHeader = fmt.Sprintf("%%%ds %%%ds %%%ds", widths[0], widths[1], widths[2])
  }
  // pring header
  if header {
    if err := printHeader(writer, formatHeader); err != nil {
      return err
    }
    if _, err := fmt.Fprintf(writer, "\n"); err != nil {
      return err
    }
  }
  // print rows
  if err := applyRows(
    func(i int) error {
      return printRow(writer, formatRow, i)
    }); err != nil {
    return err
  }
  return nil
}

/* -------------------------------------------------------------------------- */

func (granges GRanges) PrintTable(header, strand bool, args ...interface{}) string {
  var buffer bytes.Buffer
  writer := bufio.NewWriter(&buffer)

  if err := granges.WriteTable(writer, header, strand, args...); err != nil {
    return ""
  }
  writer.Flush()

  return buffer.String()
}

/* -------------------------------------------------------------------------- */

func (granges GRanges) ExportTable(filename string, header, strand, compress bool, args ...interface{}) error {
  var buffer bytes.Buffer

  w := bufio.NewWriter(&buffer)
  if err := granges.WriteTable(w, header, strand, args...); err != nil {
    return err
  }
  fmt.Fprintf(w, "\n")
  w.Flush()

  return writeFile(filename, &buffer, compress)
}

/* -------------------------------------------------------------------------- */

func (granges *GRanges) ReadTable(s io.ReadSeeker, names, types []string) error {
  var r io.Reader

  // check if file is compressed
  if g, err := gzip.NewReader(s); err != nil {
    if _, err := s.Seek(0, io.SeekStart); err != nil {
      return err
    }
    r = s
  } else {
    r = g
    defer g.Close()
  }

  scanner := bufio.NewScanner(r)

  colSeqname := -1
  colFrom    := -1
  colTo      := -1
  colStrand  := -1

  // scan header
  if scanner.Scan() {
    fields := strings.Fields(scanner.Text())
    for i := 0; i < len(fields); i++ {
      switch fields[i] {
      case "seqnames":
        colSeqname = i
      case "from":
        colFrom = i
      case "to":
        colTo = i
      case "strand":
        colStrand = i
      case "start":
        // alternative name for `from' column
        if colFrom == -1 {
          colFrom = i
        }
      case "end":
        // alternative name for `to' column
        if colTo == -1 {
          colTo = i
        }
      }
    }
  } else {
    return fmt.Errorf("reading from file failed")
  }
  if colSeqname == -1 {
    return fmt.Errorf("is missing a seqnames column")
  }
  if colFrom == -1 {
    return fmt.Errorf("is missing a from column")
  }
  if colTo == -1 {
    return fmt.Errorf("is missing a to column")
  }
  // scan data
  for i:= 2; scanner.Scan(); i++ {
    fields := strings.Fields(scanner.Text())
    if len(fields) == 0 {
      continue
    }
    // parse seqname
    if len(fields) < colSeqname {
      return fmt.Errorf("invalid table")
    }
    // parse from
    if len(fields) < colFrom {
      return fmt.Errorf("invalid table")
    }
    v1, err := strconv.ParseInt(fields[colFrom], 10, 64)
    if err != nil {
      return fmt.Errorf("parsing `from' column `%d' failed at line `%d': %v", colFrom+1, i, err)
    }
    // parse to
    if len(fields) < colTo {
      return fmt.Errorf("invalid table")
    }
    v2, err := strconv.ParseInt(fields[colTo], 10, 64)
    if err != nil {
      return fmt.Errorf("parsing `to' column `%d' failed at line `%d': %v", colTo+1, i, err)
    }
    granges.Seqnames = append(granges.Seqnames, fields[colSeqname])
    granges.Ranges   = append(granges.Ranges,   NewRange(int(v1), int(v2)))
    if colStrand != -1 {
      if len(fields) < colStrand {
        return fmt.Errorf("invalid table")
      }
      granges.Strand = append(granges.Strand,   fields[colStrand][0])
    } else {
      granges.Strand = append(granges.Strand,   '*')
    }
  }
  // rewind file
  s.Seek(0, io.SeekStart)
  // check if file is compressed
  if g, err := gzip.NewReader(s); err != nil {
    if _, err := s.Seek(0, io.SeekStart); err != nil {
      return err
    }
    r = s
  } else {
    r = g
    defer g.Close()
  }
  return granges.Meta.ReadTable(r, names, types)
}

func (granges *GRanges) ImportTable(filename string, names, types []string) error {
  // open file
  f, err := os.Open(filename)
  if err != nil {
    return err
  }
  defer f.Close()
  if err := granges.ReadTable(f, names, types); err != nil {
    return fmt.Errorf("`%s' %s", filename, err)
  } else {
    return nil
  }
}
