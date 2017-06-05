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
import "os"
import "strconv"
import "strings"

/* -------------------------------------------------------------------------- */

// Export GRanges as a table. The first line contains the header
// of the table.
func (granges GRanges) WriteTable(w io.Writer, header, strand bool) error {
  // print header
  if header {
    if strand {
      if _, err := fmt.Fprintf(w, "%14s %10s %10s %6s", "seqnames", "from", "to", "strand"); err != nil {
        return err
      }
    } else {
      if _, err := fmt.Fprintf(w, "%14s %10s %10s", "seqnames", "from", "to"); err != nil {
        return err
      }
    }
    granges.Meta.WriteTableRow(w, -1)
    if _, err := fmt.Fprintf(w, "\n"); err != nil {
      return err
    }
  }
  // print data
  for i := 0; i < granges.Length(); i++ {
    if _, err := fmt.Fprintf(w,  "%14s", granges.Seqnames[i]); err != nil {
      return err
    }
    if _, err := fmt.Fprintf(w, " %10d", granges.Ranges[i].From); err != nil {
      return err
    }
    if _, err := fmt.Fprintf(w, " %10d", granges.Ranges[i].To); err != nil {
      return err
    }
    if strand {
      if len(granges.Strand) > 0 {
        if _, err := fmt.Fprintf(w, " %6c", granges.Strand[i]); err != nil {
          return err
        }
      } else {
        if _, err := fmt.Fprintf(w, " %6c", '*'); err != nil {
          return err
        }
      }
    }
    granges.Meta.WriteTableRow(w, i)
    if _, err := fmt.Fprintf(w, "\n"); err != nil {
      return err
    }
  }
  return nil
}

func (granges GRanges) ExportTable(filename string, header, strand, compress bool) error {
  var buffer bytes.Buffer

  w := bufio.NewWriter(&buffer)
  if err := granges.WriteTable(w, header, strand); err != nil {
    return err
  }
  w.Flush()

  return writeFile(filename, &buffer, compress)
}

func (granges *GRanges) ReadTable(s io.ReadSeeker, names, types []string) error {
  var r io.Reader

  // check if file is compressed
  if g, err := gzip.NewReader(s); err != nil {
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
      }
    }
  }
  if colSeqname == -1 {
    return fmt.Errorf("is missing a seqname column")
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
  return granges.ReadTable(f, names, types)
}
