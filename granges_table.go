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
import "os"
import "strconv"
import "strings"

/* -------------------------------------------------------------------------- */

// Export GRanges as a table. The first line contains the header
// of the table.
func (granges GRanges) WriteTable(filename string, header, strand, compress bool) error {
  var buffer bytes.Buffer

  w := bufio.NewWriter(&buffer)

  // print header
  if header {
    if strand {
      fmt.Fprintf(w, "%14s %10s %10s %6s", "seqnames", "from", "to", "strand")
    } else {
      fmt.Fprintf(w, "%14s %10s %10s", "seqnames", "from", "to")
    }
    granges.Meta.WriteTableRow(w, -1)
    fmt.Fprintf(w, "\n")
  }
  // print data
  for i := 0; i < granges.Length(); i++ {
    fmt.Fprintf(w,  "%14s", granges.Seqnames[i])
    fmt.Fprintf(w, " %10d", granges.Ranges[i].From)
    fmt.Fprintf(w, " %10d", granges.Ranges[i].To)
    if strand {
      if len(granges.Strand) > 0 {
        fmt.Fprintf(w, " %6c", granges.Strand[i])
      } else {
        fmt.Fprintf(w, " %6c", '*')
      }
    }
    granges.Meta.WriteTableRow(w, i)
    fmt.Fprintf(w, "\n")
  }
  w.Flush()
  return writeFile(filename, &buffer, compress)
}

func (granges *GRanges) ReadTable(filename string, names, types []string) error {
  colSeqname := -1
  colFrom    := -1
  colTo      := -1
  colStrand  := -1

  var scanner *bufio.Scanner
  // open file
  f, err := os.Open(filename)
  if err != nil {
    return err
  }
  defer f.Close()

  // check if file is gzipped
  if isGzip(filename) {
    g, err := gzip.NewReader(f)
    if err != nil {
      return err
    }
    defer g.Close()
    scanner = bufio.NewScanner(g)
  } else {
    scanner = bufio.NewScanner(f)
  }

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
    return fmt.Errorf("ReadTable(): table `%s' is missing a seqname column", filename)
  }
  if colFrom == -1 {
    return fmt.Errorf("ReadTable(): table `%s' is missing a from column", filename)
  }
  if colTo == -1 {
    return fmt.Errorf("ReadTable(): table `%s' is missing a to column", filename)
  }
  // scan data
  for i:= 2; scanner.Scan(); i++ {
    fields := strings.Fields(scanner.Text())
    if len(fields) == 0 {
      continue
    }
    // parse seqname
    if len(fields) < colSeqname {
      return fmt.Errorf("ReadTable(): invalid table `%s'", filename)
    }
    // parse from
    if len(fields) < colFrom {
      return fmt.Errorf("ReadTable(): invalid table `%s'", filename)
    }
    v1, err := strconv.ParseInt(fields[colFrom], 10, 64)
    if err != nil {
      return fmt.Errorf("parsing `from' column `%d' failed at line `%d': %v", colFrom+1, i, err)
    }
    // parse to
    if len(fields) < colTo {
      return fmt.Errorf("ReadTable(): invalid table `%s'", filename)
    }
    v2, err := strconv.ParseInt(fields[colTo], 10, 64)
    if err != nil {
      return fmt.Errorf("parsing `to' column `%d' failed at line `%d': %v", colTo+1, i, err)
    }
    granges.Seqnames = append(granges.Seqnames, fields[colSeqname])
    granges.Ranges   = append(granges.Ranges,   NewRange(int(v1), int(v2)))
    if colStrand != -1 {
      if len(fields) < colStrand {
        return fmt.Errorf("ReadTable(): invalid table `%s'", filename)
      }
      granges.Strand = append(granges.Strand,   fields[colStrand][0])
    } else {
      granges.Strand = append(granges.Strand,   '*')
    }
  }
  return granges.Meta.ReadTable(filename, names, types)
}
