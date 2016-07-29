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
import "errors"
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

func (g *GRanges) ReadTable(filename string, names, types []string) error {
  result    := GRanges{}
  hasStrand := false

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
    // get number of columns
    if len(fields) < 4 {
      return errors.New("invalid table")
    }
    if fields[0] != "seqnames" || fields[1] != "from" || fields[2] != "to" {
      return errors.New("invalid table")
    }
    if fields[3] == "strand" {
      hasStrand = true
    }
  }
  for scanner.Scan() {
    fields := strings.Fields(scanner.Text())
    if len(fields) == 0 {
      continue
    }
    if len(fields) < 4 {
      return errors.New("invalid table")
    }
    v1, err := strconv.ParseInt(fields[1], 10, 64)
    if err != nil {
      return err
    }
    v2, err := strconv.ParseInt(fields[2], 10, 64)
    if err != nil {
      return err
    }
    result.Seqnames = append(result.Seqnames, fields[0])
    result.Ranges   = append(result.Ranges,   NewRange(int(v1), int(v2)))
    if hasStrand {
      result.Strand = append(result.Strand,   fields[3][0])
    } else {
      result.Strand = append(result.Strand,   '*')
    }
  }
  return result.Meta.ReadTable(filename, names, types)
}
