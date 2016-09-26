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

// Export GRanges object as bed file with three columns.
func (granges GRanges) WriteBed3(filename string, compress bool) error {
  var buffer bytes.Buffer

  w := bufio.NewWriter(&buffer)

  // print data
  for i := 0; i < granges.Length(); i++ {
    fmt.Fprintf(w,   "%s", granges.Seqnames[i])
    fmt.Fprintf(w, "\t%d", granges.Ranges[i].From)
    fmt.Fprintf(w, "\t%d", granges.Ranges[i].To)
    fmt.Fprintf(w, "\n")
  }
  w.Flush()

  return writeFile(filename, &buffer, compress)
}

func (granges GRanges) WriteBed6(filename string, compress bool) error {
  var buffer bytes.Buffer

  w := bufio.NewWriter(&buffer)

  name  := granges.GetMetaStr  ("name")
  score := granges.GetMetaFloat("score")

  for i := 0; i < granges.Length(); i++ {
    fmt.Fprintf(w,   "%s", granges.Seqnames[i])
    fmt.Fprintf(w, "\t%d", granges.Ranges[i].From)
    fmt.Fprintf(w, "\t%d", granges.Ranges[i].To)
    if len(name) > 0 {
      fmt.Fprintf(w, "\t%s", name[i])
    } else {
      fmt.Fprintf(w, "\t%s", ".")
    }
    if len(score) > 0 {
      fmt.Fprintf(w, "\t%f", score[i])
    } else {
      fmt.Fprintf(w, "\t%s", ".")
    }
    if len(granges.Strand) > 0 {
      fmt.Fprintf(w, "\t%c", granges.Strand[i])
    } else {
      fmt.Fprintf(w, "\t%c", '*')
    }
    fmt.Fprintf(w, "\n")
  }
  w.Flush()

  return writeFile(filename, &buffer, compress)
}

// Import GRanges from a Bed file with 3 columns.
func (g *GRanges) ReadBed3(filename string) error {
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
  // it seems that buffering the data does not increase
  // performance
  seqnames := []string{}
  from     := []int{}
  to       := []int{}

  for scanner.Scan() {
    fields := strings.Fields(scanner.Text())
    if len(fields) == 0 {
      continue
    }
    if len(fields) < 3 {
      return fmt.Errorf("ReadBed3(): Bed file must have at least three columns!")
    }
    t1, e1 := strconv.ParseInt(fields[1], 10, 64)
    t2, e2 := strconv.ParseInt(fields[2], 10, 64)
    if e1 != nil {
      return e1
    }
    if e2 != nil {
      return e2
    }
    seqnames = append(seqnames, fields[0])
    from     = append(from,     int(t1))
    to       = append(to,       int(t2))
  }
  *g = NewGRanges(seqnames, from, to, []byte{})

  return nil
}

// Import GRanges from a Bed file with 6 columns.
func (g *GRanges) ReadBed6(filename string) error {
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
      return nil
    }
    defer g.Close()
    scanner = bufio.NewScanner(g)
  } else {
    scanner = bufio.NewScanner(f)
  }
  // it seems that buffering the data does not increase
  // performance
  seqnames := []string{}
  from     := []int{}
  to       := []int{}
  name     := []string{}
  score    := []int{}
  strand   := []byte{}

  for scanner.Scan() {
    fields := strings.Fields(scanner.Text())
    if len(fields) == 0 {
      continue
    }
    if len(fields) < 6 {
      fmt.Errorf("ReadBed6(): Bed file must have at least six columns!")
    }
    t1, e1 := strconv.ParseInt(fields[1], 10, 64)
    if e1 != nil {
      return e1
    }
    t2, e2 := strconv.ParseInt(fields[2], 10, 64)
    if e2 != nil {
      return e2
    }
    t3, e3 := strconv.ParseInt(fields[4], 10, 64)
    if e3 != nil {
      return e3
    }
    seqnames = append(seqnames, fields[0])
    from     = append(from,     int(t1))
    to       = append(to,       int(t2))
    name     = append(name,     fields[3])
    score    = append(score,    int(t3))
    strand   = append(strand,   fields[5][0])
  }
  *g = NewGRanges(seqnames, from, to, strand)
  g.AddMeta("name",  name)
  g.AddMeta("score", score)

  return nil
}
