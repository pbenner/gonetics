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
      fmt.Fprintf(w, "\t%d", 0)
    }
    if len(granges.Strand) > 0 && granges.Strand[i] != '*' {
      fmt.Fprintf(w, "\t%c", granges.Strand[i])
    } else {
      fmt.Fprintf(w, "\t%s", ".")
    }
    fmt.Fprintf(w, "\n")
  }
  w.Flush()

  return writeFile(filename, &buffer, compress)
}

func (granges GRanges) WriteBed9(filename string, compress bool) error {
  var buffer bytes.Buffer

  w := bufio.NewWriter(&buffer)

  name       := granges.GetMetaStr  ("name")
  score      := granges.GetMetaFloat("score")
  thickStart := granges.GetMetaInt  ("thickStart")
  thickEnd   := granges.GetMetaInt  ("thickEnd")
  itemRgb    := granges.GetMetaStr  ("itemRgb")

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
      fmt.Fprintf(w, "\t%d", 0)
    }
    if len(granges.Strand) > 0 && granges.Strand[i] != '*' {
      fmt.Fprintf(w, "\t%c", granges.Strand[i])
    } else {
      fmt.Fprintf(w, "\t%s", ".")
    }
    if len(thickStart) > 0 {
      fmt.Fprintf(w, "\t%d", thickStart[i])
    } else {
      fmt.Fprintf(w, "\t%d", granges.Ranges[i].From)
    }
    if len(thickEnd) > 0 {
      fmt.Fprintf(w, "\t%d", thickEnd[i])
    } else {
      fmt.Fprintf(w, "\t%d", granges.Ranges[i].To)
    }
    if len(itemRgb) > 0 {
      fmt.Fprintf(w, "\t%s", itemRgb[i])
    } else {
      fmt.Fprintf(w, "\t%s", "0,0,0")
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
    t1, err := strconv.ParseInt(fields[1], 10, 64); if err != nil {
      return err
    }
    t2, err := strconv.ParseInt(fields[2], 10, 64); if err != nil {
      return err
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
      return fmt.Errorf("ReadBed6(): Bed file must have at least six columns!")
    }
    t1, err := strconv.ParseInt(fields[1], 10, 64); if err != nil {
      return err
    }
    t2, err := strconv.ParseInt(fields[2], 10, 64); if err != nil {
      return err
    }
    t3, err := strconv.ParseInt(fields[4], 10, 64); if err != nil {
      return err
    }
    seqnames = append(seqnames, fields[0])
    from     = append(from,     int(t1))
    to       = append(to,       int(t2))
    name     = append(name,     fields[3])
    score    = append(score,    int(t3))
    if fields[5][0] == '.' {
      strand   = append(strand, '*')
    } else {
      strand   = append(strand, fields[5][0])
    }
  }
  *g = NewGRanges(seqnames, from, to, strand)
  g.AddMeta("name",  name)
  g.AddMeta("score", score)

  return nil
}

func (g *GRanges) ReadBed9(filename string) error {
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
  seqnames   := []string{}
  from       := []int{}
  to         := []int{}
  name       := []string{}
  score      := []int{}
  strand     := []byte{}
  thickStart := []int{}
  thickEnd   := []int{}
  itemRgb    := []string{}

  for scanner.Scan() {
    fields := strings.Fields(scanner.Text())
    if len(fields) == 0 {
      continue
    }
    if len(fields) < 6 {
      return fmt.Errorf("ReadBed6(): Bed file must have at least six columns!")
    }
    t1, err := strconv.ParseInt(fields[1], 10, 64); if err != nil {
      return err
    }
    t2, err := strconv.ParseInt(fields[2], 10, 64); if err != nil {
      return err
    }
    t3, err := strconv.ParseInt(fields[4], 10, 64); if err != nil {
      return err
    }
    t4, err := strconv.ParseInt(fields[6], 10, 64); if err != nil {
      return err
    }
    t5, err := strconv.ParseInt(fields[7], 10, 64); if err != nil {
      return err
    }
    seqnames   = append(seqnames,   fields[0])
    from       = append(from,       int(t1))
    to         = append(to,         int(t2))
    name       = append(name,       fields[3])
    score      = append(score,      int(t3))
    thickStart = append(thickStart, int(t4))
    thickEnd   = append(thickEnd,   int(t5))
    itemRgb    = append(itemRgb,    fields[8])
    if fields[5][0] == '.' {
      strand   = append(strand, '*')
    } else {
      strand   = append(strand, fields[5][0])
    }
  }
  *g = NewGRanges(seqnames, from, to, strand)
  g.AddMeta("name",       name)
  g.AddMeta("score",      score)
  g.AddMeta("thickStart", thickStart)
  g.AddMeta("thickEnd",   thickEnd)
  g.AddMeta("itemRgb",    itemRgb)

  return nil
}
