/* Copyright (C) 2017 Philipp Benner
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

func (granges GRanges) WriteBedGraph(w io.Writer) error {
  values := granges.GetMetaFloat("values")

  if len(values) == 0 {
    return fmt.Errorf("values column required for bedGraph export")
  }

  for i := 0; i < granges.Length(); i++ {
    if _, err := fmt.Fprintf(w,   "%s", granges.Seqnames[i]); err != nil {
      return err
    }
    if _, err := fmt.Fprintf(w, "\t%d", granges.Ranges[i].From); err != nil {
      return err
    }
    if _, err := fmt.Fprintf(w, "\t%d", granges.Ranges[i].To); err != nil {
      return err
    }
    if _, err := fmt.Fprintf(w, "\t%f", values[i]); err != nil {
      return err
    }
  }
  return nil
}

func (granges GRanges) ExportBedGraph(filename string, compress bool) error {
  var buffer bytes.Buffer

  w := bufio.NewWriter(&buffer)
  if err := granges.WriteBedGraph(w); err != nil {
    return err
  }
  w.Flush()

  return writeFile(filename, &buffer, compress)
}

func (g *GRanges) ReadBedGraph(r io.Reader) error {
  scanner := bufio.NewScanner(r)
  // it seems that buffering the data does not increase
  // performance
  seqnames := []string{}
  from     := []int{}
  to       := []int{}
  values   := []float64{}

  for scanner.Scan() {
    fields := strings.Fields(scanner.Text())
    if len(fields) == 0 {
      continue
    }
    if len(fields) != 4 {
      return fmt.Errorf("ReadBedGraph(): BedGraph file must have four columns!")
    }
    t1, err := strconv.ParseInt(fields[1], 10, 64); if err != nil {
      return err
    }
    t2, err := strconv.ParseInt(fields[2], 10, 64); if err != nil {
      return err
    }
    t3, err := strconv.ParseFloat(fields[3], 64); if err != nil {
      return err
    }
    seqnames = append(seqnames, fields[0])
    from     = append(from,     int(t1))
    to       = append(to,       int(t2))
    values   = append(values,   t3)
  }
  *g = NewGRanges(seqnames, from, to, nil)
  g.AddMeta("values", values)

  return nil
}

func (g *GRanges) ImportBedGraph(filename string) error {
  var r io.Reader
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
    r = g
  } else {
    r = f
  }
  return g.ReadBedGraph(r)
}
