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
import "compress/gzip"
import "fmt"
import "io"
import "os"
import "strconv"
import "strings"

/* -------------------------------------------------------------------------- */

// Import data from wiggle files.
func (track *SimpleTrack) ReadBedGraph(reader io.Reader) error {
  fields  := []string{}
  binsize := track.BinSize
  scanner := bufio.NewScanner(reader)

  if !scanner.Scan() {
    return nil
  }
  // current sequence and name
  cur_seq  := []float64{}
  cur_name := ""
  for scanner.Scan() {
    fields = strings.Fields(scanner.Text())
    if len(fields) == 0 {
      break
    }
    if len(fields) != 4 {
      return fmt.Errorf("ReadBedGraph(): bedGraph file must have four columns!")
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
    if name := fields[0]; name != cur_name {
      cur_seq  = track.Data[name]
      cur_name = name
    }
    from  := int(t1)
    to    := int(t2)
    value := float64(t3)
    if from < 0 || to < 0 {
      continue
    }
    if from/binsize >= len(cur_seq) || to/binsize >= len(cur_seq) {
      continue
    }
    for i := from/binsize; i < to/binsize; i++ {
      cur_seq[i] = value
    }
  }
  return nil
}

func (track *SimpleTrack) ImportBedGrah(filename string) error {
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
  return track.ReadBedGraph(r)
}
