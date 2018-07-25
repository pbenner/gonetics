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

import   "fmt"
import   "math"
import   "bufio"
import   "bytes"
import   "compress/gzip"
import   "io"
import   "io/ioutil"
import   "os"

/* -------------------------------------------------------------------------- */

func getNColors(n int) []string {
  // number of levels for each color
  k := int(math.Cbrt(float64(n)))+1
  // delta for each color
  d := 255/k
  s := []string{}
  for i := 0; i < n; i++ {
    r := (i/k/k % k) * d
    g := (i/k   % k) * d
    b := (i     % k) * d
    s  = append(s, fmt.Sprintf("%d,%d,%d", r,g,b))
  }
  return s
}

/* -------------------------------------------------------------------------- */

func (track GenericTrack) exportSegmentation(granges GRanges, bedFilename, name, description string, compress bool) error {
  buffer := new(bytes.Buffer)

  w := bufio.NewWriter(buffer)
  if _, err := fmt.Fprintf(w, "track name=\"%s\" description=\"%s\" visibility=1 itemRgb=\"On\"\n", name, description); err != nil {
    return err
  }
  if err := granges.WriteBed9(w); err != nil {
    return err
  }
  w.Flush()

  if compress {
    b := new(bytes.Buffer)
    w := gzip.NewWriter(b)
    io.Copy(w, buffer)
    w.Close()
    buffer = b
  }
  return ioutil.WriteFile(bedFilename, buffer.Bytes(), 0666)
}

func (track GenericTrack) ExportSegmentation(bedFilename, bedName, bedDescription string, compress bool, stateNames, rgbChart []string) error {
  r, err := track.GRanges("state"); if err != nil {
    return err
  }
  state      := r.GetMetaFloat("state")
  name       := make([]string, len(state))
  score      := make([]int,    len(state))
  thickStart := make([]int,    len(state))
  thickEnd   := make([]int,    len(state))
  itemRgb    := make([]string, len(state))
  rgbMap     := make(map[string]string)

  if len(stateNames) == 0 {
    // determine number of states
    sMax := 0
    for i := 0; i < r.Length(); i++ {
      if s := int(state[i]); s > sMax {
        sMax = s
      }
    }
    // generate state names
    stateNames = make([]string, sMax+1)
    for i := 0; i < len(stateNames); i++ {
      stateNames[i] = fmt.Sprintf("s%d", i)
    }
  }
  { // determine number of distinct state names
    for i := 0; i < len(stateNames); i++ {
      rgbMap[stateNames[i]] = ""
    }
  }
  if len(rgbChart) == 0 {
    // get a color for each distinct state name
    rgbChart = getNColors(len(rgbMap))
  }
  if len(rgbMap) > len(rgbChart) {
    return fmt.Errorf("insufficient number of RGB colors")
  }
  { // fill rgbMap (stateName -> rgb color)
    k := 0
    for state, _ := range rgbMap {
      rgbMap[state] = rgbChart[k]; k++
    }
  }

  for i := 0; i < r.Length(); i++ {
    s := int(state[i])
    if s < 0 || math.Floor(state[i]) != state[i] {
      return fmt.Errorf("invalid state `%f' at `%s:%d-%d", state[i], r.Seqnames[i], r.Ranges[i].From, r.Ranges[i].To)
    }
    if s >= len(stateNames) {
      return fmt.Errorf("insufficient number of state names")
    }
    name      [i] = stateNames[s]
    score     [i] = 0
    thickStart[i] = r.Ranges[i].From
    thickEnd  [i] = r.Ranges[i].To
    itemRgb   [i] = rgbMap[stateNames[s]]
  }
  r.AddMeta("name",       name)
  r.AddMeta("thickStart", thickStart)
  r.AddMeta("thickEnd",   thickEnd)
  r.AddMeta("itemRgb",    itemRgb)
  // write result to file
  return track.exportSegmentation(r, bedFilename, bedName, bedDescription, compress)
}

/* -------------------------------------------------------------------------- */

func importSegmentation(filename string) (GRanges, error) {
  var r io.Reader
  var g GRanges
  // open file
  f, err := os.Open(filename)
  if err != nil {
    return g, err
  }
  defer f.Close()
  // check if file is gzipped
  if isGzip(filename) {
    gz, err := gzip.NewReader(f)
    if err != nil {
      return g, err
    }
    defer gz.Close()
    r = gz
  } else {
    r = f
  }
  // skip track line
  reader := bufio.NewReader(r)
  reader.ReadLine()
  // parse remaining file
  if err := g.ReadBed9(reader); err != nil {
    return g, err
  } else {
    return g, nil
  }
}

func (track GenericMutableTrack) ImportSegmentation(bedFilename string) ([]string, error) {
  var s TrackMutableSequence
  if r, err := importSegmentation(bedFilename); err != nil {
    return nil, err
  } else {
    binSize  := track.GetBinSize()
    states   := r.GetMetaStr("name")
    stateMap := make(map[string]int)

    if len(states) != r.Length() {
      return nil, fmt.Errorf("invalid segmentation bed file: name column is missing")
    }
    for i := 0; i < r.Length(); i++ {
      seqname := r.Seqnames[i]
      from    := r.Ranges  [i].From
      to      := r.Ranges  [i].To
      if i == 0 || r.Seqnames[i-1] != r.Seqnames[i] {
        if s_, err := track.GetMutableSequence(seqname); err != nil {
          return nil, err
        } else {
          s = s_
        }
      }
      stateIdx := -1
      // generate new state idx if this state name was observed
      // for the first time
      if k, ok := stateMap[states[i]]; ok {
        stateIdx = k
      } else {
        stateIdx = len(stateMap)
        stateMap[states[i]] = stateIdx
      }

      for k := from; k < to; k += binSize {
        s.Set(k, float64(stateIdx))
      }
    }
    stateNames := make([]string, len(stateMap))
    for name, i := range stateMap {
      stateNames[i] = name
    }
    return stateNames, nil
  }
}
