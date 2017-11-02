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
import "io"
import "os"

/* -------------------------------------------------------------------------- */

func (track GenericTrack) writeBigWig_reductionLevels(parameters BigWigParameters) []int {
  c := BbiResIncrement*track.GetBinSize()
  // reduction levels
  n := []int{}
  // length of the longest track
  l := 0
  // get length of longest track
  for _, length := range track.GetGenome().Lengths {
    if length/track.GetBinSize() > l {
      l = length/track.GetBinSize()
    }
  }
  // initial zoom level
  r := iMax(100, c)
  // compute number of zoom levels
  for len(n) <= BbiMaxZoomLevels {
    if l/r > parameters.ItemsPerSlot {
      n = append(n, r)
      r = r*c
    } else {
      break
    }
  }
  return n
}

func (track GenericTrack) WriteBigWig(writer io.WriteSeeker, args... interface{}) error {

  parameters := DefaultBigWigParameters()

  // parse arguments
  for i := 0; i < len(args); i++ {
    switch v := args[i].(type) {
    case BigWigParameters:
      parameters = v
    default:
      return fmt.Errorf("WriteBigWig(): invalid arguments")
    }
  }
  // get reduction levels for zoomed data
  if parameters.ReductionLevels == nil {
    parameters.ReductionLevels = track.writeBigWig_reductionLevels(parameters)
  }
  // create new bigWig writer
  bww, err := NewBigWigWriter(writer, track.GetGenome(), parameters)
  if err != nil {
    return err
  }
  // write data
  for _, name := range track.GetSeqNames() {
    sequence, err := track.GetSequence(name); if err != nil {
      return err
    }
    if err := bww.Write(name, sequence.sequence, track.GetBinSize()); err != nil {
      return err
    }
  }
  if err := bww.WriteIndex(); err != nil {
    return err
  }
  // write zoomed data
  for i, reductionLevel := range parameters.ReductionLevels {
    // save current offset as the beginning of zoomed data for reduction
    // level i
    if err := bww.StartZoomData(i); err != nil {
      return err
    }
    for _, name := range track.GetSeqNames() {
      sequence, err := track.GetSequence(name); if err != nil {
        return err
      }
      if err := bww.WriteZoom(name, sequence.sequence, track.GetBinSize(), reductionLevel, i); err != nil {
        return err
      }
    }
    // write index for this reduction level
    if err := bww.WriteIndexZoom(i); err != nil {
      return err
    }
  }
  bww.Close()

  return nil
}

func (track GenericTrack) ExportBigWig(filename string, args... interface{}) error {
  f, err := os.Create(filename)
  if err != nil {
    return err
  }
  defer f.Close()

  return track.WriteBigWig(f, args...)
}
