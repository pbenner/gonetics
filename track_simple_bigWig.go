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

func (track *SimpleTrack) ReadBigWig(reader io.ReadSeeker, name string, f BinSummaryStatistics, binSize, binOverlap int, init float64) error {

  bwr, err := NewBigWigReader(reader)
  if err != nil {
    return err
  }
  // extract all sequences
  sequences := [][]float64{}
  for _, seqname := range bwr.Genome.Seqnames {
    if s, b, err := bwr.QuerySequence(seqname, f, binSize, binOverlap, init); err != nil {
      return err
    } else {
      if binSize == 0 {
        binSize = b
      }
      sequences = append(sequences, s)
    }
  }

  // create new track
  if tmp, err := NewSimpleTrack(name, sequences, bwr.Genome, binSize); err != nil {
    return err
  } else {
    *track = tmp
  }
  return nil
}

func (track *SimpleTrack) ImportBigWig(filename string, name string, s BinSummaryStatistics, binSize, binOverlap int, init float64) error {
  f, err := OpenBigWigFile(filename)
  if err != nil {
    return err
  }
  defer f.Close()

  if err := track.ReadBigWig(f, name, s, binSize, binOverlap, init); err != nil {
    return fmt.Errorf("importing bigWig file from `%s' failed: %v", filename, err)
  }
  return nil
}

/* -------------------------------------------------------------------------- */

func (track SimpleTrack) WriteBigWig(writer io.WriteSeeker, args... interface{}) error {
  return GenericTrack{track}.WriteBigWig(writer, args...)
}

func (track SimpleTrack) ExportBigWig(filename string, args... interface{}) error {
  f, err := os.Create(filename)
  if err != nil {
    return err
  }
  defer f.Close()

  if err := track.WriteBigWig(f, args...); err != nil {
    return fmt.Errorf("exporting bigWig file to `%s' failed: %v", filename, err)
  }
  return nil
}
