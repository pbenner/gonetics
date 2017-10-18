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

//import "fmt"

/* -------------------------------------------------------------------------- */

func (track *SimpleTrack) ReadBigWig(filename, name string, f BinSummaryStatistics, binSize, binOverlap int, init float64) error {

  bwr, err := NewBigWigReader(filename)
  if err != nil {
    return err
  }
  // extract all sequences
  seqnames  := bwr.Genome.Seqnames
  sequences := [][]float64{}
  for i := 0; i < bwr.Genome.Length(); i++ {
    seqname := bwr.Genome.Seqnames[i]
    if s, b, err := bwr.QuerySequence(seqname, f, binSize, binOverlap, init); err != nil {
      return err
    } else {
      if binSize == 0 {
        binSize = b
      }
      sequences = append(sequences, s)
    }
  }
  bwr.Close()

  // create new track
  if tmp, err := NewSimpleTrack(name, seqnames, sequences, binSize); err != nil {
    return err
  } else {
    *track = tmp
  }
  return nil
}

/* -------------------------------------------------------------------------- */

func (track SimpleTrack) WriteBigWig(filename string, genome Genome, args... interface{}) error {
  return GenericTrack{track}.WriteBigWig(filename, genome, args...)
}
