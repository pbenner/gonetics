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

/* -------------------------------------------------------------------------- */

func (track *Track) readBigWig_block(buffer []byte, genome Genome) error {
  decoder, err := NewBbiBlockDecoder(buffer)
  if err != nil {
    return err
  }
  if idx := int(decoder.Header.ChromId); idx < 0 || idx > genome.Length() {
    return fmt.Errorf("invalid chromosome id")
  }
  // convert chromosome id to sequence name
  seqname := genome.Seqnames[int(decoder.Header.ChromId)]

  // allocate track if this is the first buffer
  if len(track.Data) == 0 && track.Binsize == 0 {
    *track = AllocTrack("", genome, int(decoder.Header.Span))
  }
  switch decoder.Header.Type {
  case 2:
    if int(decoder.Header.Span) != track.Binsize {
      return fmt.Errorf("block has invalid span `%d' for track with bin size `%d'", decoder.Header.Span, track.Binsize)
    }
  case 3:
    if int(decoder.Header.Span) != track.Binsize {
      return fmt.Errorf("block has invalid span `%d' for track with bin size `%d'", decoder.Header.Span, track.Binsize)
    }
    if int(decoder.Header.Step) != track.Binsize {
      return fmt.Errorf("block has invalid step `%d' for track with bin size `%d'", decoder.Header.Span, track.Binsize)
    }
  }
  if seq, ok := track.Data[seqname]; !ok {
    return fmt.Errorf("sequence `%s' not vailable in track", seqname)
  } else {
    for t := range decoder.Read() {
      idx := track.Index(t.From)
      if idx >= len(seq) {
        return fmt.Errorf("position `%d' on sequence `%s' is out of range (trying to access bin `%d' but track has only `%d' bins)", t.From, seqname, idx, len(seq))
      } else {
        seq[track.Index(t.From)] = t.Value
      }
    }
  }
  return nil
}

func (track *Track) ReadBigWig(filename, name string) error {

  bwr, err := NewBigWigReader(filename)
  if err != nil {
    return err
  }
  for result := range bwr.ReadBlocks() {
    if result.Error != nil {
      return fmt.Errorf("reading `%s' failed: %v", filename, result.Error)
    }
    if err := track.readBigWig_block(result.Block, bwr.Genome); err != nil {
      return fmt.Errorf("reading `%s' failed: %v", filename, err)
    }
  }
  track.Name = name

  return nil
}

/* -------------------------------------------------------------------------- */

func (track *Track) writeBigWig_reductionLevels(parameters BigWigParameters) []int {
  // reduction levels
  n := []int{}
  // length of the longest track
  l := 0
  // get length of longest track
  for _, sequence := range track.Data {
    if len(sequence) > l {
      l = len(sequence)
    }
  }
  // compute number of zoom levels
  for i, r := 0, BbiResIncrement; i <= BbiMaxZoomLevels; i++ {
    if l > parameters.ItemsPerSlot {
      l = l/BbiResIncrement
      n = append(n, r)
      r = r*BbiResIncrement
    } else {
      break
    }
  }
  return n
}

func (track *Track) writeBigWig_zoomSequence(sequence []float64, reductionLevel int) []float64 {
  result := make([]float64, divIntUp(len(sequence), reductionLevel))
  for i := 0; i < len(sequence); i += reductionLevel {
    // compute mean value
    for j := 0; j < reductionLevel && i+j < len(sequence); j++ {
      result[i/reductionLevel] += sequence[i+j]
    }
    result[i/reductionLevel] /= float64(reductionLevel)
  }
  return result
}

func (track *Track) WriteBigWig(filename, description string, args... interface{}) error {

  parameters := DefaultBigWigParameters()

  // parse arguments
  for i := 0; i < len(args); i++ {
    switch v := args[i].(type) {
    case BigWigParameters: parameters = v
    default:
      fmt.Errorf("WriteBigWig(): invalid arguments")
    }
  }
  // get reduction levels for zoomed data
  reductionLevels := track.writeBigWig_reductionLevels(parameters)
  // create new bigWig writer
  writer, err := NewBigWigWriter(filename, reductionLevels, parameters)
  if err != nil {
    return err
  }
  // write data
  for name, sequence := range track.Data {
    if err := writer.Write(name, sequence, track.Binsize); err != nil {
      return err
    }
  }
  if err := writer.WriteIndex(); err != nil {
    return err
  }
  // write zoomed data
  for i := 0; i < len(reductionLevels); i++ {
    binsize := track.Binsize*reductionLevels[i]
    // save current offset as the beginning of zoomed data for reduction
    // level i
    writer.StartZoomData(i)
    for name, sequence := range track.Data {
      s := track.writeBigWig_zoomSequence(sequence, reductionLevels[i])
      if err := writer.WriteZoom(name, s, binsize, i); err != nil {
        return err
      }
    }
    // write index for this reduction level
    if err := writer.WriteZoomIndex(i); err != nil {
      return err
    }
  }
  writer.Close()

  return nil
}
