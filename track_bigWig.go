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
      seq[track.Index(t.From)] = t.Value
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
  writer, err := NewBigWigWriter(filename, parameters)
  if err != nil {
    return err
  }
  for name, sequence := range track.Data {
    writer.Write(name, sequence, track.Binsize)
  }
  writer.Close()

  return nil
}
