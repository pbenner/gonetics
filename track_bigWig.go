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
import "encoding/binary"

/* -------------------------------------------------------------------------- */

func (track *Track) readBigWig_block(buffer []byte, genome Genome) error {
  reader, err := NewBbiBlockDecoder(buffer)
  if err != nil {
    return err
  }
  if idx := int(reader.Header.ChromId); idx < 0 || idx > genome.Length() {
    return fmt.Errorf("invalid chromosome id")
  }
  // convert chromosome id to sequence name
  seqname := genome.Seqnames[int(reader.Header.ChromId)]

  // allocate track if this is the first buffer
  if len(track.Data) == 0 && track.Binsize == 0 {
    *track = AllocTrack("", genome, int(reader.Header.Span))
  }
  switch reader.Header.Type {
  case 2:
    if int(reader.Header.Span) != track.Binsize {
      return fmt.Errorf("block has invalid span `%d' for track with bin size `%d'", reader.Header.Span, track.Binsize)
    }
  case 3:
    if int(reader.Header.Span) != track.Binsize {
      return fmt.Errorf("block has invalid span `%d' for track with bin size `%d'", reader.Header.Span, track.Binsize)
    }
    if int(reader.Header.Step) != track.Binsize {
      return fmt.Errorf("block has invalid step `%d' for track with bin size `%d'", reader.Header.Span, track.Binsize)
    }
  }
  if seq, ok := track.Data[seqname]; !ok {
    return fmt.Errorf("sequence `%s' not vailable in track", seqname)
  } else {
    for t := range reader.Read() {
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

func (track *Track) writeBigWig(bwf *BigWigFile, blockSize, itemsPerSlot int, fixedStep bool, genome Genome) error {
  tree := NewRTree()
  tree.BlockSize     = uint32(blockSize)
  tree.NItemsPerSlot = uint32(itemsPerSlot)
  // list of leaves
  leaves := []*RVertex{}
  // generate all leaves
  indices   := []int{}
  sequences := [][]float64{}
  for idx := 0; idx < genome.Length(); idx++ {
    name     := genome.Seqnames[idx]
    indices   = append(indices, idx)
    sequences = append(sequences, track.Data[name])
  }
  generator, _ := NewRVertexGenerator(indices, sequences, track.Binsize, blockSize, itemsPerSlot, fixedStep)
  encoder      := NewBbiBlockEncoder(track.Binsize, track.Binsize, fixedStep)

  for leaf := range generator.Read() {
    leaves = append(leaves, leaf)
    for i := 0; i < int(leaf.NChildren); i++ {
      if block, err := encoder.EncodeBlock(int(leaf.ChrIdxStart[i]), int(leaf.BaseStart[i]), leaf.Sequence[i]); err != nil {
        return err
      } else {
        if err := leaf.WriteBlock(&bwf.BbiFile, i, block); err != nil {
          return err
        }
      }
    }
  }
  if err := tree.BuildTree(leaves); err != nil {
    return err
  }
  bwf.Index = *tree
  bwf.WriteIndex()

  return nil
}

func (track *Track) WriteBigWig(filename, description string, args... interface{}) error {

  bwf        := NewBigWigFile()
  parameters := DefaultBigWigParameters()
  genome     := Genome{}

  // parse arguments
  for i := 0; i < len(args); i++ {
    switch v := args[i].(type) {
    case BigWigParameters: parameters = v
    default:
      fmt.Errorf("WriteBigWig(): invalid arguments")
    }
  }
  // identify ChromData key size and sequence lengths
  for name, seq := range track.Data {
    if bwf.ChromData.KeySize < uint32(len(name)+1) {
      bwf.ChromData.KeySize = uint32(len(name)+1)
    }
    genome = genome.AddSequence(name, len(seq)*track.Binsize)
  }
  // compress by default (this value is updated when writing blocks)
  bwf.Header.UncompressBufSize = 1
  // size of uint32
  bwf.ChromData.ValueSize = 8
  // fill ChromData
  for name, _ := range track.Data {
    key   := make([]byte, bwf.ChromData.KeySize)
    value := make([]byte, bwf.ChromData.ValueSize)
    copy(key, name)
    if idx, err := genome.GetIdx(name); err != nil {
      panic(err)
    } else {
      binary.LittleEndian.PutUint32(value[0:4], uint32(idx))
      binary.LittleEndian.PutUint32(value[4:8], uint32(genome.Lengths[idx]))
    }
    if err := bwf.ChromData.Add(key, value); err != nil {
      panic(err)
    }
  }
  if err := bwf.Create(filename); err != nil {
    return err
  }
  defer bwf.Close()
  // write data and index tree
  return track.writeBigWig(bwf, parameters.BlockSize, parameters.ItemsPerSlot, parameters.FixedStep, genome)
}
