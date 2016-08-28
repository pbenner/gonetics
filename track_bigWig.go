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
import "math"

/* -------------------------------------------------------------------------- */

type BigWigParameters struct {
  BlockSize    int
  ItemsPerSlot int
  FixedStep    bool
}

func DefaultBigWigParameters() BigWigParameters {
  return BigWigParameters{
    BlockSize: 256,
    ItemsPerSlot: 1024,
    FixedStep: true }
}

/* -------------------------------------------------------------------------- */

func (track *Track) readBigWig_block(buffer []byte, genome Genome) error {
  reader, err := NewBbiBlockReader(buffer)
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

func (track *Track) writeBigWig_block(vertex *RVertex, i int, genome Genome, fixedStep bool) ([]byte, error) {
  // get sequence index
  idx := vertex.ChrIdxStart[i]
  // create genomic range
  r := GRangesRow{}
  r.Seqname = genome.Seqnames[idx]
  r.Range.From = int(vertex.BaseStart[i])
  r.Range.To   = int(vertex.BaseEnd[i])
  // new block writer
  writer := NewBbiBlockWriter(int(idx), int(vertex.BaseStart[i]), track.Binsize, track.Binsize, fixedStep)

  if seq, err := track.GetSlice(r); err != nil {
    return nil, err
  } else {
    if err := writer.Write(seq); err != nil {
      return nil, err
    }
  }
  return writer.Bytes(), nil
}

func (track *Track) writeBigWig_allBlocks(bwf *BigWigFile, vertex *RVertex, genome Genome, fixedStep bool) error {

  if vertex.IsLeaf != 0 {
    for i := 0; i < int(vertex.NChildren); i++ {
      if block, err := track.writeBigWig_block(vertex, i, genome, fixedStep); err != nil {
        return err
      } else {
        if err := vertex.WriteBlock(&bwf.BbiFile, i, block); err != nil {
          return err
        }
      }
    }
  } else {
    for i := 0; i < int(vertex.NChildren); i++ {
      if err := track.writeBigWig_allBlocks(bwf, vertex.Children[i], genome, fixedStep); err != nil {
        return err
      }
    }
  }
  return nil
}

func (track *Track) WriteBigWig_buildRTreeRec(leaves []*RVertex, blockSize, level int) (*RVertex, []*RVertex) {
  v := new(RVertex)
  n := len(leaves)
  // return if there are no leaves
  if n == 0 {
    return nil, leaves
  }
  if level == 0 {
    if n > blockSize {
      n = blockSize
    }
    v.NChildren   = uint16(n)
    v.Children    = leaves[0:n]
    // update free leaf set
    leaves = leaves[n:]
  } else {
    for i := 0; i < blockSize && len(leaves) > 0; i++ {
      var vertex *RVertex
      vertex, leaves = track.WriteBigWig_buildRTreeRec(leaves, blockSize, level-1)
      v.NChildren++
      v.Children = append(v.Children, vertex)
    }
  }
  for i := 0; i < len(v.Children); i++ {
    v.ChrIdxStart = append(v.ChrIdxStart, v.Children[i].ChrIdxStart[0])
    v.ChrIdxEnd   = append(v.ChrIdxEnd,   v.Children[i].ChrIdxEnd[v.Children[i].NChildren-1])
    v.BaseStart   = append(v.BaseStart,   v.Children[i].BaseStart[0])
    v.BaseEnd     = append(v.BaseEnd,     v.Children[i].BaseEnd[v.Children[i].NChildren-1])
  }
  return v, leaves
}

func (track *Track) WriteBigWig_buildRTree(blockSize, itemsPerSlot int, fixedStep bool, genome Genome) *RTree {
  tree := NewRTree()
  tree.BlockSize     = uint32(blockSize)
  tree.NItemsPerSlot = uint32(itemsPerSlot)
  // list of leaves
  leaves := []*RVertex{}
  // current leaf
  v := new(RVertex)
  v.IsLeaf = 1
  // generate all leaves
  for idx := 0; idx < genome.Length(); idx++ {
    name := genome.Seqnames[idx]
    seq  := track.Data[name]
    for i := 0; i < len(seq); i += itemsPerSlot {
      if int(v.NChildren) == blockSize {
        // vertex is full
        leaves = append(leaves, v)
        // create new emtpy vertex
        v = new(RVertex)
        v.IsLeaf = 1
      }
      i_from := i
      i_to   := i
      if fixedStep {
        i_to = i+itemsPerSlot
        if i_to > len(seq) {
          i_to = len(seq)
        }
      } else {
        for n := 0; n < itemsPerSlot && i_to < len(seq); i_to++ {
          if seq[i_to] != 0.0 {
            n++
          }
        }
      }
      v.ChrIdxStart = append(v.ChrIdxStart, uint32(idx))
      v.ChrIdxEnd   = append(v.ChrIdxEnd, uint32(idx))
      v.BaseStart   = append(v.BaseStart, uint32(i_from*track.Binsize))
      v.BaseEnd     = append(v.BaseEnd, uint32(i_to*track.Binsize))
      v.NChildren++
      tree.NItems++
    }
  }
  if v.NChildren != 0 {
    leaves = append(leaves, v)
    // create new emtpy vertex
    v = new(RVertex)
    v.IsLeaf = 1
  }
  if len(leaves) == 0 {
    return tree
  }
  if len(leaves) == 1 {
    tree.Root = leaves[0]
  } else {
    // compute tree depth
    d := int(math.Ceil(math.Log(float64(len(leaves)))/math.Log(float64(blockSize))))
    // construct tree
    tree.Root, _ = track.WriteBigWig_buildRTreeRec(leaves, blockSize, d-1)
  }
  tree.ChrIdxStart = tree.Root.ChrIdxStart[0]
  tree.ChrIdxEnd   = tree.Root.ChrIdxEnd[tree.Root.NChildren-1]
  tree.BaseStart   = tree.Root.BaseStart[0]
  tree.BaseEnd     = tree.Root.BaseEnd[tree.Root.NChildren-1]

  return tree
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
  // construct index tree
  bwf.Index = *track.WriteBigWig_buildRTree(parameters.BlockSize, parameters.ItemsPerSlot, parameters.FixedStep, genome)

  if err := bwf.Create(filename); err != nil {
    return err
  }
  defer bwf.Close()
  // traverse index tree and write data
  return track.writeBigWig_allBlocks(bwf, bwf.Index.Root, genome, parameters.FixedStep)
}
