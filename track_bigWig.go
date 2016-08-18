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

import "bytes"
import "fmt"
import "encoding/binary"
import "math"
import "strings"

/* -------------------------------------------------------------------------- */

func (track *Track) readBigWig_block(buffer []byte, genome Genome) error {
  header := BigWigDataHeader{}
  header.ReadBuffer(buffer)
  // crop header from buffer
  buffer = buffer[24:]

  if idx := int(header.ChromId); idx < 0 || idx > genome.Length() {
    return fmt.Errorf("invalid chromosome id")
  }
  r := GRangesRow{}
  r.Seqname = genome.Seqnames[int(header.ChromId)]
  r.Range.From = int(header.Start)
  r.Range.To   = int(header.End)
  if seq, err := track.GetSlice(r); err != nil {
    return err
  } else {
    switch header.Type {
    default:
      return fmt.Errorf("unsupported block type")
    case 2:
      if len(buffer) % 8 != 0{
        return fmt.Errorf("variable step data block has invalid length")
      }
      for i := 0; i < len(buffer); i += 8 {
        position := binary.LittleEndian.Uint32(buffer[i+0:i+4])
        value1   := binary.LittleEndian.Uint32(buffer[i+4:i+8])
        value2   := math.Float32frombits(value1)
        if idx := track.Index(int(position)); idx >= len(seq) {
          return fmt.Errorf("variable step data block contains invalid index")
        } else {
          seq[idx] = float64(value2)
        }
      }
    case 3:
      if 4*len(seq) != len(buffer) {
        return fmt.Errorf("fixed step data block for sequence `%s' has invalid length (length is `%d' but should be `%d')", r.Seqname, len(seq), 4*len(buffer))
      }
      for i := 0; i < len(buffer); i += 4 {
        value1  := binary.LittleEndian.Uint32(buffer[i:i+4])
        value2  := math.Float32frombits(value1)
        seq[i/4] = float64(value2)
      }
    }
  }
  return nil
}

func (track *Track) readBigWig_allBlocks(bwf *BigWigFile, vertex *RVertex, genome Genome) error {

  if vertex.IsLeaf != 0 {
    for i := 0; i < int(vertex.NChildren); i++ {
      if block, err := vertex.ReadBlock(bwf.Fptr, bwf.Header, i); err != nil {
        return err
      } else {
        if err := track.readBigWig_block(block, genome); err != nil {
          return err
        }
      }
    }
  } else {
    for i := 0; i < int(vertex.NChildren); i++ {
      if err := track.readBigWig_allBlocks(bwf, vertex.Children[i], genome); err != nil {
        return err
      }
    }
  }
  return nil
}

func (track *Track) ReadBigWig(filename, description string, binsize int) error {

  bwf := new(BigWigFile)
  if err := bwf.Open(filename); err != nil {
    return err
  }
  defer bwf.Close()

  seqnames := make([]string, len(bwf.ChromData.Keys))
  lengths  := make([]int,    len(bwf.ChromData.Keys))

  for i := 0; i < len(bwf.ChromData.Keys); i++ {
    if len(bwf.ChromData.Values[i]) != 8 {
      return fmt.Errorf("reading `%s' failed: invalid chromosome list", filename)
    }
    idx := int(binary.LittleEndian.Uint32(bwf.ChromData.Values[i][0:4]))
    if idx >= len(bwf.ChromData.Keys) {
      return fmt.Errorf("reading `%s' failed: invalid chromosome index", filename)
    }
    seqnames[idx] = strings.TrimRight(string(bwf.ChromData.Keys[i]), "\x00")
    lengths [idx] = int(binary.LittleEndian.Uint32(bwf.ChromData.Values[i][4:8]))
  }
  genome := NewGenome(seqnames, lengths)

  *track = AllocTrack(description, genome, binsize)

  if err := track.readBigWig_allBlocks(bwf, bwf.Index.Root, genome); err != nil {
    return fmt.Errorf("reading `%s' failed: %v", filename, err)
  }
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
  // create header
  header := BigWigDataHeader{}
  header.ChromId = uint32(idx)
  header.Start   = uint32(vertex.BaseStart[i])
  header.End     = uint32(vertex.BaseEnd[i])
  header.Step    = uint32(track.Binsize)
  header.Span    = uint32(track.Binsize)
  if fixedStep {
    header.Type = 3
  } else {
    header.Type = 2
  }
  // data buffer
  var buffer bytes.Buffer

  if seq, err := track.GetSlice(r); err != nil {
    return nil, err
  } else {
    switch header.Type {
    default:
      return nil, fmt.Errorf("unsupported block type")
    case 2:
      // variable step
      tmp := make([]byte, 8)
      for i := 0; i < len(seq); i ++ {
        if seq[i] != 0.0 {
          binary.LittleEndian.PutUint32(tmp[0:4], math.Float32bits(float32(i*track.Binsize)))
          binary.LittleEndian.PutUint32(tmp[4:8], math.Float32bits(float32(seq[i])))
          if _, err := buffer.Write(tmp); err != nil {
            return nil, err
          }
          header.ItemCount++
        }
      }
    case 3:
      // fixed step
      tmp := make([]byte, 4)
      for i := 0; i < len(seq); i ++ {
        binary.LittleEndian.PutUint32(tmp, math.Float32bits(float32(seq[i])))
        if _, err := buffer.Write(tmp); err != nil {
          return nil, err
        }
        header.ItemCount++
      }
    }
  }
  block := make([]byte, 24)
  header.WriteBuffer(block)
  block = append(block, buffer.Bytes()...)

  return block, nil
}

func (track *Track) writeBigWig_allBlocks(bwf *BigWigFile, vertex *RVertex, genome Genome, fixedStep bool) error {

  if vertex.IsLeaf != 0 {
    for i := 0; i < int(vertex.NChildren); i++ {
      if block, err := track.writeBigWig_block(vertex, i, genome, fixedStep); err != nil {
        return err
      } else {
        if err := vertex.WriteBlock(bwf.Fptr, bwf.Header, i, block); err != nil {
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

func (track *Track) WriteBigWig_buildRTree(blockSize, itemsPerSlot int, genome Genome, fixedStep bool) *RTree {
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

func (track *Track) WriteBigWig(filename, description string) error {

  blockSize    := 32
  itemsPerSlot := 10
  fixedStep    := true

  bwf := NewBigWigFile()

  // store indices and sequence lengths
  genome := Genome{}

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
  bwf.Index = *track.WriteBigWig_buildRTree(blockSize, itemsPerSlot, genome, fixedStep)

  if err := bwf.Create(filename); err != nil {
    return err
  }
  defer bwf.Close()
  // traverse index tree and write data
  return track.writeBigWig_allBlocks(bwf, bwf.Index.Root, genome, fixedStep)
}
