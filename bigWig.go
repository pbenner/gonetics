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
import "math"
import "encoding/binary"
import "os"
import "strings"

/* -------------------------------------------------------------------------- */

const BIGWIG_MAGIC = 0x888FFC26

/* -------------------------------------------------------------------------- */

type BigWigParameters struct {
  BlockSize    int
  ItemsPerSlot int
}

func DefaultBigWigParameters() BigWigParameters {
  return BigWigParameters{
    BlockSize: 256,
    ItemsPerSlot: 1024 }
}

/* -------------------------------------------------------------------------- */

type BigWigFile struct {
  BbiFile
}

func NewBigWigFile() *BigWigFile {
  bbiFile := *NewBbiFile()
  bbiFile.Header.Magic = BIGWIG_MAGIC
  return &BigWigFile{bbiFile}
}

func (bwf *BigWigFile) Open(filename string) error {

  // open file
  f, err := os.Open(filename)
  if err != nil {
    return err
  }
  bwf.Fptr = f

  // parse header
  if err := bwf.Header.Read(bwf.Fptr); err != nil {
    return fmt.Errorf("reading `%s' failed: %v", filename, err)
  }
  if bwf.Header.Magic != BIGWIG_MAGIC {
    return fmt.Errorf("reading `%s' failed: not a BigWig file", filename)
  }
  // parse chromosome list, which is represented as a tree
  if _, err := bwf.Fptr.Seek(int64(bwf.Header.CtOffset), 0); err != nil {
    return fmt.Errorf("reading `%s' failed: %v", filename, err)
  }
  if err := bwf.ChromData.Read(bwf.Fptr); err != nil {
    return fmt.Errorf("reading `%s' failed: %v", filename, err)
  }
  // parse data index
  if _, err := bwf.Fptr.Seek(int64(bwf.Header.IndexOffset), 0); err != nil {
    return fmt.Errorf("reading `%s' failed: %v", filename, err)
  }
  if err := bwf.Index.Read(bwf.Fptr); err != nil {
    return fmt.Errorf("reading `%s' failed: %v", filename, err)
  }
  return nil
}

func (bwf *BigWigFile) Create(filename string) error {

  // open file
  f, err := os.Create(filename)
  if err != nil {
    return err
  }
  bwf.Fptr = f

  // write header
  if err := bwf.Header.Write(bwf.Fptr); err != nil {
    return fmt.Errorf("writing `%s' failed: %v", filename, err)
  }
  // data starts here
  if offset, err := bwf.Fptr.Seek(0, 1); err != nil {
    return fmt.Errorf("writing `%s' failed: %v", filename, err)
  } else {
    bwf.Header.DataOffset = uint64(offset)
  }
  // update offsets
  if err := bwf.Header.WriteOffsets(bwf.Fptr); err != nil {
    return err
  }
  // write number of blocks (zero at the moment)
  if err := binary.Write(bwf.Fptr, binary.LittleEndian, uint64(0)); err != nil {
    return err
  }
  return nil
}

func (bwf *BigWigFile) WriteChromList() error {
  // write chromosome list
  if offset, err := bwf.Fptr.Seek(0, 1); err != nil {
    return err
  } else {
    bwf.Header.CtOffset = uint64(offset)
  }
  if err := bwf.ChromData.Write(bwf.Fptr); err != nil {
    return err
  }
  // update offsets
  if err := bwf.Header.WriteOffsets(bwf.Fptr); err != nil {
    return err
  }
  return nil
}

func (bwf *BigWigFile) WriteIndex() error {
  // write data index offset
  if offset, err := bwf.Fptr.Seek(0, 1); err != nil {
    return err
  } else {
    bwf.Header.IndexOffset = uint64(offset)
  }
  // write data index
  if err := bwf.Index.Write(bwf.Fptr); err != nil {
    return err
  }
  // update index size
  bwf.Index.IdxSize = bwf.Header.DataOffset - bwf.Header.IndexOffset
  bwf.Index.WriteSize(bwf.Fptr)
  // update offsets
  if err := bwf.Header.WriteOffsets(bwf.Fptr); err != nil {
    return err
  }
  return nil
}

/* -------------------------------------------------------------------------- */

func IsBigWigFile(filename string) (bool, error) {

  var magic uint32

  f, err := os.Open(filename)
  if err != nil {
    return false, err
  }
  defer f.Close()
  // read magic number
  if err := binary.Read(f, binary.LittleEndian, &magic); err != nil {
    return false, err
  }
  if magic != BIGWIG_MAGIC {
    return false, nil
  }
  return true, nil

}

/* -------------------------------------------------------------------------- */

type BigWigReader struct {
  Bwf     BigWigFile
  Genome  Genome
  Channel chan BigWigReaderType
}

type BigWigReaderType struct {
  Block []byte
  Error error
}

func NewBigWigReader(filename string) (*BigWigReader, error) {
  bwr := new(BigWigReader)
  bwf := new(BigWigFile)
  if err := bwf.Open(filename); err != nil {
    return nil, err
  }
  bwr.Bwf = *bwf

  seqnames := make([]string, len(bwf.ChromData.Keys))
  lengths  := make([]int,    len(bwf.ChromData.Keys))

  for i := 0; i < len(bwf.ChromData.Keys); i++ {
    if len(bwf.ChromData.Values[i]) != 8 {
      return nil, fmt.Errorf("reading `%s' failed: invalid chromosome list", filename)
    }
    idx := int(binary.LittleEndian.Uint32(bwf.ChromData.Values[i][0:4]))
    if idx >= len(bwf.ChromData.Keys) {
      return nil, fmt.Errorf("reading `%s' failed: invalid chromosome index", filename)
    }
    seqnames[idx] = strings.TrimRight(string(bwf.ChromData.Keys[i]), "\x00")
    lengths [idx] = int(binary.LittleEndian.Uint32(bwf.ChromData.Values[i][4:8]))
  }
  bwr.Genome = NewGenome(seqnames, lengths)
  // create new channel
  bwr.Channel = make(chan BigWigReaderType)
  // fill channel with blocks
  go func() {
    bwr.fillChannel(bwf.Index.Root)
    // close channel and file
    close(bwr.Channel)
    bwr.Bwf.Close()
  }()

  return bwr, nil
}

func (reader *BigWigReader) ReadBlocks() <- chan BigWigReaderType {
  return reader.Channel
}

func (reader *BigWigReader) fillChannel(vertex *RVertex) error {

  if vertex.IsLeaf != 0 {
    for i := 0; i < int(vertex.NChildren); i++ {
      if block, err := vertex.ReadBlock(&reader.Bwf.BbiFile, i); err != nil {
        reader.Channel <- BigWigReaderType{nil, err}
      } else {
        reader.Channel <- BigWigReaderType{block, nil}
      }
    }
  } else {
    for i := 0; i < int(vertex.NChildren); i++ {
      if err := reader.fillChannel(vertex.Children[i]); err != nil {
        reader.Channel <- BigWigReaderType{nil, err}
      }
    }
  }
  return nil
}

/* -------------------------------------------------------------------------- */

type BigWigWriter struct {
  Bwf        BigWigFile
  Genome     Genome
  Parameters BigWigParameters
  Leaves     []*RVertex
}

type BigWigWriterType struct {
  Seqname    string
  Sequence []float64
}

func NewBigWigWriter(filename string, parameters BigWigParameters) (*BigWigWriter, error) {
  bww := new(BigWigWriter)
  bwf := NewBigWigFile()
  // compress by default (this value is updated when writing blocks)
  bwf.Header.UncompressBufSize = 1
  // size of uint32
  bwf.ChromData.ValueSize = 8
  // open file
  if err := bwf.Create(filename); err != nil {
    return nil, err
  }
  bww.Bwf = *bwf
  bww.Parameters = parameters

  return bww, nil
}

func (bww *BigWigWriter) useFixedStep(sequence []float64) bool {
  // number of zero or NaN values
  n := 0
  // check if sequence is dense or sparse
  for i := 0; i < len(sequence); i++ {
    if math.IsNaN(sequence[i]) || sequence[i] == 0 {
      n++
    }
  }
  return n < len(sequence)/2
}

func (bww *BigWigWriter) write(idx int, sequence []float64, binsize int) error {
  // generator for RTree vertices
  var generator *RVertexGenerator
  // data block encoder
  var encoder   *BbiBlockEncoder

  if tmp, err := NewRVertexGenerator(bww.Parameters.BlockSize, bww.Parameters.ItemsPerSlot); err != nil {
    return err
  } else {
    generator = tmp
  }
  if tmp, err := NewBbiBlockEncoder(binsize, binsize); err != nil {
    return err
  } else {
    encoder = tmp
  }
  fixedStep := bww.useFixedStep(sequence)
  // split sequence into small blocks of data and write them to file
  for leaf := range generator.Generate(idx, sequence, binsize, fixedStep) {
    // write data to file
    for i := 0; i < int(leaf.NChildren); i++ {
      if block, err := encoder.EncodeBlock(int(leaf.ChrIdxStart[i]), int(leaf.BaseStart[i]), leaf.Sequence[i], fixedStep); err != nil {
        return err
      } else {
        if err := leaf.WriteBlock(&bww.Bwf.BbiFile, i, block); err != nil {
          return err
        }
      }
    }
    // save leaf for tree construction
    bww.Leaves = append(bww.Leaves, leaf)
  }
  return nil
}

func (bww *BigWigWriter) Write(seqname string, sequence []float64, binsize int) error {
  // index of the current sequence
  var idx int
  // add sequence name and length to the genome
  if tmp, err := bww.Genome.AddSequence(seqname, len(sequence)*binsize); err != nil {
    return err
  } else {
    idx = tmp
  }
  return bww.write(idx, sequence, binsize)
}

func (bww *BigWigWriter) WriteZoom(seqname string, sequence []float64, binsize int) error {
  if idx, err := bww.Genome.GetIdx(seqname); err != nil {
    return err
  } else {
    return bww.write(idx, sequence, binsize)
  }
}

func (bww *BigWigWriter) WriteIndex() error {
  tree := NewRTree()
  tree.BlockSize     = uint32(bww.Parameters.BlockSize)
  tree.NItemsPerSlot = uint32(bww.Parameters.ItemsPerSlot)
  // construct index tree
  if err := tree.BuildTree(bww.Leaves); err != nil {
    return err
  }
  // write index to file
  bww.Bwf.Index = *tree
  if err := bww.Bwf.WriteIndex(); err != nil {
    return err
  }
  return nil
}

func (bww *BigWigWriter) WriteZoomIndex(i int) error {
  tree := NewRTree()
  tree.BlockSize     = uint32(bww.Parameters.BlockSize)
  tree.NItemsPerSlot = uint32(bww.Parameters.ItemsPerSlot)
  // construct index tree
  if err := tree.BuildTree(bww.Leaves); err != nil {
    return err
  }
  // write index to file
  bww.Bwf.IndexZoom[i] = *tree
  // TODO
  // if err := bww.Bwf.WriteIndex(i int); err != nil {
  //   return err
  // }
  return nil
}

func (bww *BigWigWriter) Close() error {
  bww.WriteIndex()
  // generate chromosome list
  for _, name := range bww.Genome.Seqnames {
    if bww.Bwf.ChromData.KeySize < uint32(len(name)+1) {
      bww.Bwf.ChromData.KeySize = uint32(len(name)+1)
    }
  }
  for _, name := range bww.Genome.Seqnames {
    key   := make([]byte, bww.Bwf.ChromData.KeySize)
    value := make([]byte, bww.Bwf.ChromData.ValueSize)
    copy(key, name)
    if idx, err := bww.Genome.GetIdx(name); err != nil {
      // this should not happen
      panic(err)
    } else {
      binary.LittleEndian.PutUint32(value[0:4], uint32(idx))
      binary.LittleEndian.PutUint32(value[4:8], uint32(bww.Genome.Lengths[idx]))
    }
    if err := bww.Bwf.ChromData.Add(key, value); err != nil {
      return err
    }
  }
  // write chromosome list to file
  if err := bww.Bwf.WriteChromList(); err != nil {
    return err
  }
  bww.Bwf.Close()

  return nil
}
