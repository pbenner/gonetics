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

// Methods for reading and writing Big Binary Indexed files such as bigWig
// and bigBed

/* -------------------------------------------------------------------------- */

import "bytes"
import "compress/zlib"
import "fmt"
import "math"
import "encoding/binary"
import "io/ioutil"
import "os"

/* -------------------------------------------------------------------------- */

const CIRTREE_MAGIC = 0x78ca8c91
const     IDX_MAGIC = 0x2468ace0

/* -------------------------------------------------------------------------- */

func fileReadAt(file *os.File, order binary.ByteOrder, offset int64, data interface{}) error {
  currentPosition, _ := file.Seek(0, 1)
  if _, err := file.Seek(offset, 0); err != nil {
    return err
  }
  if err := binary.Read(file, order, data); err != nil {
    return err
  }
  if _, err := file.Seek(currentPosition, 0); err != nil {
    return err
  }
  return nil
}

func fileWriteAt(file *os.File, order binary.ByteOrder, offset int64, data interface{}) error {
  currentPosition, _ := file.Seek(0, 1)
  if _, err := file.Seek(offset, 0); err != nil {
    return err
  }
  if err := binary.Write(file, order, data); err != nil {
    return err
  }
  if _, err := file.Seek(currentPosition, 0); err != nil {
    return err
  }
  return nil
}

func uncompressSlice(data []byte) ([]byte, error) {
  b := bytes.NewReader(data)
  z, err := zlib.NewReader(b)
  if err != nil {
    return nil, err
  }
  defer z.Close()

  return ioutil.ReadAll(z)
}

func compressSlice(data []byte) ([]byte, error) {
  var b bytes.Buffer
  z, err := zlib.NewWriterLevel(&b, zlib.BestCompression)
  if err != nil {
    panic(err)
  }
  _, err  = z.Write(data)
  if err != nil {
    return nil, err
  }
  z.Close()

  return b.Bytes(), nil
}

/* -------------------------------------------------------------------------- */

type BbiBlockReader struct {
  Header  BbiDataHeader
  Channel chan BbiBlockReaderType
}

type BbiBlockReaderType struct {
  Idx   int
  From  int
  To    int
  Value float64
}

func NewBbiBlockReader(buffer []byte) (*BbiBlockReader, error) {
  if len(buffer) < 24 {
    return nil, fmt.Errorf("block length is shorter than 24 bytes")
  }
  reader := BbiBlockReader{}
  reader.Channel = make(chan BbiBlockReaderType)
  // parse header
  reader.Header.ReadBuffer(buffer)
  // crop header from buffer
  buffer = buffer[24:]

  switch reader.Header.Type {
  default:
    return nil, fmt.Errorf("unsupported block type")
  case 2:
    if len(buffer) % 8 != 0 {
      return nil, fmt.Errorf("variable step data block has invalid length")
    }
    go func() {
      for i := 0; i < len(buffer); i += 8 {
        r := BbiBlockReaderType{}
        r.Idx   = i
        r.From  = int(binary.LittleEndian.Uint32(buffer[i+0:i+4]))
        r.To    = r.From + int(reader.Header.Span)
        r.Value = float64(math.Float32frombits(binary.LittleEndian.Uint32(buffer[i+4:i+8])))
        reader.Channel <- r
      }
      close(reader.Channel)
    }()
  case 3:
    if len(buffer) % 4 != 0 {
      return nil, fmt.Errorf("fixed step data block has invalid length")
    }
    go func() {
      for i := 0; i < len(buffer); i += 4 {
        r := BbiBlockReaderType{}
        r.Idx   = i
        r.From  = int(reader.Header.Start + uint32(i/4)*reader.Header.Step)
        r.To    = r.From + int(reader.Header.Span)
        r.Value = float64(math.Float32frombits(binary.LittleEndian.Uint32(buffer[i:i+4])))
        reader.Channel <- r
      }
      close(reader.Channel)
    }()
  }
  return &reader, nil
}

func (reader *BbiBlockReader) Read() <- chan BbiBlockReaderType {
  return reader.Channel
}

/* -------------------------------------------------------------------------- */

type BbiSequenceSplitter struct {
  Channel chan BbiSequenceSplitterType
}

type BbiSequenceSplitterType struct {
  Idx      int
  From     int
  To       int
  Sequence []float64
}

func NewBbiSequenceSplitter(indices []int, sequences [][]float64, itemsPerSlot int, fixedStep bool) (*BbiSequenceSplitter, error) {
  if len(indices) != len(sequences) {
    return nil, fmt.Errorf("NewBbiSequenceSplitter(): invalid arguments")
  }
  splitter := BbiSequenceSplitter{}
  splitter.Channel = make(chan BbiSequenceSplitterType)
  go func() {
    splitter.fillChannel(indices, sequences, itemsPerSlot, fixedStep)
    close(splitter.Channel)
  }()
  return &splitter, nil
}

func (splitter *BbiSequenceSplitter) Read() <- chan BbiSequenceSplitterType {
  return splitter.Channel
}

func (splitter *BbiSequenceSplitter) fillChannel(indices []int, sequences [][]float64, itemsPerSlot int, fixedStep bool) error {
  for k := 0; k < len(indices); k++ {
    idx := indices[k]
    seq := sequences[k]

    for i := 0; i < len(seq); i += itemsPerSlot {
      i_from := i
      i_to   := i
      if fixedStep {
        // loop over sequence and split sequence if maximum length
        // is reached or if value is NaN
        for j := 0; j < itemsPerSlot && i_to < len(seq); i_to++ {
          if math.IsNaN(seq[i_to]) {
            // split sequence if value is NaN
            break
          }
          j++
        }
      } else {
        // loop over sequence and count the number of valid data points
        // (i.e. non-zero and not NaN)
        for j := 0; j < itemsPerSlot && i_to < len(seq); i_to++ {
          if seq[i_to] != 0.0 && !math.IsNaN(seq[i_to]) {
            // increment number of valid data points
            j++
          }
        }
      }
      splitter.Channel <- BbiSequenceSplitterType{idx, i_from, i_to, seq[i_from:i_to]}
    }
  }
  return nil
}

/* -------------------------------------------------------------------------- */

type BbiBlockWriter struct {
  Header   BbiDataHeader
  Buffer   bytes.Buffer
  tmp      []byte
  position int
}

func NewBbiBlockWriter(chromId, from, step, span int, fixedStep bool) *BbiBlockWriter {
  writer := BbiBlockWriter{}
  writer.Header.ChromId = uint32(chromId)
  writer.Header.Start   = uint32(from)
  writer.Header.End     = uint32(from)
  writer.Header.Step    = uint32(step)
  writer.Header.Span    = uint32(span)
  writer.position       = from
  if fixedStep {
    writer.Header.Type = 3
    writer.tmp = make([]byte, 4)
  } else {
    writer.Header.Type = 2
    writer.tmp = make([]byte, 8)
  }
  return &writer
}


func (writer *BbiBlockWriter) Write(values []float64) error {
  switch writer.Header.Type {
  default:
    return fmt.Errorf("unsupported block type")
  case 2:
    // variable step
    for i := 0; i < len(values); i ++ {
      if values[i] != 0.0 && !math.IsNaN(values[i]) {
        binary.LittleEndian.PutUint32(writer.tmp[0:4], math.Float32bits(float32(writer.position)))
        binary.LittleEndian.PutUint32(writer.tmp[4:8], math.Float32bits(float32(values[i])))
        if _, err := writer.Buffer.Write(writer.tmp); err != nil {
          return err
        }
        writer.Header.ItemCount++
        writer.position += int(writer.Header.Step)
      }
    }
  case 3:
    // fixed step
    for i := 0; i < len(values); i ++ {
      binary.LittleEndian.PutUint32(writer.tmp, math.Float32bits(float32(values[i])))
      if _, err := writer.Buffer.Write(writer.tmp); err != nil {
        return err
      }
      writer.Header.ItemCount++
      writer.position += int(writer.Header.Step)
    }
  }
  writer.Header.End = uint32(writer.position)

  return nil
}

func (writer *BbiBlockWriter) Bytes() []byte {
  block := make([]byte, 24)
  writer.Header.WriteBuffer(block)
  block = append(block, writer.Buffer.Bytes()...)

  return block
}

/* -------------------------------------------------------------------------- */

type BTree struct {
  KeySize       uint32
  ValueSize     uint32
  ItemsPerBlock uint32
  ItemCount     uint64
  Root          BVertex
}

type BVertex struct {
  IsLeaf       uint8
  Keys     [][]byte
  Values   [][]byte
  Children   []BVertex
}

func NewBTree(data *BData) *BTree {
  tree := BTree{}
  tree.KeySize       = data.KeySize
  tree.ValueSize     = data.ValueSize
  tree.ItemsPerBlock = data.ItemsPerBlock
  tree.ItemCount     = data.ItemCount
  // compute tree depth
  d := int(math.Ceil(math.Log(float64(data.ItemCount))/math.Log(float64(data.ItemsPerBlock))))

  tree.Root.BuildTree(data, 0, data.ItemCount, d-1)

  return &tree
}

func (vertex *BVertex) BuildTree(data *BData, from, to uint64, level int) (uint64, error) {
  // number of values below this node
  i := uint64(0)
  if level == 0 {
    vertex.IsLeaf = 1
    for nVals := uint16(0); uint32(nVals) < data.ItemsPerBlock && from+i < to; nVals++ {
      if uint32(len(data.Keys[from+i])) != data.KeySize {
        return 0, fmt.Errorf("key number `%d' has invalid size", i)
      }
      if uint32(len(data.Values[from+i])) != data.ValueSize {
        return 0, fmt.Errorf("value number `%d' has invalid size", i)
      }
      vertex.Keys   = append(vertex.Keys,   data.Keys  [from+i])
      vertex.Values = append(vertex.Values, data.Values[from+i])
      i++
    }
  } else {
    vertex.IsLeaf = 0
    for nVals := uint16(0); uint32(nVals) < data.ItemsPerBlock && from+i < to; nVals++ {
      // append first key
      vertex.Keys = append(vertex.Keys, data.Keys[from+i])
      // create new child vertex
      v := BVertex{}
      if j, err := v.BuildTree(data, from+i, to, level-1); err != nil {
        return 0, err
      } else {
        i += j
      }
      // append child
      vertex.Children = append(vertex.Children, v)
    }
  }
  return i, nil
}

func (vertex *BVertex) writeLeaf(file *os.File) error {
  padding := uint8(0)
  nVals   := uint16(len(vertex.Keys))

  if err := binary.Write(file, binary.LittleEndian, vertex.IsLeaf); err != nil {
    return err
  }
  if err := binary.Write(file, binary.LittleEndian, padding); err != nil {
    return err
  }
  if err := binary.Write(file, binary.LittleEndian, nVals); err != nil {
    return err
  }
  for i := 0; i < len(vertex.Keys); i++ {
    if err := binary.Write(file, binary.LittleEndian, vertex.Keys[i]); err != nil {
      return err
    }
    if err := binary.Write(file, binary.LittleEndian, vertex.Values[i]); err != nil {
      return err
    }
  }
  return nil
}

func (vertex *BVertex) writeIndex(file *os.File) error {
  isLeaf  := uint8(0)
  padding := uint8(0)
  nVals   := uint16(len(vertex.Keys))
  offsets := make([]int64, nVals)

  if err := binary.Write(file, binary.LittleEndian, isLeaf); err != nil {
    return err
  }
  if err := binary.Write(file, binary.LittleEndian, padding); err != nil {
    return err
  }
  if err := binary.Write(file, binary.LittleEndian, nVals); err != nil {
    return err
  }
  for i := 0; i < int(nVals); i++ {
    if err := binary.Write(file, binary.LittleEndian, vertex.Keys[i]); err != nil {
      return err
    }
    // save current file offset
    offsets[i], _ = file.Seek(0, 1)
    // offset of the ith child vertex (first set to zero)
    if err := binary.Write(file, binary.LittleEndian, uint64(0)); err != nil {
      return err
    }
  }
  // write child vertices
  for i := 0; i < int(nVals); i++ {
    // get current file offset (where the ith child vertex begins)
    offset, _ := file.Seek(0, 1)
    // and write it at the expected position 
    if err := fileWriteAt(file, binary.LittleEndian, offsets[i], uint64(offset)); err != nil {
      return err
    }
    // write ith child
    if err := vertex.Children[i].write(file); err != nil {
      return err
    }
  }
  return nil
}

func (vertex *BVertex) write(file *os.File) error {
  if vertex.IsLeaf != 0 {
    return vertex.writeLeaf(file)
  } else {
    return vertex.writeIndex(file)
  }
  return nil
}

func (tree *BTree) Write(file *os.File) error {
  magic := uint32(CIRTREE_MAGIC)

  // ItemsPerBlock has 32 bits but nVals has only 16 bits, check for overflow
  if tree.ItemsPerBlock > uint32(^uint16(0)) {
    return fmt.Errorf("ItemsPerBlock too large (maximum value is `%d')", ^uint16(0))
  }

  if err := binary.Write(file, binary.LittleEndian, magic); err != nil {
    return err
  }
  if err := binary.Write(file, binary.LittleEndian, tree.ItemsPerBlock); err != nil {
    return err
  }
  if err := binary.Write(file, binary.LittleEndian, tree.KeySize); err != nil {
    return err
  }
  if err := binary.Write(file, binary.LittleEndian, tree.ValueSize); err != nil {
    return err
  }
  if err := binary.Write(file, binary.LittleEndian, tree.ItemCount); err != nil {
    return err
  }
  // padding
  if err := binary.Write(file, binary.LittleEndian, uint64(0)); err != nil {
    return err
  }
  return tree.Root.write(file)
}

/* -------------------------------------------------------------------------- */

type BData struct {
  KeySize       uint32
  ValueSize     uint32
  ItemsPerBlock uint32
  ItemCount     uint64

  Keys   [][]byte
  Values [][]byte
}

func NewBData() *BData {
  data := BData{}
  // default values
  data.KeySize       = 0
  data.ValueSize     = 0
  data.ItemsPerBlock = 0
  data.ItemCount     = 0
  return &data
}

func (data *BData) Add(key, value []byte) error {
  if uint32(len(key)) != data.KeySize {
    return fmt.Errorf("BData.Add(): key has invalid length")
  }
  if uint32(len(value)) != data.ValueSize {
    return fmt.Errorf("BData.Add(): value has invalid length")
  }
  data.Keys   = append(data.Keys,   key)
  data.Values = append(data.Values, value)
  data.ItemsPerBlock++
  data.ItemCount++
  return nil
}

func (data *BData) readVertexLeaf(file *os.File) error {
  var nVals   uint16
  var key   []byte
  var value []byte

  if err := binary.Read(file, binary.LittleEndian, &nVals); err != nil {
    return err
  }

  for i := 0; i < int(nVals); i++ {
    key   = make([]byte, data.KeySize)
    value = make([]byte, data.ValueSize)
    if err := binary.Read(file, binary.LittleEndian, &key); err != nil {
      return err
    }
    if err := binary.Read(file, binary.LittleEndian, &value); err != nil {
      return err
    }
    data.Keys   = append(data.Keys, key)
    data.Values = append(data.Values, value)
  }
  return nil
}

func (data *BData) readVertexIndex(file *os.File) error {
  var nVals     uint16
  var key     []byte
  var position  uint64

  key = make([]byte, data.KeySize)

  if err := binary.Read(file, binary.LittleEndian, &nVals); err != nil {
    return err
  }

  for i := 0; i < int(nVals); i++ {
    if err := binary.Read(file, binary.LittleEndian, &key); err != nil {
      return err
    }
    if err := binary.Read(file, binary.LittleEndian, &position); err != nil {
      return err
    }
    // save current position and jump to child vertex
    currentPosition, _ := file.Seek(0, 1)
    if _, err := file.Seek(int64(position), 0); err != nil {
      return err
    }
    data.readVertex(file)
    if _, err := file.Seek(currentPosition, 0); err != nil {
      return err
    }
  }
  return nil
}

func (data *BData) readVertex(file *os.File) error {
  var isLeaf  uint8
  var padding uint8

  if err := binary.Read(file, binary.LittleEndian, &isLeaf); err != nil {
    return err
  }
  if err := binary.Read(file, binary.LittleEndian, &padding); err != nil {
    return err
  }
  if isLeaf != 0 {
    return data.readVertexLeaf(file)
  } else {
    return data.readVertexIndex(file)
  }
}

func (data *BData) Read(file *os.File) error {

  var magic uint32

  // magic number
  if err := binary.Read(file, binary.LittleEndian, &magic); err != nil {
    return err
  }
  if magic != CIRTREE_MAGIC {
    return fmt.Errorf("invalid tree")
  }

  if err := binary.Read(file, binary.LittleEndian, &data.ItemsPerBlock); err != nil {
    return err
  }
  if err := binary.Read(file, binary.LittleEndian, &data.KeySize); err != nil {
    return err
  }
  if err := binary.Read(file, binary.LittleEndian, &data.ValueSize); err != nil {
    return err
  }
  if err := binary.Read(file, binary.LittleEndian, &data.ItemCount); err != nil {
    return err
  }
  // padding
  if err := binary.Read(file, binary.LittleEndian, &magic); err != nil {
    return err
  }
  if err := binary.Read(file, binary.LittleEndian, &magic); err != nil {
    return err
  }
  return data.readVertex(file)
}

func (data *BData) Write(file *os.File) error {
  tree := NewBTree(data)
  return tree.Write(file)
}

/* -------------------------------------------------------------------------- */

type BbiDataHeader struct {
  ChromId   uint32
  Start     uint32
  End       uint32
  Step      uint32
  Span      uint32
  Type      byte
  Reserved  byte
  ItemCount uint16
}

func (header *BbiDataHeader) ReadBuffer(buffer []byte) {

  header.ChromId   = binary.LittleEndian.Uint32(buffer[ 0: 4])
  header.Start     = binary.LittleEndian.Uint32(buffer[ 4: 8])
  header.End       = binary.LittleEndian.Uint32(buffer[ 8:12])
  header.Step      = binary.LittleEndian.Uint32(buffer[12:16])
  header.Span      = binary.LittleEndian.Uint32(buffer[16:20])
  header.Type      = buffer[20]
  header.Reserved  = buffer[21]
  header.ItemCount = binary.LittleEndian.Uint16(buffer[22:24])

}

func (header *BbiDataHeader) WriteBuffer(buffer []byte) {

  binary.LittleEndian.PutUint32(buffer[ 0: 4], header.ChromId)
  binary.LittleEndian.PutUint32(buffer[ 4: 8], header.Start)
  binary.LittleEndian.PutUint32(buffer[ 8:12], header.End)
  binary.LittleEndian.PutUint32(buffer[12:16], header.Step)
  binary.LittleEndian.PutUint32(buffer[16:20], header.Span)
  buffer[20] = header.Type
  buffer[21] = header.Reserved
  binary.LittleEndian.PutUint16(buffer[22:24], header.ItemCount)

}

/* -------------------------------------------------------------------------- */

type RTree struct {
  BlockSize     uint32
  NItems        uint64
  ChrIdxStart   uint32
  BaseStart     uint32
  ChrIdxEnd     uint32
  BaseEnd       uint32
  IdxSize       uint64
  NItemsPerSlot uint32
  Root          *RVertex
  PtrIdxSize    int64
}

func NewRTree() *RTree {
  tree := RTree{}
  // default values
  tree.BlockSize     = 256
  tree.NItemsPerSlot = 1024
  return &tree
}

func (tree *RTree) Read(file *os.File) error {

  var magic uint32

  // magic number
  if err := binary.Read(file, binary.LittleEndian, &magic); err != nil {
    return err
  }
  if magic != IDX_MAGIC {
    return fmt.Errorf("invalid bbi tree")
  }

  if err := binary.Read(file, binary.LittleEndian, &tree.BlockSize); err != nil {
    return err
  }
  if err := binary.Read(file, binary.LittleEndian, &tree.NItems); err != nil {
    return err
  }
  if err := binary.Read(file, binary.LittleEndian, &tree.ChrIdxStart); err != nil {
    return err
  }
  if err := binary.Read(file, binary.LittleEndian, &tree.BaseStart); err != nil {
    return err
  }
  if err := binary.Read(file, binary.LittleEndian, &tree.ChrIdxEnd); err != nil {
    return err
  }
  if err := binary.Read(file, binary.LittleEndian, &tree.BaseEnd); err != nil {
    return err
  }
  // get current offset
  if offset, err := file.Seek(0, 1); err != nil {
    return err
  } else {
    tree.PtrIdxSize = offset
  }
  if err := binary.Read(file, binary.LittleEndian, &tree.IdxSize); err != nil {
    return err
  }
  if err := binary.Read(file, binary.LittleEndian, &tree.NItemsPerSlot); err != nil {
    return err
  }
  // padding
  if err := binary.Read(file, binary.LittleEndian, &magic); err != nil {
    return err
  }
  tree.Root = new(RVertex)
  tree.Root.Read(file)

  return nil
}

func (tree *RTree) WriteSize(file *os.File) error {
  return fileWriteAt(file, binary.LittleEndian, tree.PtrIdxSize, tree.IdxSize)
}

func (tree *RTree) Write(file *os.File) error {
  // magic number
  if err := binary.Write(file, binary.LittleEndian, uint32(IDX_MAGIC)); err != nil {
    return err
  }
  if err := binary.Write(file, binary.LittleEndian, tree.BlockSize); err != nil {
    return err
  }
  if err := binary.Write(file, binary.LittleEndian, tree.NItems); err != nil {
    return err
  }
  if err := binary.Write(file, binary.LittleEndian, tree.ChrIdxStart); err != nil {
    return err
  }
  if err := binary.Write(file, binary.LittleEndian, tree.BaseStart); err != nil {
    return err
  }
  if err := binary.Write(file, binary.LittleEndian, tree.ChrIdxEnd); err != nil {
    return err
  }
  if err := binary.Write(file, binary.LittleEndian, tree.BaseEnd); err != nil {
    return err
  }
  // get current offset
  if offset, err := file.Seek(0, 1); err != nil {
    return err
  } else {
    tree.PtrIdxSize = offset
  }
  if err := binary.Write(file, binary.LittleEndian, tree.IdxSize); err != nil {
    return err
  }
  if err := binary.Write(file, binary.LittleEndian, tree.NItemsPerSlot); err != nil {
    return err
  }
  // padding
  if err := binary.Write(file, binary.LittleEndian, uint32(0)); err != nil {
    return err
  }
  tree.Root.Write(file)

  return nil
}

func (tree *RTree) buildTreeRec(leaves []*RVertex, level int) (*RVertex, []*RVertex) {
  v := new(RVertex)
  n := len(leaves)
  // return if there are no leaves
  if n == 0 {
    return nil, leaves
  }
  if level == 0 {
    if n > int(tree.BlockSize) {
      n = int(tree.BlockSize)
    }
    v.NChildren   = uint16(n)
    v.Children    = leaves[0:n]
    // update free leaf set
    leaves = leaves[n:]
  } else {
    for i := 0; i < int(tree.BlockSize) && len(leaves) > 0; i++ {
      var vertex *RVertex
      vertex, leaves = tree.buildTreeRec(leaves, level-1)
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

func (tree *RTree) BuildTree(leaves []*RVertex) error {
  if len(leaves) == 0 {
    return nil
  }
  if len(leaves) == 1 {
    tree.Root = leaves[0]
  } else {
    // compute tree depth
    d := int(math.Ceil(math.Log(float64(len(leaves)))/math.Log(float64(tree.BlockSize))))
    // construct tree
    if root, leaves := tree.buildTreeRec(leaves, d-1); len(leaves) != 0 {
      panic("internal error")
    } else {
      tree.Root = root
    }
  }
  tree.ChrIdxStart = tree.Root.ChrIdxStart[0]
  tree.ChrIdxEnd   = tree.Root.ChrIdxEnd[tree.Root.NChildren-1]
  tree.BaseStart   = tree.Root.BaseStart[0]
  tree.BaseEnd     = tree.Root.BaseEnd[tree.Root.NChildren-1]
  return nil
}

/* -------------------------------------------------------------------------- */

type RVertex struct {
  IsLeaf        uint8
  NChildren     uint16
  ChrIdxStart []uint32
  BaseStart   []uint32
  ChrIdxEnd   []uint32
  BaseEnd     []uint32
  DataOffset  []uint64
  Sizes       []uint64
  Children    []*RVertex
  // positions of DataOffset and Sizes values in file
  PtrDataOffset []int64
  PtrSizes      []int64
}

func (vertex *RVertex) ReadBlock(bwf *BbiFile, i int) ([]byte, error) {
  var err error
  block := make([]byte, vertex.Sizes[i])
  if err = fileReadAt(bwf.Fptr, binary.LittleEndian, int64(vertex.DataOffset[i]), &block); err != nil {
    return nil, err
  }
  if bwf.Header.UncompressBufSize != 0 {
    if block, err = uncompressSlice(block); err != nil {
      return nil, err
    }
  }
  return block, nil
}

func (vertex *RVertex) WriteBlock(bwf *BbiFile, i int, block []byte) error {
  var err error
  if bwf.Header.UncompressBufSize != 0 {
    // update header.UncompressBufSize if block length
    // exceeds size
    if uint32(len(block)) > bwf.Header.UncompressBufSize {
      bwf.Header.UncompressBufSize = uint32(len(block))
      if err = bwf.Header.WriteUncompressBufSize(bwf.Fptr); err != nil {
        return err
      }
    }
    if block, err = compressSlice(block); err != nil {
      return err
    }
  }
  // get current offset and update DataOffset[i]
  if offset, err := bwf.Fptr.Seek(0, 1); err != nil {
    return err
  } else {
    vertex.DataOffset[i] = uint64(offset)
    // write updated value to the required position in the file
    if err = fileWriteAt(bwf.Fptr, binary.LittleEndian, int64(vertex.PtrDataOffset[i]), vertex.DataOffset[i]); err != nil {
      return err
    }
  }
  // write data
  if err = binary.Write(bwf.Fptr, binary.LittleEndian, block); err != nil {
    return err
  }
  // update size of the data block
  vertex.Sizes[i] = uint64(len(block))
  // write it to the required position in the file
  if err = fileWriteAt(bwf.Fptr, binary.LittleEndian, int64(vertex.PtrSizes[i]), vertex.Sizes[i]); err != nil {
    return err
  }
  return nil
}

func (vertex *RVertex) Read(file *os.File) error {

  var padding uint8

  if err := binary.Read(file, binary.LittleEndian, &vertex.IsLeaf); err != nil {
    return err
  }
  if err := binary.Read(file, binary.LittleEndian, &padding); err != nil {
    return err
  }
  if err := binary.Read(file, binary.LittleEndian, &vertex.NChildren); err != nil {
    return err
  }
  // allocate data
  vertex.ChrIdxStart   = make([]uint32, vertex.NChildren)
  vertex.BaseStart     = make([]uint32, vertex.NChildren)
  vertex.ChrIdxEnd     = make([]uint32, vertex.NChildren)
  vertex.BaseEnd       = make([]uint32, vertex.NChildren)
  vertex.DataOffset    = make([]uint64, vertex.NChildren)
  vertex.PtrDataOffset = make([] int64, vertex.NChildren)
  if vertex.IsLeaf != 0 {
    vertex.Sizes       = make([]uint64, vertex.NChildren)
    vertex.PtrSizes    = make([] int64, vertex.NChildren)
  } else {
    vertex.Children    = make([]*RVertex, vertex.NChildren)
  }

  for i := 0; i < int(vertex.NChildren); i++ {
    if err := binary.Read(file, binary.LittleEndian, &vertex.ChrIdxStart[i]); err != nil {
      return err
    }
    if err := binary.Read(file, binary.LittleEndian, &vertex.BaseStart[i]); err != nil {
      return err
    }
    if err := binary.Read(file, binary.LittleEndian, &vertex.ChrIdxEnd[i]); err != nil {
      return err
    }
    if err := binary.Read(file, binary.LittleEndian, &vertex.BaseEnd[i]); err != nil {
      return err
    }
    if offset, err := file.Seek(0, 1); err != nil {
      return err
    } else {
      vertex.PtrDataOffset[i] = offset
    }
    if err := binary.Read(file, binary.LittleEndian, &vertex.DataOffset[i]); err != nil {
      return err
    }
    if vertex.IsLeaf != 0 {
      if offset, err := file.Seek(0, 1); err != nil {
        return err
      } else {
        vertex.PtrSizes[i] = offset
      }
      if err := binary.Read(file, binary.LittleEndian, &vertex.Sizes[i]); err != nil {
        return err
      }
    }
  }
  if vertex.IsLeaf == 0 {
    for i := 0; i < int(vertex.NChildren); i++ {
      // seek to child position
      if _, err := file.Seek(int64(vertex.DataOffset[i]), 0); err != nil {
        return err
      }
      vertex.Children[i] = new(RVertex)
      vertex.Children[i].Read(file)
    }
  }
  return nil
}

func (vertex *RVertex) Write(file *os.File) error {

  if len(vertex.DataOffset) != int(vertex.NChildren) {
    vertex.DataOffset = make([]uint64, vertex.NChildren)
  }
  if len(vertex.Sizes) != int(vertex.NChildren) {
    vertex.Sizes = make([]uint64, vertex.NChildren)
  }
  if len(vertex.PtrDataOffset) != int(vertex.NChildren) {
    vertex.PtrDataOffset = make([]int64, vertex.NChildren)
  }
  if len(vertex.PtrSizes) != int(vertex.NChildren) {
    vertex.PtrSizes = make([]int64, vertex.NChildren)
  }

  if err := binary.Write(file, binary.LittleEndian, vertex.IsLeaf); err != nil {
    return err
  }
  // padding
  if err := binary.Write(file, binary.LittleEndian, uint8(0)); err != nil {
    return err
  }
  if err := binary.Write(file, binary.LittleEndian, vertex.NChildren); err != nil {
    return err
  }

  for i := 0; i < int(vertex.NChildren); i++ {
    if err := binary.Write(file, binary.LittleEndian, vertex.ChrIdxStart[i]); err != nil {
      return err
    }
    if err := binary.Write(file, binary.LittleEndian, vertex.BaseStart[i]); err != nil {
      return err
    }
    if err := binary.Write(file, binary.LittleEndian, vertex.ChrIdxEnd[i]); err != nil {
      return err
    }
    if err := binary.Write(file, binary.LittleEndian, vertex.BaseEnd[i]); err != nil {
      return err
    }
    // save current offset
    if offset, err := file.Seek(0, 1); err != nil {
      return err
    } else {
      vertex.PtrDataOffset[i] = offset
    }
    if err := binary.Write(file, binary.LittleEndian, vertex.DataOffset[i]); err != nil {
      return err
    }
    // save current offset
    if offset, err := file.Seek(0, 1); err != nil {
      return err
    } else {
      vertex.PtrSizes[i] = offset
    }
    if vertex.IsLeaf != 0 {
      if err := binary.Write(file, binary.LittleEndian, vertex.Sizes[i]); err != nil {
        return err
      }
    }
  }
  if vertex.IsLeaf == 0 {
    for i := 0; i < int(vertex.NChildren); i++ {
      if offset, err := file.Seek(0, 1); err != nil {
        return err
      } else {
        // save current offset
        vertex.DataOffset[i] = uint64(offset)
        // and write at the required position
        fileWriteAt(file, binary.LittleEndian, vertex.PtrDataOffset[i], vertex.DataOffset[i])
        vertex.Children[i].Write(file)
      }
    }
  }
  return nil
}

/* -------------------------------------------------------------------------- */

type RVertexGenerator struct {
  Channel chan *RVertex
}

func NewRVertexGenerator(indices []int, sequences [][]float64, binsize, blockSize, itemsPerSlot int, fixedStep bool) (*RVertexGenerator, error) {
  if len(indices) != len(sequences) {
    return nil, fmt.Errorf("NewRVertexGenerator(): invalid arguments")
  }
  generator := RVertexGenerator{}
  generator.Channel = make(chan *RVertex)
  go func() {
    generator.fillChannel(indices, sequences, binsize, blockSize, itemsPerSlot, fixedStep)
    close(generator.Channel)
  }()
  return &generator, nil
}

func (generator *RVertexGenerator) Read() <- chan *RVertex {
  return generator.Channel
}

func (generator *RVertexGenerator) fillChannel(indices []int, sequences [][]float64, binsize, blockSize, itemsPerSlot int, fixedStep bool) error {
  splitter, err := NewBbiSequenceSplitter(indices, sequences, itemsPerSlot, fixedStep)
  if err != nil {
    return err
  }
  // current leaf
  v := new(RVertex)
  v.IsLeaf = 1
  // loop over sequence chunks
  for chunk := range splitter.Read() {
    if int(v.NChildren) == blockSize {
      // vertex is full
      generator.Channel <- v
      // create new emtpy vertex
      v = new(RVertex)
      v.IsLeaf = 1
    }
    v.ChrIdxStart = append(v.ChrIdxStart, uint32(chunk.Idx))
    v.ChrIdxEnd   = append(v.ChrIdxEnd,   uint32(chunk.Idx))
    v.BaseStart   = append(v.BaseStart,   uint32(chunk.From*binsize))
    v.BaseEnd     = append(v.BaseEnd,     uint32(chunk.To  *binsize))
    v.NChildren++
  }
  if v.NChildren != 0 {
    generator.Channel <- v
  }
  return nil
}

/* -------------------------------------------------------------------------- */

type BbiHeaderZoom struct {
  ReductionLevel    uint32
  Reserved          uint32
  DataOffset        uint64
  IndexOffset       uint64
}

func (zoomHeader *BbiHeaderZoom) Read(file *os.File) error {

  if err := binary.Read(file, binary.LittleEndian, &zoomHeader.ReductionLevel); err != nil {
    return err
  }
  if err := binary.Read(file, binary.LittleEndian, &zoomHeader.Reserved); err != nil {
    return err
  }
  if err := binary.Read(file, binary.LittleEndian, &zoomHeader.DataOffset); err != nil {
    return err
  }
  if err := binary.Read(file, binary.LittleEndian, &zoomHeader.IndexOffset); err != nil {
    return err
  }
  return nil
}

func (zoomHeader *BbiHeaderZoom) Write(file *os.File) error {

  if err := binary.Write(file, binary.LittleEndian, zoomHeader.ReductionLevel); err != nil {
    return err
  }
  if err := binary.Write(file, binary.LittleEndian, zoomHeader.Reserved); err != nil {
    return err
  }
  if err := binary.Write(file, binary.LittleEndian, zoomHeader.DataOffset); err != nil {
    return err
  }
  if err := binary.Write(file, binary.LittleEndian, zoomHeader.IndexOffset); err != nil {
    return err
  }
  return nil
}

/* -------------------------------------------------------------------------- */

type BbiHeader struct {
  Magic             uint32
  Version           uint16
  ZoomLevels        uint16
  CtOffset          uint64
  DataOffset        uint64
  IndexOffset       uint64
  FieldCould        uint16
  DefinedFieldCount uint16
  SqlOffset         uint64
  SummaryOffset     uint64
  UncompressBufSize uint32
  ExtensionOffset   uint64
  NBasesCovered     uint64
  MinVal            uint64
  MaxVal            uint64
  SumData           uint64
  SumSquared        uint64
  ZoomHeaders     []BbiHeaderZoom
  // offset positions
  PtrCtOffset          int64
  PtrDataOffset        int64
  PtrIndexOffset       int64
  PtrSqlOffset         int64
  PtrSummaryOffset     int64
  PtrUncompressBufSize int64
  PtrExtensionOffset   int64
}

func NewBbiHeader() *BbiHeader {
  header := BbiHeader{}
  header.Version = 4
  return &header
}

func (header *BbiHeader) Read(file *os.File) error {

  if err := binary.Read(file, binary.LittleEndian, &header.Magic); err != nil {
    return err
  }
  if err := binary.Read(file, binary.LittleEndian, &header.Version); err != nil {
    return err
  }
  if err := binary.Read(file, binary.LittleEndian, &header.ZoomLevels); err != nil {
    return err
  }
  if offset, err := file.Seek(0, 1); err != nil {
    return err
  } else {
    header.PtrCtOffset = offset
  }
  if err := binary.Read(file, binary.LittleEndian, &header.CtOffset); err != nil {
    return err
  }
  if offset, err := file.Seek(0, 1); err != nil {
    return err
  } else {
    header.PtrDataOffset = offset
  }
  if err := binary.Read(file, binary.LittleEndian, &header.DataOffset); err != nil {
    return err
  }
  if offset, err := file.Seek(0, 1); err != nil {
    return err
  } else {
    header.PtrIndexOffset = offset
  }
  if err := binary.Read(file, binary.LittleEndian, &header.IndexOffset); err != nil {
    return err
  }
  if err := binary.Read(file, binary.LittleEndian, &header.FieldCould); err != nil {
    return err
  }
  if err := binary.Read(file, binary.LittleEndian, &header.DefinedFieldCount); err != nil {
    return err
  }
  if offset, err := file.Seek(0, 1); err != nil {
    return err
  } else {
    header.PtrSqlOffset = offset
  }
  if err := binary.Read(file, binary.LittleEndian, &header.SqlOffset); err != nil {
    return err
  }
  if offset, err := file.Seek(0, 1); err != nil {
    return err
  } else {
    header.PtrSummaryOffset = offset
  }
  if err := binary.Read(file, binary.LittleEndian, &header.SummaryOffset); err != nil {
    return err
  }
  if offset, err := file.Seek(0, 1); err != nil {
    return err
  } else {
    header.PtrUncompressBufSize = offset
  }
  if err := binary.Read(file, binary.LittleEndian, &header.UncompressBufSize); err != nil {
    return err
  }
  if offset, err := file.Seek(0, 1); err != nil {
    return err
  } else {
    header.PtrExtensionOffset = offset
  }
  if err := binary.Read(file, binary.LittleEndian, &header.ExtensionOffset); err != nil {
    return err
  }
  // zoom levels
  header.ZoomHeaders = make([]BbiHeaderZoom, header.ZoomLevels)
  for i := 0; i < int(header.ZoomLevels); i++ {
    if err := header.ZoomHeaders[i].Read(file); err != nil {
      return err
    }
  }
  // summary
  if header.SummaryOffset > 0 {
    if _, err := file.Seek(int64(header.SummaryOffset), 0); err != nil {
      return err
    }

    if err := binary.Read(file, binary.LittleEndian, &header.NBasesCovered); err != nil {
      return err
    }
    if err := binary.Read(file, binary.LittleEndian, &header.MinVal); err != nil {
      return err
    }
    if err := binary.Read(file, binary.LittleEndian, &header.MaxVal); err != nil {
      return err
    }
    if err := binary.Read(file, binary.LittleEndian, &header.SumData); err != nil {
      return err
    }
    if err := binary.Read(file, binary.LittleEndian, &header.SumSquared); err != nil {
      return err
    }
  }
  return nil
}

func (header *BbiHeader) WriteOffsets(file *os.File) error {
  if err := fileWriteAt(file, binary.LittleEndian, header.PtrCtOffset, header.CtOffset); err != nil {
    return err
  }
  if err := fileWriteAt(file, binary.LittleEndian, header.PtrDataOffset, header.DataOffset); err != nil {
    return err
  }
  if err := fileWriteAt(file, binary.LittleEndian, header.PtrIndexOffset, header.IndexOffset); err != nil {
    return err
  }
  if err := fileWriteAt(file, binary.LittleEndian, header.PtrSqlOffset, header.SqlOffset); err != nil {
    return err
  }
  if err := fileWriteAt(file, binary.LittleEndian, header.PtrExtensionOffset, header.ExtensionOffset); err != nil {
    return err
  }
  return nil
}

func (header *BbiHeader) WriteUncompressBufSize(file *os.File) error {
  return fileWriteAt(file, binary.LittleEndian, header.PtrUncompressBufSize, header.UncompressBufSize)
}

func (header *BbiHeader) Write(file *os.File) error {

  if err := binary.Write(file, binary.LittleEndian, header.Magic); err != nil {
    return err
  }
  if err := binary.Write(file, binary.LittleEndian, header.Version); err != nil {
    return err
  }
  if err := binary.Write(file, binary.LittleEndian, header.ZoomLevels); err != nil {
    return err
  }
  if offset, err := file.Seek(0, 1); err != nil {
    return err
  } else {
    header.PtrCtOffset = offset
  }
  if err := binary.Write(file, binary.LittleEndian, header.CtOffset); err != nil {
    return err
  }
  if offset, err := file.Seek(0, 1); err != nil {
    return err
  } else {
    header.PtrDataOffset = offset
  }
  if err := binary.Write(file, binary.LittleEndian, header.DataOffset); err != nil {
    return err
  }
  if offset, err := file.Seek(0, 1); err != nil {
    return err
  } else {
    header.PtrIndexOffset = offset
  }
  if err := binary.Write(file, binary.LittleEndian, header.IndexOffset); err != nil {
    return err
  }
  if err := binary.Write(file, binary.LittleEndian, header.FieldCould); err != nil {
    return err
  }
  if err := binary.Write(file, binary.LittleEndian, header.DefinedFieldCount); err != nil {
    return err
  }
  if offset, err := file.Seek(0, 1); err != nil {
    return err
  } else {
    header.PtrSqlOffset = offset
  }
  if err := binary.Write(file, binary.LittleEndian, header.SqlOffset); err != nil {
    return err
  }
  if offset, err := file.Seek(0, 1); err != nil {
    return err
  } else {
    header.PtrSummaryOffset = offset
  }
  if err := binary.Write(file, binary.LittleEndian, header.SummaryOffset); err != nil {
    return err
  }
  if offset, err := file.Seek(0, 1); err != nil {
    return err
  } else {
    header.PtrUncompressBufSize = offset
  }
  if err := binary.Write(file, binary.LittleEndian, header.UncompressBufSize); err != nil {
    return err
  }
  if offset, err := file.Seek(0, 1); err != nil {
    return err
  } else {
    header.PtrExtensionOffset = offset
  }
  if err := binary.Write(file, binary.LittleEndian, header.ExtensionOffset); err != nil {
    return err
  }
  // zoom levels
  for i := 0; i < int(header.ZoomLevels); i++ {
    if err := header.ZoomHeaders[i].Write(file); err != nil {
      return err
    }
  }
  // summary
  if header.NBasesCovered > 0 {
    // get current offset
    if offset, err := file.Seek(0, 1); err != nil {
      return err
    } else {
      header.SummaryOffset = uint64(offset)
      // write curent offset to the position of SummaryOffset
      if err := fileWriteAt(file, binary.LittleEndian, header.PtrSummaryOffset, header.SummaryOffset); err != nil {
        return err
      }
    }
    if err := binary.Write(file, binary.LittleEndian, header.NBasesCovered); err != nil {
      return err
    }
    if err := binary.Write(file, binary.LittleEndian, header.MinVal); err != nil {
      return err
    }
    if err := binary.Write(file, binary.LittleEndian, header.MaxVal); err != nil {
      return err
    }
    if err := binary.Write(file, binary.LittleEndian, header.SumData); err != nil {
      return err
    }
    if err := binary.Write(file, binary.LittleEndian, header.SumSquared); err != nil {
      return err
    }
  }

  return nil
}

/* -------------------------------------------------------------------------- */

type BbiFile struct {
  Header    BbiHeader
  ChromData BData
  Index     RTree
  Fptr      *os.File
}

func NewBbiFile() *BbiFile {
  bwf := new(BbiFile)
  bwf.Header    = *NewBbiHeader()
  bwf.ChromData = *NewBData()
  bwf.Index     = *NewRTree()
  return bwf
}

func (bwf *BbiFile) WriteNBlocks(n int) error {
  return fileWriteAt(bwf.Fptr, binary.LittleEndian, int64(bwf.Header.DataOffset), uint64(n))
}

func (bwf *BbiFile) Close() error {
  return bwf.Fptr.Close()
}
