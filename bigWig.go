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
import "compress/zlib"
import "fmt"
import "math"
import "encoding/binary"
import "io/ioutil"
import "os"

/* -------------------------------------------------------------------------- */

const  BIGWIG_MAGIC = 0x888FFC26
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
  z := zlib.NewWriter(&b)
  _, err := z.Write(data)
  if err != nil {
    return nil, err
  }
  z.Close()

  return b.Bytes(), nil
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

type BigWigHeaderZoom struct {
  ReductionLevel    uint32
  Reserved          uint32
  DataOffset        uint64
  IndexOffset       uint64
}

func (zoomHeader *BigWigHeaderZoom) Read(file *os.File) error {

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

func (zoomHeader *BigWigHeaderZoom) Write(file *os.File) error {

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

type BigWigDataHeader struct {
  ChromId   uint32
  Start     uint32
  End       uint32
  Step      uint32
  Span      uint32
  Type      byte
  Reserved  byte
  ItemCount uint16
}

func (header *BigWigDataHeader) ReadBuffer(buffer []byte) {

  header.ChromId   = binary.LittleEndian.Uint32(buffer[ 0: 4])
  header.Start     = binary.LittleEndian.Uint32(buffer[ 4: 8])
  header.End       = binary.LittleEndian.Uint32(buffer[ 8:12])
  header.Step      = binary.LittleEndian.Uint32(buffer[12:16])
  header.Span      = binary.LittleEndian.Uint32(buffer[16:20])
  header.Type      = buffer[20]
  header.Reserved  = buffer[21]
  header.ItemCount = binary.LittleEndian.Uint16(buffer[22:24])

}

func (header *BigWigDataHeader) WriteBuffer(buffer []byte) {

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

type BigWigHeader struct {
  Version           uint16
  ZoomLevels        uint16
  CtOffset          uint64
  DataOffset        uint64
  IndexOffset       uint64
  FieldCould        uint16
  DefinedFieldCount uint16
  SqlOffset         uint64
  SummaryOffset     uint64
  BufSize           uint32
  ExtensionOffset   uint64
  NBasesCovered     uint64
  MinVal            uint64
  MaxVal            uint64
  SumData           uint64
  SumSquared        uint64
  ZoomHeaders     []BigWigHeaderZoom
  // offset positions
  PtrCtOffset        int64
  PtrDataOffset      int64
  PtrIndexOffset     int64
  PtrSqlOffset       int64
  PtrSummaryOffset   int64
  PtrExtensionOffset int64
}

func NewBigWigHeader() *BigWigHeader {
  header := BigWigHeader{}
  header.Version = 4
  return &header
}

func (header *BigWigHeader) Read(file *os.File) error {

  var magic uint32

  // magic number
  if err := binary.Read(file, binary.LittleEndian, &magic); err != nil {
    return err
  }
  if magic != BIGWIG_MAGIC {
    return fmt.Errorf("invalid header")
  }
  // header information
  if err := binary.Read(file, binary.LittleEndian, &header.Version); err != nil {
    return err
  }
  if err := binary.Read(file, binary.LittleEndian, &header.ZoomLevels); err != nil {
    return err
  }
  if err := binary.Read(file, binary.LittleEndian, &header.CtOffset); err != nil {
    return err
  }
  if err := binary.Read(file, binary.LittleEndian, &header.DataOffset); err != nil {
    return err
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
  if err := binary.Read(file, binary.LittleEndian, &header.SqlOffset); err != nil {
    return err
  }
  if err := binary.Read(file, binary.LittleEndian, &header.SummaryOffset); err != nil {
    return err
  }
  if err := binary.Read(file, binary.LittleEndian, &header.BufSize); err != nil {
    return err
  }
  if err := binary.Read(file, binary.LittleEndian, &header.ExtensionOffset); err != nil {
    return err
  }
  // zoom levels
  header.ZoomHeaders = make([]BigWigHeaderZoom, header.ZoomLevels)
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

func (header *BigWigHeader) WriteOffsets(file *os.File) error {
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

func (header *BigWigHeader) Write(file *os.File) error {

  // magic number
  if err := binary.Write(file, binary.LittleEndian, uint32(BIGWIG_MAGIC)); err != nil {
    return err
  }
  // header information
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
  if err := binary.Write(file, binary.LittleEndian, header.BufSize); err != nil {
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

type RTree struct {
  BlockSize     uint32
  NItems        uint64
  ChrIdxStart   uint32
  BaseStart     uint32
  ChrIdxEnd     uint32
  BaseEnd       uint32
  IdxSize       uint64
  NItemsPerSlot uint32
  Root          RVertex
}

func NewRTree() *RTree {
  tree := RTree{}
  // default values
  tree.BlockSize     = 256
  tree.NItems        = 0
  tree.ChrIdxStart   = 0
  tree.BaseStart     = 0
  tree.ChrIdxEnd     = 0
  tree.BaseEnd       = 0
  tree.IdxSize       = 0
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
    return fmt.Errorf("invalid bigWig tree")
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
  tree.Root.Read(file)

  return nil
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
  //tree.Root.Write(file)

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
  Children    []RVertex
  // positions of DataOffset and Sizes values in file
  PtrDataOffset []int64
  PtrSizes      []int64
}

func (vertex *RVertex) ReadBlock(file *os.File, header BigWigHeader, i int) ([]byte, error) {
  var err error
  block := make([]byte, vertex.Sizes[i])
  if err = fileReadAt(file, binary.LittleEndian, int64(vertex.DataOffset[i]), &block); err != nil {
    return nil, err
  }
  if header.BufSize != 0 {
    if block, err = uncompressSlice(block); err != nil {
      return nil, err
    }
  }
  return block, nil
}

func (vertex *RVertex) WriteBlock(file *os.File, header BigWigHeader, i int, block []byte) error {
  var err error
  if header.BufSize != 0 {
    if block, err = compressSlice(block); err != nil {
      return err
    }
  }
  // get current offset and update DataOffset[i]
  if offset, err := file.Seek(0, 1); err != nil {
    return err
  } else {
    vertex.DataOffset[i] = uint64(offset)
    // write updated value to the required position in the file
    if err = fileWriteAt(file, binary.LittleEndian, int64(vertex.PtrDataOffset[i]), vertex.DataOffset[i]); err != nil {
      return err
    }
  }
  // write data
  if err = binary.Write(file, binary.LittleEndian, block); err != nil {
    return err
  }
  // update size of the data block
  vertex.Sizes[i] = uint64(len(block))
  // write it to the required position in the file
  if err = fileWriteAt(file, binary.LittleEndian, int64(vertex.PtrSizes[i]), vertex.Sizes[i]); err != nil {
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
    vertex.Children    = make([]RVertex, vertex.NChildren)
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
      vertex.Children[i].Read(file)
    }
  }
  return nil
}

func (vertex *RVertex) Write(file *os.File) error {

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

type BigWigFile struct {
  Header    BigWigHeader
  ChromData BData
  Index     RTree
  Fptr      *os.File
}

func NewBigWigFile() *BigWigFile {
  bwf := new(BigWigFile)
  bwf.Header    = *NewBigWigHeader()
  bwf.ChromData = *NewBData()
  bwf.Index     = *NewRTree()
  return bwf
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
  // write chromosome list
  if offset, err := bwf.Fptr.Seek(0, 1); err != nil {
    return fmt.Errorf("writing `%s' failed: %v", filename, err)
  } else {
    bwf.Header.CtOffset = uint64(offset)
  }
  if err := bwf.ChromData.Write(bwf.Fptr); err != nil {
    return fmt.Errorf("writing `%s' failed: %v", filename, err)
  }
  // write data index
  if offset, err := bwf.Fptr.Seek(0, 1); err != nil {
    return fmt.Errorf("writing `%s' failed: %v", filename, err)
  } else {
    bwf.Header.IndexOffset = uint64(offset)
  }
  if err := bwf.Index.Write(bwf.Fptr); err != nil {
    return fmt.Errorf("writing `%s' failed: %v", filename, err)
  }
  // data starts here
  if offset, err := bwf.Fptr.Seek(0, 1); err != nil {
    return fmt.Errorf("writing `%s' failed: %v", filename, err)
  } else {
    bwf.Header.DataOffset = uint64(offset)
  }
  // update offsets
  bwf.Header.WriteOffsets(bwf.Fptr)
  return nil
}

func (bwf *BigWigFile) Close() error {
  return bwf.Fptr.Close()
}
