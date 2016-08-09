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
//import "io"
import "os"

/* -------------------------------------------------------------------------- */

const  BIGWIG_MAGIC = 0x888FFC26
const CIRTREE_MAGIC = 0x78ca8c91
const     IDX_MAGIC = 0x2468ace0

/* -------------------------------------------------------------------------- */

type BData struct {
  KeySize       uint32
  ValueSize     uint32
  ItemsPerBlock uint32

  Keys   [][]byte
  Values [][]byte
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
    currentPosition, _ := file.Seek(0, 0)
    file.Seek(int64(position), 0)
    data.readVertex(file)
    file.Seek(currentPosition, 0)
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
    data.readVertexLeaf(file)
  } else {
    data.readVertexIndex(file)
  }
  return nil
}

func (data *BData) Read(file *os.File) error {

  var magic uint32
  var itemCount uint64

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
  if err := binary.Read(file, binary.LittleEndian, &itemCount); err != nil {
    return err
  }
  // padding
  if err := binary.Read(file, binary.LittleEndian, &magic); err != nil {
    return err
  }
  if err := binary.Read(file, binary.LittleEndian, &magic); err != nil {
    return err
  }

  for i := 0; i < int(itemCount); i++ {
    data.readVertex(file)
  }
  return nil
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
    file.Seek(int64(header.SummaryOffset), 0)

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
  RootOffset    uint64

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
  if offset, err := file.Seek(0, 0); err != nil {
    return err
  } else {
    tree.RootOffset = uint64(offset)
  }
  return nil
}

/* -------------------------------------------------------------------------- */

type RTreeVertex struct {

  IsLeaf        uint8
  NChildren     uint16
  ChrIdxStart []uint32
  BaseStart   []uint32
  ChrIdxEnd   []uint32
  BaseEnd     []uint32
  DataOffset  []uint64
  Sizes       []uint64
  Children    []*RTreeVertex

}

func (vertex *RTreeVertex) Read(file *os.File) error {

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
  vertex.ChrIdxStart = make([]uint32, vertex.NChildren)
  vertex.BaseStart   = make([]uint32, vertex.NChildren)
  vertex.ChrIdxEnd   = make([]uint32, vertex.NChildren)
  vertex.BaseEnd     = make([]uint32, vertex.NChildren)
  vertex.DataOffset  = make([]uint64, vertex.NChildren)
  if vertex.IsLeaf != 0{
    vertex.Sizes  = make([]uint64, vertex.NChildren)
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
    if err := binary.Read(file, binary.LittleEndian, &vertex.DataOffset[i]); err != nil {
      return err
    }
    if vertex.IsLeaf != 0 {
      if err := binary.Read(file, binary.LittleEndian, &vertex.Sizes[i]); err != nil {
        return err
      }  
    }
  }
  
  return nil
}

/* -------------------------------------------------------------------------- */

func (track *Track) ReadBigWig(filename string) error {

  // open file
  f, err := os.Open(filename)
  if err != nil {
    return err
  }
  defer f.Close()

  header    := BigWigHeader{}
  chromData := BData{}
  tree      := RTree{}

  // parse header
  if err := header.Read(f); err != nil {
    return fmt.Errorf("reading `%s' failed: %v", filename, err)
  }
  // parse chromosome list, which is represented as a tree
  if _, err := f.Seek(int64(header.CtOffset), 0); err != nil {
    return fmt.Errorf("reading `%s' failed: %v", filename, err)
  }
  if err := chromData.Read(f); err != nil {
    return fmt.Errorf("reading `%s' failed: %v", filename, err)
  }
  seqnames := make([]string, len(chromData.Keys))
  lengths  := make([]int,    len(chromData.Keys))

  for i := 0; i < len(chromData.Keys); i++ {
    if len(chromData.Values[i]) != 8 {
      return fmt.Errorf("invalid chromosome list")
    }
    idx := int(binary.LittleEndian.Uint32(chromData.Values[i][0:4]))
    if idx >= len(chromData.Keys) {
      return fmt.Errorf("invalid chromosome index")
    }
    seqnames[idx] = string(chromData.Keys[i])
    lengths [idx] = int(binary.LittleEndian.Uint32(chromData.Values[i][4:8]))
  }
  genome := NewGenome(seqnames, lengths)

  *track = AllocTrack("", genome, 10)

  // parse data
  if _, err := f.Seek(int64(header.IndexOffset), 0); err != nil {
    return fmt.Errorf("reading `%s' failed: %v", filename, err)
  }
  if err := tree.Read(f); err != nil {
    return err
  }

  return nil
}
