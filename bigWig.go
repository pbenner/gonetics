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
import "io"
import "os"
import "regexp"
import "sort"
import "strings"

/* -------------------------------------------------------------------------- */

const BIGWIG_MAGIC = 0x888FFC26

/* -------------------------------------------------------------------------- */

type BigWigParameters struct {
  BlockSize         int
  ItemsPerSlot      int
  ReductionLevels []int
}

func DefaultBigWigParameters() BigWigParameters {
  return BigWigParameters{
    BlockSize      : 256,
    ItemsPerSlot   : 1024,
    ReductionLevels: nil }
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

func (bwf *BigWigFile) Open(reader io.ReadSeeker) error {
  // parse header
  if order, err := bwf.Header.Read(reader, BIGWIG_MAGIC); err != nil {
    return err
  } else {
    bwf.Order = order
  }
  if bwf.Header.Magic != BIGWIG_MAGIC {
    return fmt.Errorf("not a BigWig file")
  }
  // parse chromosome list, which is represented as a tree
  if _, err := reader.Seek(int64(bwf.Header.CtOffset), 0); err != nil {
    return err
  }
  if err := bwf.ChromData.Read(reader, bwf.Order); err != nil {
    return err
  }
  // parse data index
  if _, err := reader.Seek(int64(bwf.Header.IndexOffset), 0); err != nil {
    return err
  }
  if err := bwf.Index.Read(reader, bwf.Order); err != nil {
    return err
  }
  // parse zoom level indices
  bwf.IndexZoom = make([]RTree, bwf.Header.ZoomLevels)
  for i := 0; i < int(bwf.Header.ZoomLevels); i++ {
    if _, err := reader.Seek(int64(bwf.Header.ZoomHeaders[i].IndexOffset), 0); err != nil {
      return err
    }
    if err := bwf.IndexZoom[i].Read(reader, bwf.Order); err != nil {
      return err
    }
  }
  return nil
}

func (bwf *BigWigFile) Create(writer io.WriteSeeker) error {
  // write header
  if err := bwf.Header.Write(writer, bwf.Order); err != nil {
    return err
  }
  // data starts here
  if offset, err := writer.Seek(0, 1); err != nil {
    return err
  } else {
    bwf.Header.DataOffset = uint64(offset)
  }
  // update offsets
  if err := bwf.Header.WriteOffsets(writer, bwf.Order); err != nil {
    return err
  }
  // write number of blocks (zero at the moment)
  if err := binary.Write(writer, binary.LittleEndian, uint64(0)); err != nil {
    return err
  }
  return nil
}

func (bwf *BigWigFile) WriteChromList(writer io.WriteSeeker) error {
  // write chromosome list
  if offset, err := writer.Seek(0, 1); err != nil {
    return err
  } else {
    bwf.Header.CtOffset = uint64(offset)
  }
  if err := bwf.ChromData.Write(writer, bwf.Order); err != nil {
    return err
  }
  // update offsets
  if err := bwf.Header.WriteOffsets(writer, bwf.Order); err != nil {
    return err
  }
  return nil
}

func (bwf *BigWigFile) WriteIndex(writer io.WriteSeeker) error {
  // write data index offset
  if offset, err := writer.Seek(0, 1); err != nil {
    return err
  } else {
    bwf.Header.IndexOffset = uint64(offset)
  }
  // write data index
  if err := bwf.Index.Write(writer, bwf.Order); err != nil {
    return err
  }
  // update offsets
  if err := bwf.Header.WriteOffsets(writer, bwf.Order); err != nil {
    return err
  }
  return nil
}

func (bwf *BigWigFile) WriteIndexZoom(writer io.WriteSeeker, i int) error {
  // write data index offset
  if offset, err := writer.Seek(0, 1); err != nil {
    return err
  } else {
    bwf.Header.ZoomHeaders[i].IndexOffset = uint64(offset)
  }
  // write data index
  if err := bwf.IndexZoom[i].Write(writer, bwf.Order); err != nil {
    return err
  }
  // update offsets
  if err := bwf.Header.ZoomHeaders[i].WriteOffsets(writer, bwf.Order); err != nil {
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
  Reader  io.ReadSeeker
  Bwf     BigWigFile
  Genome  Genome
}

type BigWigReaderType struct {
  Block []byte
  Error error
}

func NewBigWigReader(reader io.ReadSeeker) (*BigWigReader, error) {
  bwr := new(BigWigReader)
  bwf := new(BigWigFile)
  if err := bwf.Open(reader); err != nil {
    return nil, err
  }
  bwr.Reader = reader
  bwr.Bwf    = *bwf

  seqnames := make([]string, len(bwf.ChromData.Keys))
  lengths  := make([]int,    len(bwf.ChromData.Keys))

  for i := 0; i < len(bwf.ChromData.Keys); i++ {
    if len(bwf.ChromData.Values[i]) != 8 {
      return nil, fmt.Errorf("invalid chromosome list")
    }
    idx := int(binary.LittleEndian.Uint32(bwf.ChromData.Values[i][0:4]))
    if idx >= len(bwf.ChromData.Keys) {
      return nil, fmt.Errorf("invalid chromosome index")
    }
    seqnames[idx] = strings.TrimRight(string(bwf.ChromData.Keys[i]), "\x00")
    lengths [idx] = int(binary.LittleEndian.Uint32(bwf.ChromData.Values[i][4:8]))
  }
  bwr.Genome = NewGenome(seqnames, lengths)

  return bwr, nil
}

func (reader *BigWigReader) ReadBlocks() <- chan BigWigReaderType {
  // create new channel
  channel := make(chan BigWigReaderType, 10)
  // fill channel with blocks
  go func() {
    reader.fillChannel(channel, reader.Bwf.Index.Root)
    // close channel and file
    close(channel)
  }()
  return channel
}

func (reader *BigWigReader) fillChannel(channel chan BigWigReaderType, vertex *RVertex) error {

  if vertex.IsLeaf != 0 {
    for i := 0; i < int(vertex.NChildren); i++ {
      if block, err := vertex.ReadBlock(reader.Reader, &reader.Bwf.BbiFile, i); err != nil {
        channel <- BigWigReaderType{nil, err}
      } else {
        channel <- BigWigReaderType{block, nil}
      }
    }
  } else {
    for i := 0; i < int(vertex.NChildren); i++ {
      if err := reader.fillChannel(channel, vertex.Children[i]); err != nil {
        channel <- BigWigReaderType{nil, err}
      }
    }
  }
  return nil
}

func (reader *BigWigReader) Query(seqRegex string, from, to, binSize int) <- chan BbiQueryType {
  channel := make(chan BbiQueryType, 100)
  done    := make(chan bool)
  go func() {
    defer close(channel)
    defer close(done)
    if r, err := regexp.Compile("^"+seqRegex+"$"); err != nil {
      return
    } else {
      for _, seqname := range reader.Genome.Seqnames {
        if !r.MatchString(seqname) {
          continue
        }
        if idx, err := reader.Genome.GetIdx(seqname); err != nil {
          return
        } else {
          if ok := reader.Bwf.query(reader.Reader, channel, done, idx, from, to, binSize); !ok {
            return
          }
        }
      }
    }
  }()
  return channel
}

func (reader *BigWigReader) QuerySlice(seqregex string, from, to int, f BinSummaryStatistics, binSize, binOverlap int, init float64) ([]float64, int, error) {
  // first collect all records
  r := []BbiSummaryRecord{}
  // a binSize of 0 means that the raw data is returned as is
  if binSize == 0 {
    for record := range reader.Query(seqregex, from, to, binSize) {
      if record.Error != nil {
        return nil, -1, record.Error
      }
      // try to determine binSize from the first record (this most likely
      // fails for bedGraph files)
      if binSize == 0 {
        if record.DataType == BbiTypeBedGraph {
          return nil, -1, fmt.Errorf("failed determine bin-size for bigWig file: data has type bedGraph")
        }
        binSize = record.To - record.From
        r = make([]BbiSummaryRecord, divIntDown(to-from, binSize))
      }
      for idx := record.From/binSize; idx < record.To/binSize; idx++ {
        if idx >= 0 && idx < len(r) {
          r[idx] = record.BbiSummaryRecord
        }
      }
    }
  } else {
    r = make([]BbiSummaryRecord, divIntDown(to-from, binSize))
    for record := range reader.Query(seqregex, from, to, binSize) {
      if record.Error != nil {
        return nil, -1, record.Error
      }
      for idx := (record.From - from)/binSize; idx < (record.To - from)/binSize; idx++ {
        if idx >= 0 && idx < len(r) {
          r[idx] = record.BbiSummaryRecord
        }
      }
    }
  }
  // convert summary records to sequence
  s := make([]float64, len(r))
  if init != 0.0 {
    for i := 0; i < len(s); i++ {
      s[i] = init
    }
  }
  if binOverlap != 0 {
    t := BbiSummaryRecord{}
    for i := 0; i < len(s); i++ {
      t.Reset()
      for j := i-binOverlap; j <= i+binOverlap; j++ {
        if j < 0 || j >= len(s) {
          continue
        }
        if r[j].Valid > 0 {
          t.AddRecord(r[j])
        }
      }
      if t.Valid > 0 {
        s[i] = f(t.Sum, t.SumSquares, t.Min, t.Max, t.Valid, float64(binSize))
      }
    }
  } else {
    for i, t := range r {
      if t.Valid > 0 {
        s[i] = f(t.Sum, t.SumSquares, t.Min, t.Max, t.Valid, float64(binSize))
      }
    }
  }
  return s, binSize, nil
}

func (reader *BigWigReader) QuerySequence(seqregex string, f BinSummaryStatistics, binSize, binOverlap int, init float64) ([]float64, int, error) {
  if seqlength, err := reader.Genome.SeqLength(seqregex); err != nil {
    return nil, -1, err
  } else {
    return reader.QuerySlice(seqregex, 0, seqlength, f, binSize, binOverlap, init)
  }
}

func (reader *BigWigReader) GetBinSize() (int, error) {
  binSize := 0
  for record := range reader.Query(".*", 0, math.MaxInt64, binSize) {
    if record.Error != nil {
      return 0, record.Error
    }
    if record.DataType == BbiTypeBedGraph {
      return 0, fmt.Errorf("failed determine bin-size for bigWig file: data has type bedGraph")
    }
    record.Quit()
    binSize = record.To - record.From
  }
  return binSize, nil
}

/* -------------------------------------------------------------------------- */

type BigWigWriter struct {
  Writer     io.WriteSeeker
  Bwf        BigWigFile
  Genome     Genome
  Parameters BigWigParameters
  generator  *RVertexGenerator
  Leaves     map[int][]*RVertex
}

type BigWigWriterType struct {
  Seqname    string
  Sequence []float64
}

func NewBigWigWriter(writer io.WriteSeeker, genome Genome, parameters BigWigParameters) (*BigWigWriter, error) {
  bww := new(BigWigWriter)
  bwf := NewBigWigFile()
  // create new leaf map
  bww.resetLeafMap()
  // create vertex generator
  if tmp, err := NewRVertexGenerator(parameters.BlockSize, parameters.ItemsPerSlot, bwf.Order); err != nil {
    return nil, err
  } else {
    bww.generator = tmp
  }
  // add zoom headers
  for i := 0; i < len(parameters.ReductionLevels); i++ {
    bwf.Header.ZoomHeaders = append(bwf.Header.ZoomHeaders,
      BbiHeaderZoom{ReductionLevel: uint32(parameters.ReductionLevels[i])})
  }
  bwf.Header.ZoomLevels = uint16(len(parameters.ReductionLevels))
  // allocate space for zoom indices
  bwf.IndexZoom = make([]RTree, len(parameters.ReductionLevels))
  // compress by default (this value is updated when writing blocks)
  bwf.Header.UncompressBufSize = 1
  // size of uint32
  bwf.ChromData.ValueSize = 8
  // open file
  if err := bwf.Create(writer); err != nil {
    return nil, err
  }
  bww.Writer = writer
  bww.Bwf    = *bwf
  bww.Genome = genome
  bww.Parameters = parameters

  return bww, nil
}

func (bww *BigWigWriter) useFixedStep(sequence []float64) bool {
  // number of zero or NaN values
  n := 0
  // check if sequence is dense or sparse
  for i := 0; i < len(sequence); i++ {
    if math.IsNaN(sequence[i]) {
      n++
    }
  }
  return n < len(sequence)/2
}

func (bww *BigWigWriter) write(idx int, sequence []float64, binSize int) (int, error) {
  // number of blocks written
  n := 0
  // determine if fixed step sizes should be used
  // (this is false if data is sparse)
  fixedStep := bww.useFixedStep(sequence)
  // split sequence into small blocks of data and write them to file
  for tmp := range bww.generator.Generate(idx, sequence, binSize, 0, fixedStep) {
    // write data to file
    for i := 0; i < int(tmp.Vertex.NChildren); i++ {
      if err := tmp.Vertex.WriteBlock(bww.Writer, &bww.Bwf.BbiFile, i, tmp.Blocks[i]); err != nil {
        return n, err
      }
      // increment number of blocks
      n++
    }
    // save leaf for tree construction
    bww.Leaves[idx] = append(bww.Leaves[idx], tmp.Vertex)
  }
  return n, nil
}

func (bww *BigWigWriter) Write(seqname string, sequence []float64, binSize int) error {
  if idx, err := bww.Genome.GetIdx(seqname); err != nil {
    return err
  } else {
    if n, err := bww.write(idx, sequence, binSize); err != nil {
      return err
    } else {
      bww.Bwf.Header.NBlocks += uint64(n)
    }
  }
  return nil
}

func (bww *BigWigWriter) writeZoom(idx int, sequence []float64, binSize, reductionLevel int) (int, error) {
  // number of blocks written
  n := 0
  // split sequence into small blocks of data and write them to file
  for tmp := range bww.generator.Generate(idx, sequence, binSize, reductionLevel, true) {
    // write data to file
    for i := 0; i < int(tmp.Vertex.NChildren); i++ {
      if err := tmp.Vertex.WriteBlock(bww.Writer, &bww.Bwf.BbiFile, i, tmp.Blocks[i]); err != nil {
        return n, err
      }
      // increment number of blocks
      n++
    }
    // save leaf for tree construction
    bww.Leaves[idx] = append(bww.Leaves[idx], tmp.Vertex)
  }
  return n, nil
}

func (bww *BigWigWriter) WriteZoom(seqname string, sequence []float64, binSize, reductionLevel, i int) error {
  if idx, err := bww.Genome.GetIdx(seqname); err != nil {
    return err
  } else {
    if n, err := bww.writeZoom(idx, sequence, binSize, reductionLevel); err != nil {
      return err
    } else {
      bww.Bwf.Header.ZoomHeaders[i].NBlocks += uint32(n)
    }
  }
  return nil
}

func (bww *BigWigWriter) getLeavesSorted() []*RVertex {
  var indices []int
  var leaves  []*RVertex
  for k := range bww.Leaves {
    indices = append(indices, k)
  }
  sort.Ints(indices)

  for _, idx := range indices {
    leaves = append(leaves, bww.Leaves[idx]...)
  }
  return leaves
}

func (bww *BigWigWriter) resetLeafMap() {
  bww.Leaves = make(map[int][]*RVertex)
}

func (bww *BigWigWriter) WriteIndex() error {
  tree := NewRTree()
  tree.BlockSize     = uint32(bww.Parameters.BlockSize)
  tree.NItemsPerSlot = uint32(bww.Parameters.ItemsPerSlot)
  // get a sorted list of leaves
  leaves := bww.getLeavesSorted()
  // construct index tree
  if err := tree.BuildTree(leaves); err != nil {
    return err
  }
  // delete leaves
  bww.resetLeafMap()
  // write index to file
  bww.Bwf.Index = *tree
  if err := bww.Bwf.WriteIndex(bww.Writer); err != nil {
    return err
  }
  return nil
}

func (bww *BigWigWriter) WriteIndexZoom(i int) error {
  tree := NewRTree()
  tree.BlockSize     = uint32(bww.Parameters.BlockSize)
  tree.NItemsPerSlot = uint32(bww.Parameters.ItemsPerSlot)
  // get a sorted list of leaves
  leaves := bww.getLeavesSorted()
  // construct index tree
  if err := tree.BuildTree(leaves); err != nil {
    return err
  }
  // delete leaves
  bww.resetLeafMap()
  // write index to file
  bww.Bwf.IndexZoom[i] = *tree
  if err := bww.Bwf.WriteIndexZoom(bww.Writer, i); err != nil {
    return err
  }
  return nil
}

func (bww *BigWigWriter) StartZoomData(i int) error {
  if offset, err := bww.Writer.Seek(0, 1); err != nil {
    return err
  } else {
    bww.Bwf.Header.ZoomHeaders[i].DataOffset = uint64(offset)
  }
  // write NBlocks
  if err := binary.Write(bww.Writer, binary.LittleEndian, bww.Bwf.Header.ZoomHeaders[i].NBlocks); err != nil {
    return err
  }
  // update offsets
  if err := bww.Bwf.Header.ZoomHeaders[i].WriteOffsets(bww.Writer, bww.Bwf.Order); err != nil {
    return err
  }
  return nil
}

func (bww *BigWigWriter) Close() error {
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
  if err := bww.Bwf.WriteChromList(bww.Writer); err != nil {
    return err
  }
  // write nblocks
  if err := bww.Bwf.Header.WriteNBlocks(bww.Writer, bww.Bwf.Order); err != nil {
    return err
  }
  for i := 0; i < len(bww.Bwf.Header.ZoomHeaders); i++ {
    if err := bww.Bwf.Header.ZoomHeaders[i].WriteNBlocks(bww.Writer, bww.Bwf.Order); err != nil {
      return err
    }
  }
  // go to the end of the file
  if _, err := bww.Writer.Seek(0, 2); err != nil {
    return err
  }
  // write magic number
  if err := binary.Write(bww.Writer, binary.LittleEndian, bww.Bwf.Header.Magic); err != nil {
    return err
  }
  return nil
}

/* utility
 * -------------------------------------------------------------------------- */

func BigWigReadGenome(reader io.ReadSeeker) (Genome, error) {
  r, err := NewBigWigReader(reader)
  if err != nil {
    return Genome{}, err
  }
  return r.Genome, nil
}

func BigWigImportGenome(filename string) (Genome, error) {
  f, err := os.Open(filename)
  if err != nil {
    return Genome{}, err
  }
  defer f.Close()

  if genome, err := BigWigReadGenome(f); err != nil {
    return genome, fmt.Errorf("importing genome from `%s' failed: %v", filename, err)
  } else {
    return genome, nil
  }
}
