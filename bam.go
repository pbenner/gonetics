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
import "bufio"
import "bytes"
import "encoding/binary"
import "io"
import "os"

/* -------------------------------------------------------------------------- */

type BamSeq []byte

func (seq BamSeq) String() string {
  var buffer bytes.Buffer
  writer := bufio.NewWriter(&buffer)

  t := []byte{'=', 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'}

  for i := 0; i < len(seq); i++ {
    b1 := seq[i] >> 4
    b2 := seq[i] & 0xf
    fmt.Fprintf(writer, "%c", t[b1])
    if b2 != 0 {
      fmt.Fprintf(writer, "%c", t[b2])
    }
  }
  writer.Flush()

  return buffer.String()
}

/* -------------------------------------------------------------------------- */

type BamAuxiliary struct {
  Tag   [2]byte
  Value interface{}
}

func (aux *BamAuxiliary) Read(reader io.Reader) (int, error) {
  var valueType byte
  // number of read bytes
  n := 0
  // read data
  if err := binary.Read(reader, binary.LittleEndian, &aux.Tag[0]); err != nil {
    return n, err
  }
  n += 1
  if err := binary.Read(reader, binary.LittleEndian, &aux.Tag[1]); err != nil {
    return n, err
  }
  n += 1
  if err := binary.Read(reader, binary.LittleEndian, &valueType); err != nil {
    return n, err
  }
  // three bytes read so far
  n += 1
  // read value
  switch valueType {
  case 'c':
    value := int8(0)
    if err := binary.Read(reader, binary.LittleEndian, &value); err != nil {
      return n, err
    }
    aux.Value = value; n += 1
  case 'C':
    value := uint8(0)
    if err := binary.Read(reader, binary.LittleEndian, &value); err != nil {
      return n, err
    }
    aux.Value = value; n += 1
  case 's':
    value := int16(0)
    if err := binary.Read(reader, binary.LittleEndian, &value); err != nil {
      return n, err
    }
    aux.Value = value; n += 2
  case 'S':
    value := uint16(0)
    if err := binary.Read(reader, binary.LittleEndian, &value); err != nil {
      return n, err
    }
    aux.Value = value; n += 2
  case 'i':
    value := int32(0)
    if err := binary.Read(reader, binary.LittleEndian, &value); err != nil {
      return n, err
    }
    aux.Value = value; n += 4
  case 'I':
    value := uint32(0)
    if err := binary.Read(reader, binary.LittleEndian, &value); err != nil {
      return n, err
    }
    aux.Value = value; n += 4
  case 'B':
    var t byte
    var k int32
    if err := binary.Read(reader, binary.LittleEndian, &t); err != nil {
      return n, err
    }
    n += 1
    if err := binary.Read(reader, binary.LittleEndian, &k); err != nil {
      return n, err
    }
    n += 4
    switch t {
    case 'c':
      tmp := make([]int8, k)
      for i := 0; i < int(k); i++ {
        if err := binary.Read(reader, binary.LittleEndian, &tmp[i]); err != nil {
          return n, err
        }
        n += 1
      }
      aux.Value = tmp
    case 'C':
      tmp := make([]uint8, k)
      for i := 0; i < int(k); i++ {
        if err := binary.Read(reader, binary.LittleEndian, &tmp[i]); err != nil {
          return n, err
        }
        n += 1
      }
      aux.Value = tmp
    case 's':
      tmp := make([]int16, k)
      for i := 0; i < int(k); i++ {
        if err := binary.Read(reader, binary.LittleEndian, &tmp[i]); err != nil {
          return n, err
        }
        n += 2
      }
      aux.Value = tmp
    case 'S':
      tmp := make([]uint16, k)
      for i := 0; i < int(k); i++ {
        if err := binary.Read(reader, binary.LittleEndian, &tmp[i]); err != nil {
          return n, err
        }
        n += 2
      }
      aux.Value = tmp
    case 'i':
      tmp := make([]int32, k)
      for i := 0; i < int(k); i++ {
        if err := binary.Read(reader, binary.LittleEndian, &tmp[i]); err != nil {
          return n, err
        }
        n += 4
      }
      aux.Value = tmp
    case 'I':
      tmp := make([]uint32, k)
      for i := 0; i < int(k); i++ {
        if err := binary.Read(reader, binary.LittleEndian, &tmp[i]); err != nil {
          return n, err
        }
        n += 4
      }
      aux.Value = tmp
    case 'f':
      tmp := make([]float32, k)
      for i := 0; i < int(k); i++ {
        if err := binary.Read(reader, binary.LittleEndian, &tmp[i]); err != nil {
          return n, err
        }
        n += 4
      }
      aux.Value = tmp
    default:
      return n, fmt.Errorf("invalid auxiliary array value type `%v'", t)
    }
  default:
    return n, fmt.Errorf("invalid auxiliary value type `%v'", valueType)
  }
  return n, nil
}

/* -------------------------------------------------------------------------- */

type Flag uint16

func (flag Flag) Bit(i uint) bool {
  if (flag >> i) & 1 == 1 {
    return true
  } else {
    return false
  }
}

/* -------------------------------------------------------------------------- */

type Cigar []uint32

func (cigar Cigar) String() string {
  var buffer bytes.Buffer
  writer := bufio.NewWriter(&buffer)

  t := []byte{'M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X'}

  for i := 0; i < len(cigar); i++ {
    fmt.Fprintf(writer, "%d%c", cigar[i] >> 4, t[cigar[i] & 0xf])
  }
  writer.Flush()
  return buffer.String()
}

/* -------------------------------------------------------------------------- */

type BamHeader struct {
  TextLength int32
  Text       string
  NRef       int32
}

type BamBlock struct {
  RefID        int32
  Position     int32
  Bin          uint16
  MapQ         uint8
  RNLength     uint8
  Flag         Flag
  NCigarOp     uint16
  LSeq         int32
  NextRefID    int32
  NextPosition int32
  TLength      int32
  ReadName     string
  Cigar        Cigar
  Seq          BamSeq
  Qual         []byte
  Auxiliary    []BamAuxiliary
  Error        error
}

/* -------------------------------------------------------------------------- */

func IsBamFile(filename string) (bool, error) {
  file   := new(os.File)
  reader := new(BamReader)
  magic  := make([]byte, 4)
  if f, err := os.Open(filename);  err != nil {
    return false, err
  } else {
    file = f
  }
  defer file.Close()
  if tmp, err := NewBgzfReader(file); err != nil {
    return false, nil
  } else {
    reader.BgzfReader = *tmp
  }
  if _, err := reader.Read(magic); err != nil {
    return false, nil
  }
  if string(magic) != "BAM\001" {
    return false, nil
  }
  return true, nil
}

/* -------------------------------------------------------------------------- */

type BamReader struct {
  BgzfReader
  Header  BamHeader
  Genome  Genome
  Channel chan BamBlock
}

func NewBamReader(filename string) (*BamReader, error) {
  file   := new(os.File)
  reader := new(BamReader)
  magic  := make([]byte, 4)
  genome := new(Genome)
  // temporary space for reading bytes
  var tmp []byte

  if f, err := os.Open(filename);  err != nil {
    return nil, err
  } else {
    file = f
  }
  if tmp, err := NewBgzfReader(file); err != nil {
    return nil, err
  } else {
    reader.BgzfReader = *tmp
  }
  if _, err := reader.Read(magic); err != nil {
    return nil, err
  }
  if string(magic) != "BAM\001" {
    return nil, fmt.Errorf("not a BAM file")
  }
  if err := binary.Read(reader, binary.LittleEndian, &reader.Header.TextLength); err != nil {
    return nil, err
  } else {
    tmp = make([]byte, reader.Header.TextLength)

    if _, err := reader.Read(tmp); err != nil {
      return nil, err
    }
    reader.Header.Text = string(tmp)
  }

  if err := binary.Read(reader, binary.LittleEndian, &reader.Header.NRef); err != nil {
    return nil, err
  }
  for i := 0; i < int(reader.Header.NRef); i++ {
    lengthName := int32(0)
    lengthSeq  := int32(0)
    // read length of sequence name
    if err := binary.Read(reader, binary.LittleEndian, &lengthName); err != nil {
      return nil, err
    }
    // allocate memory for reading the name
    tmp = make([]byte, lengthName)
    // read sequence name
    if _, err := reader.Read(tmp); err != nil {
      return nil, err
    }
    // read sequence length
    if err := binary.Read(reader, binary.LittleEndian, &lengthSeq); err != nil {
      return nil, err
    }
    genome.AddSequence(string(tmp[0:lengthName-1]), int(lengthSeq))
  }
  reader.Genome  = *genome
  reader.Channel = make(chan BamBlock)
  // fill channel with blocks
  go func() {
    reader.fillChannel()
    // close channel and file
    close(reader.Channel)
    file.Close()
  }()

  return reader, nil
}

func (reader *BamReader) ReadBlocks() <- chan BamBlock {
  return reader.Channel
}

func (reader *BamReader) fillChannel() {
  var blockSize int32
  var flagNc    uint32
  var binMqNl   uint32
  for {
    buf := bytes.NewBuffer([]byte{})
    // read block size
    if err := binary.Read(reader, binary.LittleEndian, &blockSize); err != nil {
      if err == io.EOF {
        return
      }
      reader.Channel <- BamBlock{Error: err}
      return
    }
    // allocate new block
    block := BamBlock{}
    // read block data
    if err := binary.Read(reader, binary.LittleEndian, &block.RefID); err != nil {
      reader.Channel <- BamBlock{Error: err}
      return
    }
    if err := binary.Read(reader, binary.LittleEndian, &block.Position); err != nil {
      reader.Channel <- BamBlock{Error: err}
      return
    }
    if err := binary.Read(reader, binary.LittleEndian, &binMqNl); err != nil {
      reader.Channel <- BamBlock{Error: err}
      return
    }
    block.Bin      = uint16((binMqNl >> 16) & 0xffff)
    block.MapQ     = uint8 ((binMqNl >>  8) & 0xff)
    block.RNLength = uint8 ((binMqNl >>  0) & 0xff)
    if err := binary.Read(reader, binary.LittleEndian, &flagNc); err != nil {
      reader.Channel <- BamBlock{Error: err}
      return
    }
    // get Flag and NCigarOp from FlagNc
    block.Flag     = Flag(flagNc >> 16)
    block.NCigarOp = uint16(flagNc & 0xffff)
    if err := binary.Read(reader, binary.LittleEndian, &block.LSeq); err != nil {
      reader.Channel <- BamBlock{Error: err}
      return
    }
    if err := binary.Read(reader, binary.LittleEndian, &block.NextRefID); err != nil {
      reader.Channel <- BamBlock{Error: err}
      return
    }
    if err := binary.Read(reader, binary.LittleEndian, &block.NextPosition); err != nil {
      reader.Channel <- BamBlock{Error: err}
      return
    }
    if err := binary.Read(reader, binary.LittleEndian, &block.TLength); err != nil {
      reader.Channel <- BamBlock{Error: err}
      return
    }
    // parse the read name
    var b byte
    for {
      if err := binary.Read(reader, binary.LittleEndian, &b); err != nil {
        reader.Channel <- BamBlock{Error: err}
        return
      }
      if b == 0 {
        block.ReadName = buf.String()
        break
      }
      buf.WriteByte(b)
    }
    // parse cigar block
    block.Cigar = make(Cigar, block.NCigarOp)
    for i := 0; i < int(block.NCigarOp); i++ {
      if err := binary.Read(reader, binary.LittleEndian, &block.Cigar[i]); err != nil {
        reader.Channel <- BamBlock{Error: err}
        return
      }
    }
    // parse seq
    block.Seq = make([]byte, (block.LSeq+1)/2)
    for i := 0; i < int((block.LSeq+1)/2); i++ {
      if err := binary.Read(reader, binary.LittleEndian, &block.Seq[i]); err != nil {
        reader.Channel <- BamBlock{Error: err}
        return
      }
    }
    // parse qual block
    block.Qual = make([]byte, block.LSeq)
    for i := 0; i < int(block.LSeq); i++ {
      if err := binary.Read(reader, binary.LittleEndian, &block.Qual[i]); err != nil {
        reader.Channel <- BamBlock{Error: err}
        return
      }
    }
    // read auxiliary data
    position := 8*4 + int(block.RNLength) + 4*int(block.NCigarOp) + int((block.LSeq + 1)/2) + int(block.LSeq)
    for i := 0; position + i < int(blockSize); {
      aux := BamAuxiliary{}
      if n, err := aux.Read(reader); err != nil {
        reader.Channel <- BamBlock{Error: err}
        return
      } else {
        i += n
      }
      block.Auxiliary = append(block.Auxiliary, aux)
    }
    // send block to reading thread
    reader.Channel <- block
  }
}
