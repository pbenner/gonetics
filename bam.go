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
    fmt.Fprintf(writer, "%c", t[b2])
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
  BinMqNl      uint32
  FlagNc       uint32
  Flag         uint16
  NCigarOp     uint16
  LSeq         int32
  NextRefID    int32
  NextPosition int32
  TLength      int32
  ReadName     string
  Cigar        []uint32
  Seq          BamSeq
  Qual         []byte
  Auxiliary    []BamAuxiliary
  Error        error
}

type BamAuxiliary struct {
  Tag     [2]byte
  ValType byte
  Value   int
}

/* -------------------------------------------------------------------------- */

type BamReader struct {
  BgzfReader
  File    *os.File
  Header  BamHeader
  Genome  Genome
  Channel chan BamBlock
}

func NewBamReader(filename string) (*BamReader, error) {
  reader := new(BamReader)
  magic  := make([]byte, 4)
  genome := new(Genome)
  // temporary space for reading bytes
  var tmp []byte

  if f, err := os.Open(filename);  err != nil {
    return nil, err
  } else {
    reader.File = f
  }
  if tmp, err := NewBgzfReader(reader.File); err != nil {
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
    genome.AddSequence(string(tmp), int(lengthSeq))
  }
  reader.Genome  = *genome
  reader.Channel = make(chan BamBlock)
  // fill channel with blocks
  go func() {
    reader.fillChannel()
    // close channel and file
    close(reader.Channel)
    reader.File.Close()
  }()

  return reader, nil
}

func (reader *BamReader) ReadBlocks() <- chan BamBlock {
  return reader.Channel
}

func (reader *BamReader) fillChannel() {
  var blockSize int32
  var flagNc    uint32
  for {
    buf := bytes.NewBuffer([]byte{})
    // read block size
    if err := binary.Read(reader, binary.LittleEndian, &blockSize); err != nil {
      reader.Channel <- BamBlock{Error: err}
      return
    }
    blockStart, _ := reader.File.Seek(0, 1)
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
    if err := binary.Read(reader, binary.LittleEndian, &block.BinMqNl); err != nil {
      reader.Channel <- BamBlock{Error: err}
      return
    }
    if err := binary.Read(reader, binary.LittleEndian, &flagNc); err != nil {
      reader.Channel <- BamBlock{Error: err}
      return
    }
    // get Flag and NCigarOp from FlagNc
    block.Flag     = uint16(flagNc >> 16)
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
    block.Cigar = make([]uint32, block.NCigarOp)
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
    fmt.Printf("block: %+v\n", block)
    // read auxiliary data
    for {
      position, _ := reader.File.Seek(0, 1)
      if int32(position - blockStart) >= blockSize {
        break
      }
      aux := BamAuxiliary{}
      // read data
      if err := binary.Read(reader, binary.LittleEndian, &aux.Tag[0]); err != nil {
        reader.Channel <- BamBlock{Error: err}
        return
      }
      if err := binary.Read(reader, binary.LittleEndian, &aux.Tag[1]); err != nil {
        reader.Channel <- BamBlock{Error: err}
        return
      }
      if err := binary.Read(reader, binary.LittleEndian, &aux.ValType); err != nil {
        reader.Channel <- BamBlock{Error: err}
        return
      }
      switch aux.ValType {
      case 'c':
        value := int8(0)
        if err := binary.Read(reader, binary.LittleEndian, &value); err != nil {
          reader.Channel <- BamBlock{Error: err}
          return
        }
        aux.Value = int(value)
      case 'C':
        value := uint8(0)
        if err := binary.Read(reader, binary.LittleEndian, &value); err != nil {
          reader.Channel <- BamBlock{Error: err}
          return
        }
        aux.Value = int(value)
      case 's':
        value := int16(0)
        if err := binary.Read(reader, binary.LittleEndian, &value); err != nil {
          reader.Channel <- BamBlock{Error: err}
          return
        }
        aux.Value = int(value)
      case 'S':
        value := uint16(0)
        if err := binary.Read(reader, binary.LittleEndian, &value); err != nil {
          reader.Channel <- BamBlock{Error: err}
          return
        }
        aux.Value = int(value)
      case 'i':
        value := int32(0)
        if err := binary.Read(reader, binary.LittleEndian, &value); err != nil {
          reader.Channel <- BamBlock{Error: err}
          return
        }
        aux.Value = int(value)
      case 'I':
        value := uint32(0)
        if err := binary.Read(reader, binary.LittleEndian, &value); err != nil {
          reader.Channel <- BamBlock{Error: err}
          return
        }
        aux.Value = int(value)
      default:
        reader.Channel <- BamBlock{Error: fmt.Errorf("invalid auxiliary value type `%v'", aux.ValType)}
        return
      }
      block.Auxiliary = append(block.Auxiliary, aux)
    }
    // send block to reading thread
    reader.Channel <- block
  }
}
