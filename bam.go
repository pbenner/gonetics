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
import "os"

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
  LSeq         int32
  NextRefID    int32
  NextPosition int32
  TLength      int32
  ReadName     string
  Cigar        uint32
  Seq          []byte
  Qual         string
}

type BamAuxiliary struct {
  Tag     [2]byte
  ValType byte
  Value   int
}

/* -------------------------------------------------------------------------- */

type BamReader struct {
  BgzfReader
  Header  BamHeader
  Genome  Genome
  Channel chan BamBlock
}

func NewBamReader(filename string) (*BamReader, error) {
  reader := new(BamReader)
  magic  := make([]byte, 4)
  genome := new(Genome)
  file   := new(os.File)
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
    genome.AddSequence(string(tmp), int(lengthSeq))
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

func (reader *BamReader) fillChannel() error {
  var blockSize int32
  for {
    // read block size
    if err := binary.Read(reader, binary.LittleEndian, &blockSize); err != nil {
      return err
    }
    // allocate new block
    block := BamBlock{}
    // read block data
    if err := binary.Read(reader, binary.LittleEndian, &block.RefID); err != nil {
      return err
    }
    if err := binary.Read(reader, binary.LittleEndian, &block.Position); err != nil {
      return err
    }
    if err := binary.Read(reader, binary.LittleEndian, &block.BinMqNl); err != nil {
      return err
    }
    if err := binary.Read(reader, binary.LittleEndian, &block.FlagNc); err != nil {
      return err
    }
    if err := binary.Read(reader, binary.LittleEndian, &block.LSeq); err != nil {
      return err
    }
    if err := binary.Read(reader, binary.LittleEndian, &block.NextRefID); err != nil {
      return err
    }
    if err := binary.Read(reader, binary.LittleEndian, &block.NextPosition); err != nil {
      return err
    }
    if err := binary.Read(reader, binary.LittleEndian, &block.TLength); err != nil {
      return err
    }
    if err := binary.Read(reader, binary.LittleEndian, &block.Cigar); err != nil {
      return err
    }
    // send block to reading thread
    reader.Channel <- block
  }
}
