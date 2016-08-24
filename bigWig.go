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

const BIGWIG_MAGIC = 0x888FFC26

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
  // write chromosome list
  if offset, err := bwf.Fptr.Seek(0, 1); err != nil {
    return fmt.Errorf("writing `%s' failed: %v", filename, err)
  } else {
    bwf.Header.CtOffset = uint64(offset)
  }
  if err := bwf.ChromData.Write(bwf.Fptr); err != nil {
    return fmt.Errorf("writing `%s' failed: %v", filename, err)
  }
  // write data index offset
  if offset, err := bwf.Fptr.Seek(0, 1); err != nil {
    return fmt.Errorf("writing `%s' failed: %v", filename, err)
  } else {
    bwf.Header.IndexOffset = uint64(offset)
  }
  // write data index
  if err := bwf.Index.Write(bwf.Fptr); err != nil {
    return fmt.Errorf("writing `%s' failed: %v", filename, err)
  }
  // data starts here
  if offset, err := bwf.Fptr.Seek(0, 1); err != nil {
    return fmt.Errorf("writing `%s' failed: %v", filename, err)
  } else {
    bwf.Header.DataOffset = uint64(offset)
  }
  // update index size
  bwf.Index.IdxSize = bwf.Header.DataOffset - bwf.Header.IndexOffset
  bwf.Index.WriteSize(bwf.Fptr)
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
