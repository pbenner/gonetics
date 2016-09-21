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
import "io"
import "compress/gzip"
import "encoding/binary"

/* -------------------------------------------------------------------------- */

type BgzfExtra struct {
  SI1   uint8
  SI2   uint8
  SLen  uint16
  BSize uint16
}

type BgzfReader struct {
  gzip.Reader
}

/* -------------------------------------------------------------------------- */

func NewBgzfReader(r io.Reader) (*BgzfReader, error) {
  if reader, err := gzip.NewReader(r); err != nil {
    return nil, err
  } else {
    return &BgzfReader{Reader: *reader}, nil
  }
}

/* -------------------------------------------------------------------------- */

func (reader *BgzfReader) GetExtra() (*BgzfExtra, error) {
  if len(reader.Header.Extra) != 6 {
    return nil, fmt.Errorf("no extra information available")
  }
  extra := BgzfExtra{}
  extra.SI1   = reader.Header.Extra[0]
  extra.SI2   = reader.Header.Extra[1]
  extra.SLen  = binary.LittleEndian.Uint16(reader.Header.Extra[2:4])
  extra.BSize = binary.LittleEndian.Uint16(reader.Header.Extra[4:6])
  return &extra, nil
}
