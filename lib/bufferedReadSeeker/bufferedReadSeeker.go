/* Copyright (C) 2017 Philipp Benner
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

package bufferedReadSeeker

/* -------------------------------------------------------------------------- */

import   "fmt"
import   "io"

/* -------------------------------------------------------------------------- */

type BufferedReadSeeker struct {
  reader     io.ReadSeeker
  position   int64
  offset     int64
  bufsize    int64
  buffer   []byte
}

/* -------------------------------------------------------------------------- */

func New(reader io.ReadSeeker, bufsize int) (*BufferedReadSeeker, error) {
  if bufsize <= 0 {
    return nil, fmt.Errorf("invalid buffer size")
  }
  return &BufferedReadSeeker{reader, 0, 0, 0, make([]byte, bufsize)}, nil
}

/* -------------------------------------------------------------------------- */

func (reader *BufferedReadSeeker) fillBuffer() error {
  if _, err := reader.reader.Seek(reader.position+reader.bufsize, io.SeekStart); err != nil {
    return err
  }
  if n, err := reader.reader.Read(reader.buffer); err != nil {
    return err
  } else {
    reader.position = reader.position+reader.bufsize
    reader.bufsize  = int64(n)
    reader.offset   = 0
  }
  return nil
}

func (reader *BufferedReadSeeker) Read(p []byte) (n int, err error) {
  if len(p) > len(reader.buffer) {
    // more bytes requested than can be buffered
    reader.position += int64(len(p))
    reader.bufsize   = 0
    reader.offset    = 0
    return reader.reader.Read(p)
  } else {
    if k := int64(len(p)); k <= reader.bufsize - reader.offset {
      // buffer contains requested bytes
      copy(p, reader.buffer[reader.offset:reader.offset+k])
      reader.offset += k
    } else {
      // need to receive more bytes
      n := reader.bufsize - reader.offset
      m := k - n
      // copy what is left in the buffer
      copy(p, reader.buffer[reader.offset:reader.offset+n])
      // fill buffer with new data
      if err := reader.fillBuffer(); err != nil {
        return 0, err
      }
      // copy remaining bytes
      copy(p[n:], reader.buffer[0:m])
      reader.offset += m
    }
    return len(p), nil
  }
}

func (reader *BufferedReadSeeker) Seek(offset int64, whence int) (int64, error) {
  if whence == io.SeekCurrent {
    if n, err := reader.reader.Seek(offset, whence); err != nil {
      return n, err
    }
    whence = io.SeekStart
    offset = reader.position + reader.offset + offset
  }
  if n, err := reader.reader.Seek(offset, whence); err != nil {
    return n, err
  } else {
    if n < reader.position || n >= reader.position + reader.bufsize {
      reader.bufsize  = 0
      reader.offset   = 0
      reader.position = n
    } else {
      reader.offset   = n - reader.position
    }
    return reader.position+reader.offset, err
  }
}

func (reader *BufferedReadSeeker) SetBufSize(n int) error {
  if n <= 0 {
    return fmt.Errorf("invalid buffer size")
  }
  if n <= len(reader.buffer) {
    reader.buffer = reader.buffer[0:n]
  } else {
    reader.buffer = make([]byte, n)
  }
  reader.bufsize = 0
  reader.offset  = 0
  return nil
}
