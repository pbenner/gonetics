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

import "bufio"
import "compress/gzip"
import "fmt"
import "os"
import "strings"
import "unicode"

/* -------------------------------------------------------------------------- */

// Structure containing genomic sequences.
type StringSet map[string][]byte

/* -------------------------------------------------------------------------- */

func NewStringSet(seqnames []string, sequences []byte) StringSet {
  if len(seqnames) != len(sequences) {
    panic("NewStringSet(): invalid parameters")
  }
  s := make(StringSet)

  for i := 0; i < len(sequences); i++ {
    s[seqnames[i]] = sequences
  }
  return s
}

func EmptyStringSet() StringSet {
  return make(StringSet)
}

/* -------------------------------------------------------------------------- */

func (s StringSet) GetSlice(name string, r Range) ([]byte, error) {
  result, ok := s[name]
  if !ok {
    return nil, fmt.Errorf("GetSlice(): invalid sequence name")
  }
  return result[r.From:r.To], nil
}

/* -------------------------------------------------------------------------- */

func (s StringSet) ReadFasta(filename string) error {

  var scanner *bufio.Scanner
  // open file
  f, err := os.Open(filename)
  if err != nil {
    return err
  }
  defer f.Close()

  // check if file is gzipped
  if isGzip(filename) {
    g, err := gzip.NewReader(f)
    if err != nil {
      return err
    }
    defer g.Close()
    scanner = bufio.NewScanner(g)
  } else {
    scanner = bufio.NewScanner(f)
  }
  // current sequence
  name := ""
  seq  := []byte{}

  for scanner.Scan() {
    line := scanner.Text()
    if len(line) == 0 {
      continue
    }
    if line[0] == '>' {
      // save data
      if name != "" {
        s[name] = seq
      }
      // header
      fields := strings.FieldsFunc(line, func(c rune) bool {
        return unicode.IsSpace(c) || c == '>' || c == '|'
      })
      if len(fields) == 0 {
        return fmt.Errorf("ReadFasta(): invalid fasta file")
      }
      name = fields[0]
      seq  = []byte{}
    } else {
      // data
      if name == "" {
        return fmt.Errorf("ReadFasta(): invalid fasta file")
      }
      // append sequence
      seq = append(seq, line...)
    }
  }
  if name != "" {
    s[name] = seq
  }
  return nil
}
