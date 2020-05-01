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
import "compress/gzip"
import "io"
import "os"
import "strings"
import "unicode"

/* -------------------------------------------------------------------------- */

// Structure containing genomic sequences.
type StringSet map[string][]byte

/* -------------------------------------------------------------------------- */

func NewStringSet(seqnames []string, sequences [][]byte) StringSet {
  if len(seqnames) != len(sequences) {
    panic("NewStringSet(): invalid parameters")
  }
  s := make(StringSet)

  for i := 0; i < len(sequences); i++ {
    s[seqnames[i]] = sequences[i]
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
  from := iMax(          0, r.From)
  from  = iMin(len(result), from)
  to   := iMin(len(result), r.To)
  if from != r.From || to != r.To {
    return result[from:to], fmt.Errorf("GetSlice(): range `(%d,%d)' out of bounds for sequence `%s'", r.From, r.To, name)
  } else {
    return result[from:to], nil
  }
}

/* -------------------------------------------------------------------------- */

func (s StringSet) Scan(query []byte) GRanges {
  if len(s) == 0 {
    return GRanges{}
  }
  seqnames := []string{}
  from     := []int{}
  to       := []int{}
  strand   := []byte{}
  for name, sequence := range s {
    for i := 0; i+len(query)-1 < len(sequence); i++ {
      match := true
      for j := 0; j < len(query); j++ {
        if unicode.ToLower(rune(sequence[i+j])) != unicode.ToLower(rune(query[j])) {
          match = false; break
        }
      }
      if match {
        seqnames = append(seqnames, name)
        from     = append(from, i)
        to       = append(to,   i+len(query))
      }
    }
  }
  return NewGRanges(seqnames, from, to, strand)
}

/* -------------------------------------------------------------------------- */

func (s StringSet) ReadFasta(reader io.Reader) error {
  scanner := bufio.NewScanner(reader)

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

func (s StringSet) ImportFasta(filename string) error {

  var reader io.Reader
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
    reader = g
  } else {
    reader = f
  }
  return s.ReadFasta(reader)
}

func (s StringSet) WriteFasta(writer io.Writer) error {
  for name, seq := range s {
    if _, err := fmt.Fprintf(writer,  ">%s\n", name); err != nil {
      return err
    }
    for i := 0; i < len(seq); i += 80 {
      from := i
      to   := i+80
      if to >= len(seq) {
        to = len(seq)
      }
      if _, err := fmt.Fprintf(writer, "%s\n", seq[from:to]); err != nil {
        return err
      }
    }
  }
  return nil
}

func (s StringSet) ExportFasta(filename string, compress bool) error {
  var buffer bytes.Buffer

  writer := bufio.NewWriter(&buffer)
  if err := s.WriteFasta(writer); err != nil {
    return err
  }
  writer.Flush()

  return writeFile(filename, &buffer, compress)
}
