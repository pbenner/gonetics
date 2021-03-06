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
import "bytes"
import "fmt"
import "io"
import "os"
import "strconv"
import "strings"

/* -------------------------------------------------------------------------- */

// Structure containing chromosome sizes.
type Genome struct {
  Seqnames []string
  Lengths  []int
}

/* constructor
 * -------------------------------------------------------------------------- */

func NewGenome(seqnames []string, lengths []int) Genome {
  if len(seqnames) != len(lengths) {
    panic("NewGenome(): invalid parameters")
  }
  return Genome{seqnames, lengths}
}

/* -------------------------------------------------------------------------- */

func (genome Genome) Clone() Genome {
  seqnames := make([]string, len(genome.Seqnames))
  lengths  := make([]int,    len(genome.Lengths))
  for i := 0; i < len(genome.Seqnames); i++ {
    seqnames[i] = genome.Seqnames[i]
    lengths [i] = genome.Lengths [i]
  }
  return NewGenome(seqnames, lengths)
}

/* -------------------------------------------------------------------------- */

// Number of chromosomes in the structure.
func (genome Genome) Length() int {
  return len(genome.Seqnames)
}

// Length of the given chromosome. Returns an error if the chromosome
// is not found.
func (genome Genome) SeqLength(seqname string) (int, error) {
  for i, s := range genome.Seqnames {
    if seqname == s {
      return genome.Lengths[i], nil
    }
  }
  return 0, fmt.Errorf("sequence `%s' not found in genome", seqname)
}

func (genome Genome) SumLengths() int {
  r := 0
  for _, l := range genome.Lengths {
    r += l
  }
  return r
}

func (genome *Genome) AddSequence(seqname string, length int) (int, error) {
  if idx, err := genome.GetIdx(seqname); err != nil {
    genome.Seqnames = append(genome.Seqnames, seqname)
    genome.Lengths  = append(genome.Lengths,  length)
    return genome.Length()-1, nil
  } else {
    return idx, fmt.Errorf("sequence `%s' already exists", seqname)
  }
}

func (genome Genome) GetIdx(seqname string) (int, error) {
  for i := 0; i < genome.Length(); i++ {
    if genome.Seqnames[i] == seqname {
      return i, nil
    }
  }
  return -1, fmt.Errorf("sequence `%s' not found in genome", seqname)
}

func (genome Genome) Filter(f func(name string, length int) bool) Genome {
  seqnames := []string{}
  lengths  := []int   {}
  for i := 0; i < genome.Length(); i++ {
    if f(genome.Seqnames[i], genome.Lengths[i]) {
      seqnames = append(seqnames, genome.Seqnames[i])
      lengths  = append(lengths,  genome.Lengths [i])
    }
  }
  return NewGenome(seqnames, lengths)
}

/* -------------------------------------------------------------------------- */

func (genome Genome) Equals(g Genome) bool {
  if genome.Length() != g.Length() {
    return false
  }
  for i, seqname := range genome.Seqnames {
    l1      := genome.Lengths[i]
    l2, err := g.SeqLength(seqname); if err != nil {
      return false
    }
    if l1 != l2 {
      return false
    }
  }
  return true
}

/* convert to string
 * -------------------------------------------------------------------------- */

func (genome Genome) String() string {
  var buffer bytes.Buffer

  printRow := func(i int) {
    if i != 0 {
      buffer.WriteString("\n")
    }
    buffer.WriteString(
      fmt.Sprintf("%10s %10d",
        genome.Seqnames[i],
        genome.Lengths [i]))
  }

  // pring header
  buffer.WriteString(
    fmt.Sprintf("%10s %10s\n", "seqnames", "lengths"))

  for i := 0; i < genome.Length(); i++ {
    printRow(i)
  }
  return buffer.String()
}

/* i/o
 * -------------------------------------------------------------------------- */

// Import chromosome sizes from a UCSC text file. The format is a whitespace
// separated table where the first column is the name of the chromosome and
// the second column the chromosome length.
func (genome *Genome) Read(r io.Reader) error {
  scanner := bufio.NewScanner(r)
  // it seems that buffering the data does not increase
  // performance
  seqnames := []string{}
  lengths  := []int{}

  for scanner.Scan() {
    fields := strings.Fields(scanner.Text())
    if len(fields) == 0 {
      continue
    }
    if len(fields) < 2 {
      return fmt.Errorf("invalid genome file")
    }
    t1, e1 := strconv.ParseInt(fields[1], 10, 64)
    if e1 != nil {
      return e1
    }
    seqnames = append(seqnames, fields[0])
    lengths  = append(lengths,  int(t1))
  }
  *genome = NewGenome(seqnames, lengths)

  return nil
}

func (genome *Genome) Import(filename string) error {

  f, err := os.Open(filename)
  if err != nil {
    return err
  }
  defer f.Close()

  if err := genome.Read(f); err != nil {
    return fmt.Errorf("reading genome from `%s' failed: %v", filename, err)
  }
  return nil
}
