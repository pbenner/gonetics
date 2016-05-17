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
import "errors"
import "fmt"
import "os"
import "strconv"
import "strings"

import . "github.com/pbenner/pshape/Utility"

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
    panic("NewGenome(): Invalid parameters!")
  }
  return Genome{seqnames, lengths}
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
  return 0, errors.New("sequence not found")
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
func ReadGenome(filename string) Genome {

  f, err := os.Open(filename)
  Check(err)

  // it seems that buffering the data does not increase
  // performance
  seqnames := []string{}
  lengths  := []int{}

  scanner := bufio.NewScanner(f)
  for scanner.Scan() {
    fields := strings.Fields(scanner.Text())
    if len(fields) == 0 {
      continue
    }
    if len(fields) < 2 {
      panic("Invalid genome file!")
    }
    t1, e1 := strconv.ParseInt(fields[1], 10, 64)
    Check(e1)
    seqnames = append(seqnames, fields[0])
    lengths  = append(lengths,  int(t1))
  }
  return NewGenome(seqnames, lengths)
}
