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

package main

/* -------------------------------------------------------------------------- */

import   "fmt"
import   "log"
import   "os"

import   "github.com/pborman/getopt"

import . "github.com/pbenner/gonetics"

/* -------------------------------------------------------------------------- */

func reg2bin(beg, end int) int {
  end -= 1
  if beg>>14 == end>>14 {
    return ((1<<15)-1)/7 + (beg>>14)
  }
  if beg>>17 == end>>17 {
    return ((1<<12)-1)/7 + (beg>>17)
  }
  if beg>>20 == end>>20 {
    return ((1<<9)-1)/7  + (beg>>20)
  }
  if beg>>23 == end>>23 {
    return ((1<<6)-1)/7  + (beg>>23)
  }
  if beg>>26 == end>>26 {
    return ((1<<3)-1)/7  + (beg>>26)
  }
  return 0;
}

/* -------------------------------------------------------------------------- */

func alignmentLength(cigar BamCigar) int {
  length := 0
  for cigarBlock := range ParseCigar(cigar) {
    switch cigarBlock.Type {
    case 'M': fallthrough
    case 'D': fallthrough
    case 'N': fallthrough
    case '=': fallthrough
    case 'X':
      length += cigarBlock.N
    }
  }
  return length
}

/* -------------------------------------------------------------------------- */

func checkBin(filenameIn string, verbose bool) {
  var reader *BamReader
  // options for the bam reader
  options := BamReaderOptions{}
  // default options
  options.ReadName      = true
  options.ReadCigar     = true
  options.ReadSequence  = true
  options.ReadAuxiliary = false
  options.ReadQual      = false

  if f, err := os.Open(filenameIn); err != nil {
    log.Fatal(err)
  } else {
    defer f.Close()
    if r, err := NewBamReader(f, options); err != nil {
      log.Fatal(err)
    } else {
      reader = r
    }
  }
  i := 0

  for block := range reader.ReadSingleEnd() {
    if block.Error != nil {
      log.Fatal(block.Error)
    }
    from := int(block.Position)
    to   := int(block.Position) + alignmentLength(block.Cigar)
    bin  := reg2bin(from, to)

    if bin != int(block.Bin) {
      fmt.Printf("record `%d' has invalid bin value (expected `%d' but bin value is `%d')\n", i, bin, block.Bin)
      fmt.Printf(" -> %+v\n\n", block)
    } else {
      if verbose {
        fmt.Printf("record `%d' has correct bin value `%d'\n", i, bin)
        fmt.Printf(" -> %+v\n\n", block)
      }
    }
    i++
  }
}

func main() {
  options := getopt.New()
  options.SetProgram(fmt.Sprintf("%s", os.Args[0]))

  optHelp       := options.  BoolLong("help",        'h',     "print help")
  optVerbose    := options.  BoolLong("verbose",     'v',     "be verbose")

  options.SetParameters("<input.bam>")
  options.Parse(os.Args)

  if *optHelp {
    options.PrintUsage(os.Stdout)
    os.Exit(0)
  }
  if len(options.Args()) != 1 {
    options.PrintUsage(os.Stderr)
    os.Exit(1)
  }
  filenameIn  := options.Args()[0]

  checkBin(filenameIn, *optVerbose)
}
