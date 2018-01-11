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
import   "strconv"
import   "os"

import   "github.com/pborman/getopt"

import . "github.com/pbenner/gonetics"

/* -------------------------------------------------------------------------- */

type Config struct {
  PrintReadName bool
  PrintCigar    bool
  PrintSequence bool
}

/* -------------------------------------------------------------------------- */

func bamView(config Config, filenameIn string) {
  var reader *BamReader
  // options for the bam reader
  options := BamReaderOptions{}
  // default options
  options.ReadName      = config.PrintReadName
  options.ReadCigar     = config.PrintCigar
  options.ReadSequence  = config.PrintSequence
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
  // print header
  fmt.Printf("%10s %15s %17s %4s",
    "Seqname", "Position", "Flag", "MapQ")
  if options.ReadCigar {
    fmt.Printf(" %20s", "Cigar")
  }
  if options.ReadName {
    fmt.Printf(" %40s", "ReadName")
  }
  if options.ReadSequence {
    fmt.Printf(" %s", "Sequence")
  }
  fmt.Println()

  for block := range reader.ReadSingleEnd() {
    if block.Error != nil {
      log.Fatal(block.Error)
    }
    fmt.Printf("%10s %15d %5d:%011s %4d",
      reader.Genome.Seqnames[block.RefID],
      block.Position,
      block.Flag,
      strconv.FormatInt(int64(block.Flag), 2),
      block.MapQ)
    
    if options.ReadCigar {
      if s := block.Cigar.String(); s == "" {
        fmt.Printf(" %20s", "-")
      } else {
        fmt.Printf(" %20s", s)
      }
    }
    if options.ReadName {
      if s := block.ReadName; s == "" {
        fmt.Printf(" %40s", "-")
      } else {
        fmt.Printf(" %40s", s)
      }
    }
    if options.ReadSequence {
      if s := block.Seq.String(); s == "" {
        fmt.Printf(" %s", "-")
      } else {
        fmt.Printf(" %s", s)
      }
    }
    fmt.Println()
  }
}

func main() {

  config  := Config{}
  config.PrintReadName = true
  config.PrintCigar    = true
  config.PrintSequence = true

  options := getopt.New()
  options.SetProgram(fmt.Sprintf("%s", os.Args[0]))

  optNoReadName := options.  BoolLong("no-read-name",  0 ,     "do not print read names")
  optNoCigar    := options.  BoolLong("no-cigar",      0 ,     "do not print cigar strings")
  optNoSequence := options.  BoolLong("no-sequence",   0 ,     "do not print read sequences")
  optHelp       := options.  BoolLong("help",         'h',     "print help")

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
  config.PrintReadName = !*optNoReadName
  config.PrintCigar    = !*optNoCigar
  config.PrintSequence = !*optNoSequence

  bamView(config, options.Args()[0])
}
