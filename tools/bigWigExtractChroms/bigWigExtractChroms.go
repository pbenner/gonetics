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
import   "math"
import   "os"
import   "strings"

import   "github.com/pborman/getopt"

import . "github.com/pbenner/gonetics"

/* -------------------------------------------------------------------------- */

func divIntUp(a, b int) int {
  return (a+b-1)/b
}

func extract(chromNames []string, filenameIn, filenameOut string, verbose bool) {
  f, err := os.Open(filenameIn); if err != nil {
    log.Fatal(err)
  }
  defer f.Close()

  if reader, err := NewBigWigReader(f); err != nil {
    log.Fatal(err)
  } else {
    binSize   := 0
    sequences := [][]float64{}
    genome    := Genome{}
    for _, chromName := range chromNames {
      if length, err := reader.Genome.SeqLength(chromName); err != nil {
        log.Fatal(err)
      } else {
        if verbose {
          log.Printf("reading sequence %s", chromName)
        }
        if s, b, err := reader.QuerySequence(chromName, BinMean, binSize, 0, math.NaN()); err != nil {
          log.Fatal (err)
        } else {
          binSize   = b
          sequences = append(sequences, s)
        }
        genome.AddSequence(chromName, length)
      }
    }
    t, _ := NewSimpleTrack("", sequences, genome, binSize)
    if verbose {
      log.Printf("writing track `%s'", filenameOut)
    }
    t.ExportBigWig(filenameOut)
  }
}

func main() {
  options := getopt.New()
  options.SetProgram(fmt.Sprintf("%s", os.Args[0]))

  optHelp    := options.BoolLong("help",    'h',  "print help")
  optVerbose := options.BoolLong("verbose", 'v',  "be verbose")

  options.SetParameters("<chrom1,chrom2,...> <input.bw> <output.bw>")
  options.Parse(os.Args)

  if *optHelp {
    options.PrintUsage(os.Stdout)
    os.Exit(0)
  }
  if len(options.Args()) != 3 {
    options.PrintUsage(os.Stderr)
    os.Exit(1)
  }
  chromNames  := strings.Split(options.Args()[0], ",")
  filenameIn  := options.Args()[1]
  filenameOut := options.Args()[2]

  extract(chromNames, filenameIn, filenameOut, *optVerbose)
}
