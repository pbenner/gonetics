/* Copyright (C) 2018 Philipp Benner
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

import . "github.com/pbenner/gonetics"

import   "github.com/pborman/getopt"

/* -------------------------------------------------------------------------- */

type SessionConfig struct {
  OutputFormat string
  Verbose      int
}

/* -------------------------------------------------------------------------- */

func PrintStderr(config SessionConfig, level int, format string, args ...interface{}) {
  if config.Verbose >= level {
    fmt.Fprintf(os.Stderr, format, args...)
  }
}

/* -------------------------------------------------------------------------- */

func importBed(config SessionConfig, filename string) GRanges {
  granges := GRanges{}
  if filename == "" {
    if err := granges.ReadBed3(os.Stdin); err != nil {
      log.Fatal(err)
    }
  } else {
    PrintStderr(config, 1, "Reading bed file `%s'... ", filename)
    if err := granges.ImportBed3(filename); err != nil {
      PrintStderr(config, 1, "failed\n")
      log.Fatal(err)
    }
    PrintStderr(config, 1, "done\n")
  }
  return granges
}

func importFasta(config SessionConfig, filename string) StringSet {
  s := StringSet{}
  PrintStderr(config, 1, "Reading fasta file `%s'... ", filename)
  if err := s.ImportFasta(filename); err != nil {
    PrintStderr(config, 1, "failed\n")
    log.Fatal(err)
  }
  PrintStderr(config, 1, "done\n")
  return s
}

/* -------------------------------------------------------------------------- */

func observedOverExpectedCpG(config SessionConfig, filenameFasta, filenameRegions string) {
  granges := importBed(config, filenameRegions)
  genomicSequences := importFasta(config, filenameFasta)

  r := make([]float64, granges.Length())
  
  PrintStderr(config, 1, "Scanning regions... ")
  r, err := granges.ObservedOverExpectedCpG(genomicSequences)
  if err != nil {
    PrintStderr(config, 1, "failed\n")
    log.Fatal(err)
  } else {
    PrintStderr(config, 1, "done\n")
  }
  granges.AddMeta("CpG", r)

  switch config.OutputFormat {
  default: fallthrough
  case "table":
    if err := granges.WriteTable(os.Stdout, true, false); err != nil {
      log.Fatal(err)
    }
    fmt.Fprintf(os.Stdout, "\n")
  case "vector":
    for i := 0; i < len(r); i++ {
      fmt.Println(r[i])
    }
  }
}

/* -------------------------------------------------------------------------- */

func main() {

  config  := SessionConfig{}
  options := getopt.New()

  optOutputFormat := options. StringLong("output-format",  0 , "table", "table (default), vector")

  optVerbose      := options.CounterLong("verbose",       'v',          "verbose level [-v or -vv]")
  optHelp         := options.   BoolLong("help",          'h',          "print help")

  options.SetParameters("<GENOME.fa> [REGIONS.bed]")
  options.Parse(os.Args)

  if *optHelp {
    options.PrintUsage(os.Stdout)
    os.Exit(0)
  }
  if len(options.Args()) != 1 && len(options.Args()) != 2 {
    options.PrintUsage(os.Stderr)
    os.Exit(1)
  }
  filenameFasta   := options.Args()[0]
  filenameRegions := ""

  if len(options.Args()) == 2 {
    filenameRegions = options.Args()[1]
  }
  config.OutputFormat = *optOutputFormat
  config.Verbose      = *optVerbose

  observedOverExpectedCpG(config, filenameFasta, filenameRegions)
}
