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
import   "strings"

import   "github.com/pborman/getopt"

import . "github.com/pbenner/gonetics"

/* -------------------------------------------------------------------------- */

type Config struct {
  OutFormat  string
  Verbose    int
}

/* -------------------------------------------------------------------------- */

func PrintStderr(config Config, level int, format string, args ...interface{}) {
  if config.Verbose >= level {
    fmt.Fprintf(os.Stderr, format, args...)
  }
}

/* -------------------------------------------------------------------------- */

func importBed3(config Config, filename string) GRanges {
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
    } else {
      PrintStderr(config, 1, "done\n")
    }
  }
  return granges
}

func importFasta(config Config, filename string) StringSet {
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

func exportFasta(config Config, genomicSequences StringSet, filename string) {
  if filename == "" {
    if err := genomicSequences.WriteFasta(os.Stdout); err != nil {
      log.Fatal(err)
    }
  } else {
    PrintStderr(config, 1, "Writing table `%s'... ", filename)
    if err := genomicSequences.ExportFasta(filename, false); err != nil {
      PrintStderr(config, 1, "failed\n")
      log.Fatal(err)
    } else {
      PrintStderr(config, 1, "done\n")
    }
  }
}

func exportTable(config Config, granges GRanges, filename string) {
  if filename == "" {
    if err := granges.WriteTable(os.Stdout, true, false); err != nil {
      log.Fatal(err)
    }
  } else {
    PrintStderr(config, 1, "Writing table `%s'... ", filename)
    if err := granges.ExportTable(filename, true, false, false); err != nil {
      PrintStderr(config, 1, "failed\n")
      log.Fatal(err)
    } else {
      PrintStderr(config, 1, "done\n")
    }
  }
}

/* -------------------------------------------------------------------------- */

func exportAsFasta(config Config, regions GRanges, sequences []string, filenameOutput string) {
  ss := EmptyStringSet()
  for i := 0; i < regions.Length(); i++ {
    name := fmt.Sprintf("%s|%d|%d", regions.Seqnames[i], regions.Ranges[i].From, regions.Ranges[i].To)
    ss[name] = []byte(sequences[i])
  }
  exportFasta(config, ss, filenameOutput)
}

func exportAsTable(config Config, regions GRanges, sequences []string, filenameOutput string) {
  regions.AddMeta("sequences", sequences)
  exportTable(config, regions, filenameOutput)
}

/* -------------------------------------------------------------------------- */

func extract(config Config, filenameFasta, filenameInput, filenameOutput string) {
  genomicSequence := importFasta(config, filenameFasta)
  // import bed file first to check if it exists
  regions   := importBed3(config, filenameInput)
  sequences := make([]string, regions.Length())

  for i := 0; i < regions.Length(); i++ {
    sequence, err := genomicSequence.GetSlice(regions.Seqnames[i], regions.Ranges[i])
    // if sequence is nil, it means the fasta file is missing a chromosome
    if sequence == nil {
      log.Fatalf("sequence `%s' not found in fasta file", regions.Seqnames[i])
    }
    // if squence is not nil but there is an error then the region is out of bounds,
    // GetSlice() then returns only that part which actually exists
    if err != nil {
      log.Fatal(err.Error())
    }
    sequences[i] = string(sequence)
  }

  switch config.OutFormat {
  case "fasta":
    exportAsFasta(config, regions, sequences, filenameOutput)
  case "table":
    exportAsTable(config, regions, sequences, filenameOutput)
  default:
    panic("internal error")
  }
}

/* -------------------------------------------------------------------------- */

func main() {

  config  := Config{}

  options := getopt.New()

  optOutFormat  := options. StringLong("output-format",  0 , "fasta", "output file format [fasta [default], table]")
  optHelp       := options.   BoolLong("help",          'h',          "print help")
  optVerbose    := options.CounterLong("verbose",       'v',          "be verbose")

  options.SetParameters("<SEQUENCES.fa> [INPUT.bed [OUTPUT.(table|fa)]]\n")
  options.Parse(os.Args)

  if *optHelp {
    options.PrintUsage(os.Stdout)
    os.Exit(0)
  }
  if len(options.Args()) != 1 && len(options.Args()) != 2 && len(options.Args()) != 3 {
    options.PrintUsage(os.Stderr)
    os.Exit(1)
  }
  config.Verbose    = *optVerbose
  config.OutFormat  = strings.ToLower(*optOutFormat)
  switch config.OutFormat {
  case "fasta":
  case "table":
  default:
    options.PrintUsage(os.Stderr)
    os.Exit(1)
  }
  filenameFasta  := ""
  filenameInput  := ""
  filenameOutput := ""
  if len(options.Args()) > 0 {
    filenameFasta  = options.Args()[0]
  }
  if len(options.Args()) > 1 {
    filenameInput  = options.Args()[1]
  }
  if len(options.Args()) > 2 {
    filenameOutput = options.Args()[1]
  }
  extract(config, filenameFasta, filenameInput, filenameOutput)
}
