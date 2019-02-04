/* Copyright (C) 2019 Philipp Benner
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

func exportBed3(config Config, granges GRanges, filename string) {
  if filename == "" {
    if err := granges.WriteBed3(os.Stdout); err != nil {
      log.Fatal(err)
    }
  } else {
    PrintStderr(config, 1, "Writing table `%s'... ", filename)
    if err := granges.ExportBed3(filename, false); err != nil {
      PrintStderr(config, 1, "failed\n")
      log.Fatal(err)
    } else {
      PrintStderr(config, 1, "done\n")
    }
  }
}

/* -------------------------------------------------------------------------- */

func fastaUnresolvedRegions(config Config, filenameFasta, filenameOutput string) {
  genomicSequence := importFasta(config, filenameFasta)

  seqnames := []string{}
  from     := []int{}
  to       := []int{}

  for name, sequence := range genomicSequence {
    for i := 0; i < len(sequence); i++ {
      if sequence[i] == 'n' || sequence[i] == 'N' {
        i_from := i
        i_to   := i+1
        for i++; i < len(sequence); i++ {
          if sequence[i] == 'n' || sequence[i] == 'N' {
            i_to = i+1
          } else {
            break
          }
        }
        seqnames = append(seqnames, name)
        from     = append(from,     i_from)
        to       = append(to,       i_to)
      }
    }
  }
  exportBed3(config, NewGRanges(seqnames, from, to, nil), filenameOutput)
}

/* -------------------------------------------------------------------------- */

func main() {

  config  := Config{}

  options := getopt.New()

  optHelp       := options.   BoolLong("help",    'h', "print help")
  optVerbose    := options.CounterLong("verbose", 'v', "be verbose")

  options.SetParameters("<SEQUENCES.fa> [OUTPUT.bed]\n")
  options.Parse(os.Args)

  if *optHelp {
    options.PrintUsage(os.Stdout)
    os.Exit(0)
  }
  if len(options.Args()) != 1 && len(options.Args()) != 2 {
    options.PrintUsage(os.Stderr)
    os.Exit(1)
  }
  config.Verbose    = *optVerbose
  filenameFasta  := options.Args()[0]
  filenameOutput := ""
  if len(options.Args()) == 2 {
    filenameOutput = options.Args()[1]
  }
  fastaUnresolvedRegions(config, filenameFasta, filenameOutput)
}
