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

package main

/* -------------------------------------------------------------------------- */

import   "fmt"
import   "log"
import   "os"
import   "strconv"

import   "github.com/pborman/getopt"

import . "github.com/pbenner/gonetics"

/* -------------------------------------------------------------------------- */

type Config struct {
  Exclude string
  Verbose int
}

/* -------------------------------------------------------------------------- */

func PrintStderr(config Config, level int, format string, args ...interface{}) {
  if config.Verbose >= level {
    fmt.Fprintf(os.Stderr, format, args...)
  }
}

/* -------------------------------------------------------------------------- */

func importGenome(config Config, filename string) Genome {
  genome := Genome{}
  PrintStderr(config, 1, "Reading genome `%s'... ", filename)
  if err := genome.Import(filename); err != nil {
    PrintStderr(config, 1, "failed\n")
    log.Fatal(err)
  }
  PrintStderr(config, 1, "done\n")

  return genome
}

func importBed3(config Config, filename string) GRanges {
  granges := GRanges{}
  PrintStderr(config, 1, "Reading bed file `%s'... ", filename)
  if err := granges.ImportBed3(filename); err != nil {
    PrintStderr(config, 1, "failed\n")
    log.Fatal(err)
  } else {
    PrintStderr(config, 1, "done\n")
  }
  return granges
}

func exportBed3(config Config, granges GRanges, filename string) {
  if filename == "" {
    if err := granges.WriteBed3(os.Stdout); err != nil {
      log.Fatal(err)
    }
  } else {
    PrintStderr(config, 1, "Writing bed file `%s'... ", filename)
    if err := granges.ExportBed3(filename, false); err != nil {
      PrintStderr(config, 1, "failed\n")
      log.Fatal(err)
    }
    PrintStderr(config, 1, "done\n")
  }
}

/* -------------------------------------------------------------------------- */

func drawGRanges(config Config, genome Genome, exclude []GRanges, n, l int) GRanges {
  r := GRanges{}
  // generate random regions until we found n that do not overlap
  // with excluded regions
  for r.Length() < n {
    PrintStderr(config, 2, "Generating %d regions...\n", n - r.Length())
    // generate random regions and remove overlaps with the exclude
    s := RandomGRanges(n - r.Length(), l, genome, false)
    for i := 0; i < len(exclude); i++ {
      queryHits, _ := FindOverlaps(s, exclude[i])
      if len(queryHits) > 0 {
        PrintStderr(config, 2, "Deleting %d regions that overlap with excluded regions:\n", len(queryHits))
        PrintStderr(config, 2, "%v\n", s.Subset(queryHits))
        // remove all regions that overlap with excluded regions
        s = s.Remove(queryHits)
      }
    }
    // append remaining regions to result
    r = r.Append(s)
  }
  return r
}

func draw(config Config, filenameGenome string, filenameOut string, n, l int) {
  exclude := GRanges{}
  if config.Exclude != "" {
    exclude = importBed3(config, config.Exclude)
  }
  genome  := importGenome(config, filenameGenome)
  granges := drawGRanges(config, genome, []GRanges{exclude}, n, l)

  exportBed3(config, granges, filenameOut)
}

/* -------------------------------------------------------------------------- */

func main() {

  var n, l int

  config  := Config{}

  options := getopt.New()

  optExclude := options. StringLong("exclude", 'e', "", "comma separated list of bed files containing regions to exclude")
  optVerbose := options.CounterLong("verbose", 'v',     "verbose level [-v or -vv]")
  optHelp    := options.   BoolLong("help",    'h',     "print help")

  options.SetParameters("<GENOME> <N> <LENGTH> [<OUTPUT.bed>]")
  options.Parse(os.Args)

  if *optHelp {
    options.PrintUsage(os.Stdout)
    os.Exit(0)
  }
  if len(options.Args()) != 3 && len(options.Args()) != 4 {
    options.PrintUsage(os.Stderr)
    os.Exit(1)
  }
  config.Exclude = *optExclude
  config.Verbose = *optVerbose

  if t, err := strconv.ParseInt(options.Args()[1], 10, 64); err != nil {
    log.Fatal(err)
  } else {
    n = int(t)
  }
  if t, err := strconv.ParseInt(options.Args()[2], 10, 64); err != nil {
    log.Fatal(err)
  } else {
    l = int(t)
  }

  filenameGenome := options.Args()[0]
  filenameOut    := ""

  if len(options.Args()) == 4 {
    filenameOut = options.Args()[3]
  }

  draw(config, filenameGenome, filenameOut, n, l)
}
