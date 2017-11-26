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
  Verbose   int
  BinSize   int
  BinOver   int
  TrackInit float64
  BinStat   BinSummaryStatistics
}

/* -------------------------------------------------------------------------- */

func PrintStderr(config Config, level int, format string, args ...interface{}) {
  if config.Verbose >= level {
    fmt.Fprintf(os.Stderr, format, args...)
  }
}

/* -------------------------------------------------------------------------- */

func getBinSummaryStatistics(str string) BinSummaryStatistics {
  switch str {
  case "mean":
    return BinMean
  case "discrete mean":
    return BinDiscreteMean
  case "min":
    return BinMin
  case "max":
    return BinMax
  }
  log.Fatal("invalid bin summary statistics: %s", str)
  return nil
}

/* -------------------------------------------------------------------------- */

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

func exportTable(config Config, granges GRanges, filename string) {
  PrintStderr(config, 1, "Writing table `%s'... ", filename)
  if err := granges.ExportTable(filename, true, true, false); err != nil {
    PrintStderr(config, 1, "failed\n")
    log.Fatal(err)
  } else {
    PrintStderr(config, 1, "done\n")
  }
}

/* -------------------------------------------------------------------------- */

func extract(config Config, filenameBw, filenameBed, filenameOut string) {

  // import bed file first to check if it exists
  granges := importBed3(config, filenameBed)
  // if this is a BigWig file, we do not have to import
  // the full track
  if f, err := os.Open(filenameBw); err != nil {
    log.Fatal(err)
  } else {
    if reader, err := NewBigWigReader(f); err != nil {
      log.Fatal(err)
    } else {
      granges.ImportBigWig(reader, "counts", config.BinStat, config.BinSize, config.BinOver, config.TrackInit, false)
    }
    f.Close()
  }
  exportTable(config, granges, filenameOut)
}

/* -------------------------------------------------------------------------- */

func main() {

  config  := Config{}

  options := getopt.New()

  optBinSize   := options.    IntLong("bin-size",     0 ,      0, "bin size")
  optBinStat   := options. StringLong("bin-summary",  0 , "mean", "bin summary statistic [mean (default), max, min, discrete mean]")
  optBinOver   := options.    IntLong("bin-overlap",  0 ,      0, "number of overlapping bins when computing the summary")
  optTrackInit := options. StringLong("initial-value",  0 , "", "track initial value [default: 0]")
  optHelp      := options.   BoolLong("help",        'h',         "print help")
  optVerbose   := options.CounterLong("verbose",     'v',         "be verbose")

  options.SetParameters("<INPUT.bw> <INPUT.bed> <OUTPUT.table>")
  options.Parse(os.Args)

  if *optHelp {
    options.PrintUsage(os.Stdout)
    os.Exit(0)
  }
  if len(options.Args()) != 3 {
    options.PrintUsage(os.Stderr)
    os.Exit(1)
  }
  if *optTrackInit != "" {
    v, err := strconv.ParseFloat(*optTrackInit, 64)
    if err != nil {
      log.Fatalf("parsing initial value failed: %v", err)
    }
    config.TrackInit = v
  } else {
    config.TrackInit = 0.0
  }
  config.Verbose = *optVerbose
  config.BinSize = *optBinSize
  config.BinOver = *optBinOver
  config.BinStat = getBinSummaryStatistics(*optBinStat)

  filenameBw  := options.Args()[0]
  filenameBed := options.Args()[1]
  filenameOut := options.Args()[2]

  extract(config, filenameBw, filenameBed, filenameOut)
}
