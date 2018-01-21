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
import   "math"
import   "os"
import   "strconv"

import   "github.com/pborman/getopt"

import . "github.com/pbenner/gonetics"

/* -------------------------------------------------------------------------- */

type Config struct {
  Verbose    int
  BinSize    int
  TrackInit  float64
  BinStat    BinSummaryStatistics
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

func importTrack(config Config, filename string) SimpleTrack {
  track := SimpleTrack{}
  PrintStderr(config, 1, "Importing track `%s'... ", filename)
  if err := track.ImportBigWig(filename, "", config.BinStat, config.BinSize, 0, config.TrackInit); err != nil {
    PrintStderr(config, 1, "failed\n")
    log.Fatal(err)
  } else {
    PrintStderr(config, 1, "done\n")
  }
  return track
}

func exportTrack(config Config, filename string, track SimpleTrack) {
  PrintStderr(config, 1, "Writing track `%s'... ", filename)
  if err := track.ExportBigWig(filename); err != nil {
    PrintStderr(config, 1, "failed\n")
    log.Fatal(err)
  } else {
    PrintStderr(config, 1, "done\n")
  }
}

/* -------------------------------------------------------------------------- */

func bigWigQuantileNormalize(config Config, filenameRef, filenameIn, filenameOut string) {
  trackRef := importTrack(config, filenameRef)
  trackIn  := importTrack(config, filenameIn)

  PrintStderr(config, 1, "Quantile normalizing track... ")
  if err := (GenericMutableTrack{trackIn}).QuantileNormalize(trackRef); err != nil {
    PrintStderr(config, 1, "failed\n")
    log.Fatal(err)
  }
  PrintStderr(config, 1, "done\n")

  exportTrack(config, filenameOut, trackIn)
}

/* -------------------------------------------------------------------------- */

func main() {

  config  := Config{}

  options := getopt.New()

  optBinSize    := options.    IntLong("bin-size",       0 ,      0, "bin size")
  optBinStat    := options. StringLong("bin-summary",    0 , "mean", "bin summary statistic [mean (default), max, min, discrete mean]")
  optTrackInit  := options. StringLong("initial-value",  0 ,     "", "track initial value [default: NaN]")
  optHelp       := options.   BoolLong("help",          'h',         "print help")
  optVerbose    := options.CounterLong("verbose",       'v',         "be verbose")

  options.SetParameters("<REFERENCE.bw> <INPUT.bw> <OUTPUT.bw>")
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
    config.TrackInit = math.NaN()
  }
  config.Verbose    = *optVerbose
  config.BinSize    = *optBinSize
  config.BinStat    = getBinSummaryStatistics(*optBinStat)

  filenameRef := options.Args()[0]
  filenameIn  := options.Args()[1]
  filenameOut := options.Args()[2]

  bigWigQuantileNormalize(config, filenameRef, filenameIn, filenameOut)
}
