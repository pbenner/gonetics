/* Copyright (C) 2016-2017 Philipp Benner
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

import   "github.com/pborman/getopt"
import . "github.com/pbenner/gonetics"

/* -------------------------------------------------------------------------- */

type Config struct {
  Bins       int
  BinSize    int
  BinStat    BinSummaryStatistics
  Cumulative bool
  Verbose    int
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

func track_histogram(config Config, filename string) {
  track := SimpleTrack{}
  PrintStderr(config, 1, "Importing track `%s'... ", filename)
  if err := track.ImportBigWig(filename, "", config.BinStat, config.BinSize, 0, math.NaN()); err != nil {
    PrintStderr(config, 1, "failed\n")
    log.Fatal(err)
  } else {
    PrintStderr(config, 1, "done\n")
  }
  statistics := GenericTrack{track}.SummaryStatistics()

  var histogram TrackHistogram
  if config.Cumulative {
    histogram = GenericTrack{track}.CumulativeHistogram(statistics.Min, statistics.Max, config.Bins)
  } else {
    histogram = GenericTrack{track}.Histogram(statistics.Min, statistics.Max, config.Bins)
  }

  fmt.Printf("%15s\t%15s\n", "x", "y")
  for i := 0; i < len(histogram.X); i++ {
    fmt.Printf("%15e\t%15f\n", histogram.X[i], histogram.Y[i])
  }
}

/* -------------------------------------------------------------------------- */

func main() {

  config  := Config{}

  options := getopt.New()

  optBins       := options.    IntLong("bins",         'b',    100, "number of histogram bins")
  optBinSize    := options.    IntLong("bin-size",      0 ,      0, "bin size")
  optBinStat    := options. StringLong("bin-summary",   0 , "mean", "bin summary statistic [mean (default), max, min, discrete mean]")
  optCumulative := options.   BoolLong("cumulative",   'c',         "compute cumulative histogram")
  optHelp       := options.   BoolLong("help",         'h',         "print help")
  optVerbose    := options.CounterLong("verbose",      'v',         "verbose level [-v or -vv]")

  options.SetParameters("<INPUT.bw>")
  options.Parse(os.Args)

  if *optHelp {
    options.PrintUsage(os.Stdout)
    os.Exit(0)
  }
  if len(options.Args()) != 1 {
    options.PrintUsage(os.Stderr)
    os.Exit(1)
  }
  config.Bins       = *optBins
  config.BinSize    = *optBinSize
  config.BinStat    = getBinSummaryStatistics(*optBinStat)
  config.Cumulative = *optCumulative
  config.Verbose    = *optVerbose

  filenameIn := options.Args()[0]

  track_histogram(config, filenameIn)
}
