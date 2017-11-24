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
import   "plugin"

import   "github.com/pborman/getopt"

import . "github.com/pbenner/gonetics"

/* -------------------------------------------------------------------------- */

type Config struct {
  Verbose int
  BinSize int
  BinOver int
  BinStat BinSummaryStatistics
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

func importTracks(config Config, filenames []string) []Track {
  tracks := make([]Track, len(filenames))
  for i, filename := range filenames {
    t := LazyTrackFile{}
    PrintStderr(config, 1, "Lazy importing track `%s'... ", filename)
    if err := t.ImportBigWig(filename, "", config.BinStat, config.BinSize, config.BinOver, math.NaN()); err != nil {
      PrintStderr(config, 1, "failed\n")
      log.Fatal(err)
    } else {
      PrintStderr(config, 1, "done\n")
    }
    tracks[i] = t
  }
  return tracks
}

/* -------------------------------------------------------------------------- */

func bigWigMap(config Config, filenamePlugin, filenameOutput string, filenameInputs []string) {

  var f func(string, int, ...float64) float64

  p, err := plugin.Open(filenamePlugin)
  if err != nil {
    log.Fatalf("opening plugin `%s' failed: %v", filenamePlugin, err)
  }

  g, err := p.Lookup("f")
  if err != nil {
    log.Fatal(err)
  }
  switch t := g.(type) {
  case func(string, int, ...float64) float64:
    f = t
  default:
    log.Fatal("error while executing plugin: function has invalid type")
  }

  tracks := importTracks(config, filenameInputs)
  trackr := AllocSimpleTrack("", tracks[0].GetGenome(), config.BinSize)

  if err := (GenericMutableTrack{trackr}).MapList(tracks, f); err != nil {
    log.Fatal(err)
  }

  PrintStderr(config, 1, "Writing track `%s'... ", filenameOutput)
  parameters := DefaultBigWigParameters()
  if err := (GenericTrack{trackr}).ExportBigWig(filenameOutput, parameters); err != nil {
    PrintStderr(config, 1, "failed\n")
    log.Fatal(err)
  } else {
    PrintStderr(config, 1, "done\n")
  }
}

/* -------------------------------------------------------------------------- */

func main() {
  options := getopt.New()
  options.SetProgram(fmt.Sprintf("%s", os.Args[0]))

  config := Config{}

  optBinSize := options.    IntLong("bin-size",     0 ,      0, "bin size")
  optBinStat := options. StringLong("bin-summary",  0 , "mean", "bin summary statistic [mean (default), max, min, discrete mean]")
  optBinOver := options.    IntLong("bin-overlap",  0 ,      0, "number of overlapping bins when computing the summary")
  optHelp    := options.   BoolLong("help",        'h',         "print help")
  optVerbose := options.CounterLong("verbose",     'v',         "be verbose")

  options.SetParameters("<plugin.so> <output.bw> <input.bw>...")
  options.Parse(os.Args)

  if *optHelp {
    options.PrintUsage(os.Stdout)
    os.Exit(0)
  }
  if len(options.Args()) < 3 {
    options.PrintUsage(os.Stderr)
    os.Exit(1)
  }
  config.Verbose = *optVerbose
  config.BinSize = *optBinSize
  config.BinOver = *optBinOver
  config.BinStat = getBinSummaryStatistics(*optBinStat)

  filenamePlugin := options.Args()[0]
  filenameOutput := options.Args()[1]
  filenameInputs := options.Args()[2:]

  bigWigMap(config, filenamePlugin, filenameOutput, filenameInputs)
}
