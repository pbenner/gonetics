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
import   "math"
import   "log"
import   "os"
import   "sort"
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

type cumDist struct {
  x []float64
  y []int
  n   int
}

func (obj cumDist) Len() int {
  return len(obj.x)
}

func (obj cumDist) Less(i, j int) bool {
  return obj.x[i] < obj.x[j]
}

func (obj cumDist) Swap(i, j int) {
  obj.x[i], obj.x[j] = obj.x[j], obj.x[i]
  obj.y[i], obj.y[j] = obj.y[j], obj.y[i]
}

func newCumDist(m map[float64]int) cumDist {
  x := make([]float64, len(m))
  y := make([]int,     len(m))
  i := 0
  for k, v := range m {
    x[i] = k
    y[i] = v
    i++
  }
  c := cumDist{x, y, 0}
  sort.Sort(c)

  // compute cumulative distribution
  n := 0
  for i := 0; i < len(x); i++ {
    n += y[i]; y[i] = n
  }
  c.n = n

  return c
}

/* -------------------------------------------------------------------------- */

func bigWigQuantileNormalize(config Config, filenameRef, filenameIn, filenameOut string) {
  trackRef := importTrack(config, filenameRef)
  trackIn  := importTrack(config, filenameIn)
  mapRef := make(map[float64]int)
  mapIn  := make(map[float64]int)
  mapTr  := make(map[float64]float64)

  PrintStderr(config, 1, "Quantile normalizing track... ")

  if err := (GenericMutableTrack{}).Map(trackRef, func(seqname string, position int, value float64) float64 {
    if !math.IsNaN(value) {
      mapRef[value] += 1
    }
    return 0.0
  }); err != nil {
    PrintStderr(config, 1, "failed\n")
    log.Fatal(err)
  }
  if err := (GenericMutableTrack{}).Map(trackIn, func(seqname string, position int, value float64) float64 {
    if !math.IsNaN(value) {
      mapIn[value] += 1
    }
    return 0.0
  }); err != nil {
    PrintStderr(config, 1, "failed\n")
    log.Fatal(err)
  }
  distRef := newCumDist(mapRef)
  distIn  := newCumDist(mapIn)

  for i, j := 0, 0; i < len(distRef.x); i++ {
    pRef := float64(distRef.y[i])/float64(distRef.n)
    for ; j < len(distIn.x); j++ {
      pIn := float64(distIn.y[j])/float64(distIn.n)
      if pIn > pRef {
        break
      }
      // map input x_j to reference x_i
      mapTr[distIn.x[j]] = distRef.x[i]
    }
  }

  if err := (GenericMutableTrack{trackIn}).Map(trackIn, func(seqname string, position int, value float64) float64 {
    if math.IsNaN(value) {
      return value
    } else {
      return mapTr[value]
    }
  }); err != nil {
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
  optTrackInit  := options. StringLong("initial-value",  0 ,     "", "track initial value [default: 0]")
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
    config.TrackInit = 0.0
  }
  config.Verbose    = *optVerbose
  config.BinSize    = *optBinSize
  config.BinStat    = getBinSummaryStatistics(*optBinStat)

  filenameRef := options.Args()[0]
  filenameIn  := options.Args()[1]
  filenameOut := options.Args()[2]

  bigWigQuantileNormalize(config, filenameRef, filenameIn, filenameOut)
}
