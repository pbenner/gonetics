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
import   "strings"

import   "github.com/pborman/getopt"

import . "github.com/pbenner/gonetics"

/* -------------------------------------------------------------------------- */

type Config struct {
  FilterMaxSize int
  RegionSize    int
  Verbose       int
}

/* i/o
 * -------------------------------------------------------------------------- */

func printStderr(config Config, level int, format string, args ...interface{}) {
  if config.Verbose >= level {
    fmt.Fprintf(os.Stderr, format, args...)
  }
}

/* -------------------------------------------------------------------------- */

func importBed6(config Config, filename string) GRanges {
  granges := GRanges{}
  printStderr(config, 1, "Reading bed file `%s'... ", filename)
  if err := granges.ImportBed6(filename); err != nil {
    printStderr(config, 1, "failed\n")
    log.Fatal(err)
  } else {
    printStderr(config, 1, "done\n")
  }
  return granges
}

/* -------------------------------------------------------------------------- */

func filterRegions(states []string, s GRanges) GRanges {
  idx  := []int{}
  name := s.GetMetaStr("name")
  for i := 0; i < len(name); i++ {
    for _, state := range states {
      if name[i] == state {
        idx = append(idx, i)
        break
      }
    }
  }
  return s.Subset(idx)
}

/* -------------------------------------------------------------------------- */

func callDifferentialRegions(config Config, states string, segmentationFilenames []string) {
  regions := make([]GRanges, len(segmentationFilenames))
  for i, filename := range segmentationFilenames {
    regions[i] = importBed6(config, filename)
    regions[i] = filterRegions(strings.Split(states, ","), regions[i])
  }
  r := GRanges{}
  r  = r.Merge(regions...)
  labels := make([][]int, r.Length())
  for i := 0; i < len(regions); i++ {
    queryIdx, _ := FindOverlaps(r, regions[i])
    for j := 0; j < len(queryIdx); j++ {
      qj := queryIdx[j]
      // append j if it was not already appended before
      if n := len(labels[qj]); n == 0 || labels[qj][n-1] != i {
        labels[qj] = append(labels[qj], i)
      }
    }
  }
  names := make([]string, r.Length())
  for i := 0; i < len(names); i++ {
    names[i] = states
  }
  // add meta information
  r.AddMeta("name"  , names)
  r.AddMeta("labels", labels)
  // filter regions
  if config.FilterMaxSize != 0 {
    idx := []int{}
    for i := 0; i < r.Length(); i++ {
      if r.Ranges[i].To - r.Ranges[i].From > config.FilterMaxSize {
        idx = append(idx, i)
      }
    }
    r = r.Remove(idx)
  }
  // resize regions
  if config.RegionSize != 0 {
    for i := 0; i < r.Length(); i++ {
      x := (r.Ranges[i].From + r.Ranges[i].To)/2
      r.Ranges[i].From = x - config.RegionSize/2
      r.Ranges[i].To   = r.Ranges[i].From + config.RegionSize
    }
  }
  r.WriteTable(os.Stdout, true, false)
}

/* -------------------------------------------------------------------------- */

func main() {

  config  := Config{}
  options := getopt.New()

  optFilterMaxSize := options.    IntLong("filter-max-size",  0 , 0,  "filter out regions that are longer than the given threshold")
  optRegionSize    := options.    IntLong("region-size",      0 , 0,  "if not zero, all regions are resized to the given length")
  optVerbose       := options.CounterLong("verbose",         'v',     "verbose level [-v or -vv]")
  optHelp          := options.   BoolLong("help",            'h',     "print help")

  options.SetParameters("<STATE_1,STATE_2,...> <SEGMENTATION_1.bed> <SEGMENTATION_2.bed> [SEGMENTATION_3.bed]...")
  options.Parse(os.Args)

  // parse options
  //////////////////////////////////////////////////////////////////////////////
  if *optHelp {
    options.PrintUsage(os.Stdout)
    os.Exit(0)
  }
  if *optVerbose != 0 {
    config.Verbose = *optVerbose
  }
  if len(options.Args()) < 3 {
    options.PrintUsage(os.Stderr)
    os.Exit(1)
  }
  config.FilterMaxSize = *optFilterMaxSize
  config.RegionSize    = *optRegionSize

  callDifferentialRegions(config, options.Args()[0], options.Args()[1:])
}
