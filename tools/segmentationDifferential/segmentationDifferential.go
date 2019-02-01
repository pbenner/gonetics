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
  Verbose int
}

/* i/o
 * -------------------------------------------------------------------------- */

func printStderr(config Config, level int, format string, args ...interface{}) {
  if config.Verbose >= level {
    fmt.Fprintf(os.Stderr, format, args...)
  }
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

func callDifferentialRegions(states string, segmentationFilenames []string) {
  regions := make([]GRanges, len(segmentationFilenames))
  for i, filename := range segmentationFilenames {
    if err := regions[i].ImportBed6(filename); err != nil {
      log.Fatal(err)
    }
    regions[i] = filterRegions(strings.Split(states, ","), regions[i])
  }
  r := GRanges{}
  r  = r.Merge(regions...)
  occursIn := make([][]int, r.Length())
  for i := 0; i < len(regions); i++ {
    queryIdx, _ := FindOverlaps(r, regions[i])
    for j := 0; j < len(queryIdx); j++ {
      qj := queryIdx[j]
      // append j if it was not already appended before
      if n := len(occursIn[qj]); n == 0 || occursIn[qj][n-1] != j {
        occursIn[qj] = append(occursIn[qj], j)
      }
    }
  }
  names := make([]string, r.Length())
  for i := 0; i < len(names); i++ {
    names[i] = states
  }
  r.AddMeta("name", names)
  r.AddMeta("occurrence", occursIn)
  r.WriteBed6(os.Stdout)
}

/* -------------------------------------------------------------------------- */

func main() {

  config  := Config{}
  options := getopt.New()

  optVerbose           := options.CounterLong("verbose",   'v',     "verbose level [-v or -vv]")
  optHelp              := options.   BoolLong("help",      'h',     "print help")

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
  callDifferentialRegions(options.Args()[0], options.Args()[1:])
}
