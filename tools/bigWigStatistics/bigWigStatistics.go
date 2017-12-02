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

import   "github.com/pborman/getopt"

import . "github.com/pbenner/gonetics"

/* -------------------------------------------------------------------------- */

func PrintStderr(verbose int, level int, format string, args ...interface{}) {
  if verbose >= level {
    fmt.Fprintf(os.Stderr, format, args...)
  }
}

/* -------------------------------------------------------------------------- */

func bigWigStatistics(filename string, binSize, verbose int) {
  track := LazyTrackFile{}
  PrintStderr(verbose, 1, "Lazy importing track `%s'... ", filename)
  if err := track.ImportBigWig(filename, "", BinMean, binSize, 0, math.NaN()); err != nil {
    PrintStderr(verbose, 1, "failed\n")
    log.Fatal(err)
  } else {
    defer track.Close()
    PrintStderr(verbose, 1, "done\n")
  }
  s := GenericTrack{track}.SummaryStatistics()

  fmt.Println(s)
}

/* -------------------------------------------------------------------------- */

func main() {
  options := getopt.New()
  options.SetProgram(fmt.Sprintf("%s", os.Args[0]))

  optBinSize    := options.    IntLong("bin-size",     0 ,  0, "bin size")
  optHelp       := options.   BoolLong("help",        'h',     "print help")
  optVerbose    := options.CounterLong("verbose",     'v',     "be verbose")

  options.SetParameters("<input.bw>")
  options.Parse(os.Args)

  if *optHelp {
    options.PrintUsage(os.Stdout)
    os.Exit(0)
  }
  if len(options.Args()) != 1 {
    options.PrintUsage(os.Stderr)
    os.Exit(1)
  }
  filename := options.Args()[0]

  bigWigStatistics(filename, *optBinSize, *optVerbose)
}
