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
import   "strconv"

import   "github.com/pborman/getopt"

import . "github.com/pbenner/gonetics"

/* -------------------------------------------------------------------------- */

func getBinSummaryStatisticse(str string) (BinSummaryStatistics, error) {
  switch str {
  case "mean":
    return BinMean, nil
  case "discrete mean":
    return BinDiscreteMean, nil
  case "min":
    return BinMin, nil
  case "max":
    return BinMax, nil
  }
  return nil, fmt.Errorf("invalid bin summary statistics")
}

func query(filenameIn, chrom string, from, to, binSize, binOverlap int, summary BinSummaryStatistics, verbose bool) {
  f, err := os.Open(filenameIn); if err != nil {
    log.Fatal(err)
  }
  defer f.Close()

  if reader, err := NewBigWigReader(f); err != nil {
    log.Fatal(err)
  } else {
    if r, _, err := reader.QuerySlice(chrom, from, to, summary, binSize, binOverlap, math.NaN()); err != nil {
      log.Fatal(err)
    } else {
      fmt.Println(r)
    }
  }
}

func main() {
  options := getopt.New()
  options.SetProgram(fmt.Sprintf("%s", os.Args[0]))

  optBinOverlap := options.   IntLong("bin-overlap", 'b',  0, "number of overlapping bins when computing the summary")
  optSummary    := options.StringLong("summary",     's', "", "bin summary statistic [mean (default), max, min]")
  optHelp       := options.  BoolLong("help",        'h',     "print help")
  optVerbose    := options.  BoolLong("verbose",     'v',     "be verbose")

  options.SetParameters("<input.bw> <chrom> <from> <to> <binsize>")
  options.Parse(os.Args)

  binSummary := BinMean

  if *optHelp {
    options.PrintUsage(os.Stdout)
    os.Exit(0)
  }
  if len(options.Args()) != 5 {
    options.PrintUsage(os.Stderr)
    os.Exit(1)
  }
  if *optSummary != "" {
    if s, err := getBinSummaryStatisticse(*optSummary); err != nil {
      log.Fatal(err)
    } else {
      binSummary = s
    }
  }
  filenameIn  := options.Args()[0]
  chrom       := options.Args()[1]

  from,    err := strconv.ParseInt(options.Args()[2], 10, 64); if err != nil {
    log.Fatal(err)
  }
  to,      err := strconv.ParseInt(options.Args()[3], 10, 64); if err != nil {
    log.Fatal(err)
  }
  binSize, err := strconv.ParseInt(options.Args()[4], 10, 64); if err != nil {
    log.Fatal(err)
  }
  query(filenameIn, chrom, int(from), int(to), int(binSize), *optBinOverlap, binSummary, *optVerbose)
}
