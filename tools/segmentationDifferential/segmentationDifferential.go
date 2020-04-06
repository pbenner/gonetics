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
import   "math"
import   "os"
import   "strconv"
import   "strings"

import   "github.com/pborman/getopt"

import . "github.com/pbenner/gonetics"

/* -------------------------------------------------------------------------- */

type Config struct {
  FilterMaxSize int
  RegionSize    int
  Summary       BinSummaryStatistics
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

func filterBw(config Config, filename string, thr float64, s GRanges) GRanges {
  printStderr(config, 1, "Filtering regions with data from `%s'... ", filename)
  f, err := OpenBigWigFile(filename)
  if err != nil {
    printStderr(config, 1, "failed\n")
    log.Fatal(err)
  }
  defer f.Close()

  bwr, err := NewBigWigReader(f)
  if err != nil {
    printStderr(config, 1, "failed\n")
    log.Fatal(err)
  }

  del := []int{}

  for i := 0; i < s.Length(); i++ {
    seqname := s.Seqnames[i]
    from    := s.Ranges[i].From
    to      := s.Ranges[i].To
    if s, _, err := bwr.QuerySlice(seqname, from, to, config.Summary, to-from, 0, math.NaN()); err != nil {
      printStderr(config, 1, "failed\n")
      log.Fatal(err)
    } else {
      if len(s) != 1 {
        panic("internal error")
      }
      if s[0] != math.NaN() && s[0] < thr {
        del = append(del, i)
      }
    }
  }
  printStderr(config, 1, "done\n")
  return s.Remove(del)
}

func filterOutBw(config Config, filename string, thr float64, s GRanges) GRanges {
  printStderr(config, 1, "Filtering regions with data from `%s'... ", filename)
  f, err := OpenBigWigFile(filename)
  if err != nil {
    printStderr(config, 1, "failed\n")
    log.Fatal(err)
  }
  defer f.Close()

  bwr, err := NewBigWigReader(f)
  if err != nil {
    printStderr(config, 1, "failed\n")
    log.Fatal(err)
  }
  del := []int{}

  for i := 0; i < s.Length(); i++ {
    seqname := s.Seqnames[i]
    from    := s.Ranges[i].From
    to      := s.Ranges[i].To
    if s, _, err := bwr.QuerySlice(seqname, from, to, config.Summary, to-from, 0, math.NaN()); err != nil {
      printStderr(config, 1, "failed\n")
      log.Fatal(err)
    } else {
      if len(s) != 1 {
        panic("internal error")
      }
      if s[0] != math.NaN() && s[0] > thr {
        del = append(del, i)
      }
    }
  }
  printStderr(config, 1, "done\n")
  return s.Remove(del)
}

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

func computeScores(config Config, r GRanges, labels [][]int, bwFiles []string) []float64 {
  readers := make([]*BigWigReader, len(bwFiles))
  // open all bigWig files
  for i, _ := range bwFiles {
    printStderr(config, 1, "Importing scores from `%s'... ", bwFiles[i])
    f, err := OpenBigWigFile(bwFiles[i])
    if err != nil {
      printStderr(config, 1, "failed\n")
      log.Fatal(err)
    }
    defer f.Close()

    bwr, err := NewBigWigReader(f)
    if err != nil {
      printStderr(config, 1, "failed\n")
      log.Fatal(err)
    } else {
      printStderr(config, 1, "done\n")
    }
    readers[i] = bwr
  }
  scores_pos := make([]float64, r.Length())
  scores_neg := make([]float64, r.Length())
  for i := 0; i < r.Length(); i++ {
    overlaps := make([]bool, len(bwFiles))
    n        := 0
    for _, j := range labels[i] {
      overlaps[j] = true
      n          += 1
    }
    for j := 0; j < len(bwFiles); j++ {
      if s, _, err := readers[j].QuerySlice(r.Seqnames[i], r.Ranges[i].From, r.Ranges[i].To, BinMax, r.Ranges[i].To-r.Ranges[i].From, 0, 0.0); err != nil {
        log.Fatal(err)
      } else {
        if overlaps[j] {
          scores_pos[i] += s[0]
        } else {
          scores_neg[i] += s[0]
        }
      }
    }
    switch n {
    case 0:
      scores_pos[i] = math.Inf(-1)
    case len(bwFiles):
      scores_pos[i] = math.Inf(-1)
    default:
      scores_pos[i] = math.Log(scores_pos[i]/float64(n)) - math.Log(scores_neg[i]/float64(len(bwFiles)-n))
    }
  }
  return scores_pos
}

/* -------------------------------------------------------------------------- */

func callDifferentialRegions(config Config, states string, segmentationFilenames, bwFiles1, bwFiles2, bwFiles3 []string, bwThr1, bwThr2 []float64) {
  regions := make([]GRanges, len(segmentationFilenames))
  for i, filename := range segmentationFilenames {
    regions[i] = importBed6(config, filename)
    regions[i] = filterRegions(strings.Split(states, ","), regions[i])
    if len(bwFiles1) > 0 {
      regions[i] = filterBw(config, bwFiles1[i], bwThr1[i], regions[i])
    }
    if len(bwFiles2) > 0 {
      regions[i] = filterOutBw(config, bwFiles2[i], bwThr2[i], regions[i])
    }
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
  if len(bwFiles3) > 0 {
    r.AddMeta("score", computeScores(config, r, labels, bwFiles3))
    r, _ = r.Sort("score", true)
  }
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

  optFilterBw      := options. StringLong("filter-bigwig",     0 , "", "remove regions that have a value smaller than the given " +
                                                                       "threshold [format: file1.bw:threshold1,file2.bw:threshold2,...]")
  optFilterOutBw   := options. StringLong("filter-out-bigwig", 0 , "", "remove regions that have a value larger than the given " +
                                                                       "threshold [format: file1.bw:threshold1,file2.bw:threshold2,...]")
  optFilterMaxSize := options.    IntLong("filter-max-size",   0 ,  0, "filter out regions that are longer than the given threshold")
  optScores        := options. StringLong("scores",            0 , "", "bigWig files with scores")
  optRegionSize    := options.    IntLong("region-size",       0 ,  0, "if not zero, all regions are resized to the given length")
  optVerbose       := options.CounterLong("verbose",          'v',     "verbose level [-v or -vv]")
  optHelp          := options.   BoolLong("help",             'h',     "print help")

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
  bwFiles1 := []string {}
  bwFiles2 := []string {}
  bwFiles3 := []string {}
  bwThr1   := []float64{}
  bwThr2   := []float64{}

  if len(*optFilterBw) > 0 {
    for _, str := range strings.Split(*optFilterBw, ",") {
      s := strings.Split(str, ":")
      if len(s) != 2 {
        log.Fatalf("invalid optional argument `%s'", str)
      }
      v, err := strconv.ParseFloat(s[1], 64)
      if err != nil {
        log.Fatal(err)
      }
      bwFiles1 = append(bwFiles1, s[0])
      bwThr1   = append(bwThr1  , v)
    }
  }
  if len(*optFilterOutBw) > 0 {
    for _, str := range strings.Split(*optFilterOutBw, ",") {
      s := strings.Split(str, ":")
      if len(s) != 2 {
        log.Fatalf("invalid optional argument `%s'", str)
      }
      v, err := strconv.ParseFloat(s[1], 64)
      if err != nil {
        log.Fatal(err)
      }
      bwFiles2 = append(bwFiles2, s[0])
      bwThr2   = append(bwThr2  , v)
    }
  }
  if len(*optScores) > 0 {
    for _, str := range strings.Split(*optScores, ",") {
      bwFiles3 = append(bwFiles3, str)
    }
  }
  if len(bwFiles1) != 0 && len(bwFiles1) != len(options.Args())-1 {
    log.Fatal("optional argument `--filter-biwgwig' has invalid number of fields")
  }
  if len(bwFiles2) != 0 && len(bwFiles2) != len(options.Args())-1 {
    log.Fatal("optional argument `--filter-out-biwgwig' has invalid number of fields")
  }
  if len(bwFiles3) != 0 && len(bwFiles3) != len(options.Args())-1 {
    log.Fatal("optional argument `--scores' has invalid number of fields")
  }

  config.FilterMaxSize = *optFilterMaxSize
  config.RegionSize    = *optRegionSize
  config.Summary       = BinSummaryStatisticsFromString("mean")

  callDifferentialRegions(config, options.Args()[0], options.Args()[1:], bwFiles1, bwFiles2, bwFiles3, bwThr1, bwThr2)
}
