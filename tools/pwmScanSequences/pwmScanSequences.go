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

import   "github.com/pborman/getopt"

import . "github.com/pbenner/gonetics"
import . "github.com/pbenner/threadpool"

/* -------------------------------------------------------------------------- */

type Config struct {
  BinSize int
  BinStat BinSummaryStatistics
  Threads int
  Verbose int
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

/* i/o
 * -------------------------------------------------------------------------- */

func PrintStderr(config Config, level int, format string, args ...interface{}) {
  if config.Verbose >= level {
    fmt.Fprintf(os.Stderr, format, args...)
  }
}

/* -------------------------------------------------------------------------- */

func ImportPWM(config Config, filename string) PWM {
  p := PWM{}
  PrintStderr(config, 1, "Reading PWM `%s'... ", filename)
  if err := p.ImportMatrix(filename); err != nil {
    PrintStderr(config, 1, "failed\n")
    log.Fatal(err)
  }
  PrintStderr(config, 1, "done\n")
  return p
}

func ImportFasta(config Config, filename string) StringSet {
  s := StringSet{}
  PrintStderr(config, 1, "Reading fasta file `%s'... ", filename)
  if err := s.ImportFasta(filename); err != nil {
    PrintStderr(config, 1, "failed\n")
    log.Fatal(err)
  }
  PrintStderr(config, 1, "done\n")
  return s
}

/* -------------------------------------------------------------------------- */

func genomeFromStringSet(ss StringSet) Genome {
  seqnames := []string{}
  lengths  := []int{}
  for k, v := range ss {
    seqnames = append(seqnames, k)
    lengths  = append(lengths,  len(v))
  }
  return NewGenome(seqnames, lengths)
}

/* -------------------------------------------------------------------------- */

func pwmScanSequence(config Config, pwm PWM, sequence []byte, r TrackMutableSequence, pool ThreadPool) error {

  jobGroup := pool.NewJobGroup()

  if err := pool.AddRangeJob(0, r.NBins(), jobGroup, func(i int, pool ThreadPool, erf func() error) error {
    j_from := (i+0)*config.BinSize
    j_to   := (i+1)*config.BinSize
    // trim range
    if j_to+pwm.Length() > len(sequence) {
      j_to = len(sequence)-pwm.Length()
    }
    s := BbiSummaryStatistics{}
    // scan subsequence
    for j := j_from; j < j_to; j++ {
      v1, err := pwm.Score(sequence[j:j+pwm.Length()], false); if err != nil {
        panic(err)
      }
      v2, err := pwm.Score(sequence[j:j+pwm.Length()], true); if err != nil {
        panic(err)
      }
      s.AddValue(v1)
      s.AddValue(v2)
    }
    r.SetBin(i, config.BinStat(s.Sum, s.SumSquares, s.Min, s.Max, s.Valid))
    return nil
  }); err != nil {
    return err
  }
  if err := pool.Wait(jobGroup); err != nil {
    return err
  }
  return nil
}

/* -------------------------------------------------------------------------- */

func pwmScanSequences(config Config, filenamePWM, filenameFasta, filenameOut string) {
  pool:= NewThreadPool(config.Threads, 100*config.Threads)
  pwm := ImportPWM(config, filenamePWM)
  genomicSequences := ImportFasta(config, filenameFasta)
  genome := genomeFromStringSet(genomicSequences)
  result := AllocSimpleTrack("", genome, config.BinSize)

  for name, seq := range genomicSequences {
    r, err := result.GetMutableSequence(name)
    if err != nil {
      log.Fatal(err)
    }
    PrintStderr(config, 1, "Scanning sequence `%s'... ", name)
    if err := pwmScanSequence(config, pwm, seq, r, pool); err != nil {
      PrintStderr(config, 1, "failed\n")
      log.Fatal(err)
    } else {
      PrintStderr(config, 1, "done\n")
    }
  }
}

/* -------------------------------------------------------------------------- */

func main() {

  config  := Config{}
  options := getopt.New()

  optBinStat := options. StringLong("bin-summary", 0 , "mean", "bin summary statistic [mean (default), max, min, discrete mean]")
  optBinSize := options.    IntLong("bin-size",    0 ,     10, "track bin size [default: 10]")
  optThreads := options.    IntLong("threads",     0 ,      1, "number of threads [default: 1]")
  optVerbose := options.CounterLong("verbose",    'v',         "verbose level [-v or -vv]")
  optHelp    := options.   BoolLong("help",       'h',         "print help")

  options.SetParameters("<PWM.table> <INPUT.fa> <OUTPUT.bw>")
  options.Parse(os.Args)

  if *optHelp {
    options.PrintUsage(os.Stdout)
    os.Exit(0)
  }
  if len(options.Args()) != 3 {
    options.PrintUsage(os.Stderr)
    os.Exit(1)
  }
  config.BinSize = *optBinSize
  config.BinStat = getBinSummaryStatistics(*optBinStat)
  config.Threads = *optThreads
  config.Verbose = *optVerbose
  // check required arguments
  filenamePWM   := options.Args()[0]
  filenameFasta := options.Args()[1]
  filenameOut   := options.Args()[2]

  pwmScanSequences(config, filenamePWM, filenameFasta, filenameOut)
}
