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

import . "github.com/pbenner/gonetics"
import . "github.com/pbenner/gonetics/lib/logarithmetic"
import   "github.com/pbenner/gonetics/lib/progress"

import   "github.com/pbenner/threadpool"
import   "github.com/pborman/getopt"

/* -------------------------------------------------------------------------- */

type SessionConfig struct {
  Columns int
  Status  bool
  Summary string
  Threads int
  Verbose int
}

/* -------------------------------------------------------------------------- */

func PrintStderr(config SessionConfig, level int, format string, args ...interface{}) {
  if config.Verbose >= level {
    fmt.Fprintf(os.Stderr, format, args...)
  }
}

/* -------------------------------------------------------------------------- */

func importPWMList(config SessionConfig, filenames []string) []PWM {
  result := []PWM{}
	for _, filename := range filenames {
    p := PWM{}
    PrintStderr(config, 1, "Reading PWM `%s'... ", filename)
    if err := p.ImportMatrix(filename); err != nil {
      PrintStderr(config, 1, "failed\n")
      log.Fatal(err)
    }
    PrintStderr(config, 1, "done\n")
    result = append(result, p)
	}
  return result
}

func importBed(config SessionConfig, filename string) GRanges {
  granges := GRanges{}
  if filename == "" {
    if err := granges.ReadBed(os.Stdin, config.Columns); err != nil {
      log.Fatal(err)
    }
  } else {
    PrintStderr(config, 1, "Reading bed file `%s'... ", filename)
    if err := granges.ImportBed(filename, config.Columns); err != nil {
      PrintStderr(config, 1, "failed\n")
      log.Fatal(err)
    }
    PrintStderr(config, 1, "done\n")
  }
  return granges
}

func importFasta(config SessionConfig, filename string) StringSet {
  s := StringSet{}
  PrintStderr(config, 1, "Reading fasta file `%s'... ", filename)
  if err := s.ImportFasta(filename); err != nil {
    PrintStderr(config, 1, "failed\n")
    log.Fatal(err)
  }
  PrintStderr(config, 1, "done\n")
  return s
}

func exportTable(config SessionConfig, granges GRanges, filename string, header, strand, compress bool, args ...interface{}) {
  if filename == "" {
    if err := granges.WriteTable(os.Stdout, header, strand, args...); err != nil {
      log.Fatal(err)
    } else {
      fmt.Fprintf(os.Stdout, "\n")
    }
  } else {
    PrintStderr(config, 1, "Writing table to `%s'... ", filename)
    if err := granges.ExportTable(filename, header, strand, false, args...); err != nil {
      PrintStderr(config, 1, "failed\n")
      log.Fatal(err)
    }
    PrintStderr(config, 1, "done\n")
  }
}

/* -------------------------------------------------------------------------- */

func scanRegion(config SessionConfig, pwmList []PWM, genomicSequence StringSet, r GRangesRow) []float64 {
  counts := make([]float64, len(pwmList))

  sequence, err := genomicSequence.GetSlice(r.Seqname, r.Range)
  // if sequence is nil, it means the fasta file is missing a chromosome
  if sequence == nil {
    log.Fatalf("sequence `%s' not found in fasta file", r.Seqname)
  }
  // if squence is not nil but there is an error then the region is out of bounds,
  // GetSlice() then returns only that part which actually exists
  if err != nil {
    log.Fatal(err.Error())
  }
  // loop over pwm list and scan the ith region
  for j := 0; j < len(pwmList); j++ {
    switch config.Summary {
    case "max":
      tmp1 := pwmList[j].MaxScore(sequence, true)
      // reverse complement
      tmp2 := pwmList[j].MaxScore(sequence, false)
      // take the maximum
      counts[j] = math.Max(tmp1, tmp2)
    case "mean":
      tmp1 := pwmList[j].MeanScore(sequence, true)
      // reverse complement
      tmp2 := pwmList[j].MeanScore(sequence, false)
      // take the maximum
      counts[j] = LogAdd(tmp1, tmp2) - math.Log(2.0)
    default:
      log.Fatalf("invalid summary statistic: %s", config.Summary)
    }
  }

  return counts
}

func scanRegions(config SessionConfig, granges GRanges, pwmList []PWM, genomicSequence StringSet) GRanges {

  pool   := threadpool.New(config.Threads, 100*config.Threads)
  counts := make([][]float64, granges.Length())
  
  if !config.Status {
    PrintStderr(config, 1, "Scanning for PWM hits... ")
  }
  g := pool.NewJobGroup()

  for n, i := granges.Length(), 0; i < n; i++ {
    // make a thread safe copy of i
    j := i
    // add task to the thread pool
    pool.AddJob(g, func(pool threadpool.ThreadPool, erf func() error) error {
      counts[j] = scanRegion(config, pwmList, genomicSequence, granges.Row(j))
      return nil
    })
    if config.Status {
      progress.New(n, 1000).PrintStderr(i+1)
    }
  }
  pool.Wait(g)
  if !config.Status {
    PrintStderr(config, 1, "done\n")
  }
  granges.AddMeta("counts", counts)

  return granges
}

/* -------------------------------------------------------------------------- */

func tfbsScan(config SessionConfig, filenameGRanges, filenameOut, filenameFasta string, filenamesPWM []string) {
  pwmList := importPWMList(config, filenamesPWM)
  granges := importBed(config, filenameGRanges)
  genomicSequence := importFasta(config, filenameFasta)

  granges = scanRegions(config, granges, pwmList, genomicSequence)
  exportTable(config, granges, filenameOut, true, true, false)
}

/* -------------------------------------------------------------------------- */

func main() {

  config  := SessionConfig{}
  options := getopt.New()

  optInput    := options. StringLong("input",         0 , "",    "read regions from file")
  optColumns  := options.    IntLong("input-columns", 0 ,  3,    "number of bed input columns [3 (default), 6, 9]")
  optOutput   := options. StringLong("output",        0 , "",    "write results to file")
  optStatus   := options.   BoolLong("status",        0 ,        "show status bar")
  optSummary  := options. StringLong("summary",       0 , "max", "summary statistics [mean, max]")
  optThreads  := options.    IntLong("threads",       0 ,  1,    "number of threads")

  optVerbose  := options.CounterLong("verbose",    'v',         "verbose level [-v or -vv]")
  optHelp     := options.   BoolLong("help",       'h',         "print help")

  options.SetParameters("<GENOME.fa> <PWM.table ...>")
  options.Parse(os.Args)

  if *optHelp {
    options.PrintUsage(os.Stdout)
    os.Exit(0)
  }
  if len(options.Args()) < 2 {
    options.PrintUsage(os.Stderr)
    os.Exit(1)
  }
  n := len(options.Args())
  filenameFasta   := options.Args()[0]
  filenamesPWM    := options.Args()[1:n]

  config.Columns = *optColumns
  config.Status  = *optStatus
  config.Summary = *optSummary
  config.Threads = *optThreads
  config.Verbose = *optVerbose

  tfbsScan(config, *optInput, *optOutput, filenameFasta, filenamesPWM)
}
