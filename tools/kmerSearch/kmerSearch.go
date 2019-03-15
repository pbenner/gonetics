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
import   "bufio"
import   "log"
import   "math"
import   "io"
import   "os"
import   "strconv"

import   "github.com/pborman/getopt"

import . "github.com/pbenner/gonetics"
import   "github.com/pbenner/threadpool"

/* -------------------------------------------------------------------------- */

type Config struct {
  Complement bool
  Reverse    bool
  Revcomp    bool
  Human      bool
  Header     bool
  Threads    int
  Verbose    int
}

/* i/o
 * -------------------------------------------------------------------------- */

func PrintStderr(config Config, level int, format string, args ...interface{}) {
  if config.Verbose >= level {
    fmt.Fprintf(os.Stderr, format, args...)
  }
}

/* -------------------------------------------------------------------------- */

func ipow(x, k int) int {
  return int(math.Pow(float64(x), float64(k)))
}

/* -------------------------------------------------------------------------- */

func ImportBed3(config Config, filename string) GRanges {
  granges := GRanges{}
  PrintStderr(config, 1, "Reading bed file `%s'... ", filename)
  if err := granges.ImportBed3(filename); err != nil {
    PrintStderr(config, 1, "failed\n")
    log.Fatal(err)
  } else {
    PrintStderr(config, 1, "done\n")
  }
  return granges
}

func ImportFasta(config Config, filename string) OrderedStringSet {
  s := OrderedStringSet{}
  if filename == "" {
    if err := s.ReadFasta(os.Stdin); err != nil {
      log.Fatal(err)
    }
  } else {
    PrintStderr(config, 1, "Reading fasta file `%s'... ", filename)
    if err := s.ImportFasta(filename); err != nil {
      PrintStderr(config, 1, "failed\n")
      log.Fatal(err)
    }
    PrintStderr(config, 1, "done\n")
  }
  return s
}

func WriteResult(config Config, kmersCounter KmersCounter, granges GRanges, filenameOut string) {
  var writer io.Writer

  if filenameOut == "" {
    writer = os.Stdout
  } else {
    f, err := os.Create(filenameOut)
    if err != nil {
      log.Fatal(err)
    }
    buffer := bufio.NewWriter(f)
    writer  = buffer
    defer f.Close()
    defer buffer.Flush()
  }
  // write header
  if config.Header {
    for i := 0; i < kmersCounter.CodeLength(); i++ {
      fmt.Fprintf(writer, "# %s\n", kmersCounter.KmerName(i))
    }
  }
  // convert kmer counts to human readable string
  if config.Human {
    kmers    := granges.GetMeta("kmers").([][]int)
    kmersNew := make([]string, len(kmers))
    for i, _ := range kmers {
      for j, _ := range kmers[i] {
        if len(kmersNew[i]) == 0 {
          kmersNew[i] = fmt.Sprintf("%s=%d", kmersCounter.KmerName(j), kmers[i][j])
        } else {
          kmersNew[i] = fmt.Sprintf("%s,%s=%d", kmersNew[i], kmersCounter.KmerName(j), kmers[i][j])
        }
      }
    }
    granges.AddMeta("kmers", kmersNew)
  }
  if err := granges.WriteTable(writer, true, false); err != nil {
    log.Fatal(err)
  }
}

/* -------------------------------------------------------------------------- */

func ImportData(config Config, filenameRegions, filenameFasta string) (GRanges, [][]byte) {
  ss := ImportFasta(config, filenameFasta)
  if filenameRegions == "" {
    seqnames  := ss.Seqnames
    from      := make([]int,    len(seqnames))
    to        := make([]int,    len(seqnames))
    sequences := make([][]byte, len(seqnames))
    for i := 0; i < len(seqnames); i++ {
      sequences[i] = ss.Sequences[seqnames[i]]
      to       [i] = len(sequences[i])
    }
    return NewGRanges(seqnames, from, to, nil), sequences
  } else {
    regions   := ImportBed3(config, filenameRegions)
    sequences := make([][]byte, regions.Length())
    for i := 0; i < regions.Length(); i++ {
      sequence, err := ss.GetSlice(regions.Seqnames[i], regions.Ranges[i])
      // if sequence is nil, it means the fasta file is missing a chromosome
      if sequence == nil {
        log.Fatalf("sequence `%s' not found in fasta file", regions.Seqnames[i])
      }
      // if squence is not nil but there is an error then the region is out of bounds,
      // GetSlice() then returns only that part which actually exists
      if err != nil {
        log.Fatal(err.Error())
      }
      sequences[i] = sequence
    }
    return regions, sequences
  }
}

/* -------------------------------------------------------------------------- */

func scanSequence(config Config, kmersCounter KmersCounter, sequence []byte) []int {
  result := make([]int, kmersCounter.CodeLength())
  if err := kmersCounter.CountKmers(result, sequence); err != nil {
    log.Fatal(err)
  }
  return result
}

/* -------------------------------------------------------------------------- */

func kmerSearch(config Config, n, m int, filenameRegions, filenameFasta, filenameOut string) {
  pool := threadpool.New(config.Threads, 100*config.Threads)
  jg   := pool.NewJobGroup()
  granges, sequences := ImportData(config, filenameRegions, filenameFasta)

  result := make([][]int, len(sequences))
  kmersCounter, err := NewKmersCounter(n, m, config.Complement, config.Reverse, config.Revcomp); if err != nil {
    log.Fatal(err)
  }

  pool.AddRangeJob(0, len(sequences), jg, func(i int, pool threadpool.ThreadPool, erf func() error) error {
    result[i] = scanSequence(config, kmersCounter, sequences[i])
    return nil
  })

  granges.AddMeta("kmers", result)
  WriteResult(config, kmersCounter, granges, filenameOut)
}

/* -------------------------------------------------------------------------- */

func main() {

  config  := Config{}
  options := getopt.New()

  optRegions    := options. StringLong("regions",    0 , "", "bed with with regions")
  optHeader     := options.   BoolLong("header",     0 ,     "print kmer header")
  optHuman      := options.   BoolLong("human",      0 ,     "print human readable kmer statistics")
  optThreads    := options.    IntLong("threads",    0 ,  1, "number of threads [default: 1]")
  optComplement := options.   BoolLong("complement", 0 ,     "consider complement sequences")
  optReverse    := options.   BoolLong("reverse",    0 ,     "consider reverse sequences")
  optRevcomp    := options.   BoolLong("revcomp",    0 ,     "consider reverse complement sequences")
  optVerbose    := options.CounterLong("verbose",   'v',     "verbose level [-v or -vv]")
  optHelp       := options.   BoolLong("help",      'h',     "print help")

  options.SetParameters("<MIN-KMER-LENGTH> <MAX-KMER-LENGTH> [<INPUT.fasta> [OUTPUT.table]]")
  options.Parse(os.Args)

  if *optHelp {
    options.PrintUsage(os.Stdout)
    os.Exit(0)
  }
  if len(options.Args()) < 2 || len(options.Args()) > 4 {
    options.PrintUsage(os.Stderr)
    os.Exit(1)
  }
  config.Complement = *optComplement
  config.Reverse    = *optReverse
  config.Revcomp    = *optRevcomp
  config.Header     = *optHeader
  config.Human      = *optHuman
  config.Threads    = *optThreads
  config.Verbose    = *optVerbose
  // check required arguments
  n, err := strconv.ParseInt(options.Args()[0], 10, 64); if err != nil {
    options.PrintUsage(os.Stderr)
    os.Exit(1)
  }
  m, err := strconv.ParseInt(options.Args()[1], 10, 64); if err != nil {
    options.PrintUsage(os.Stderr)
    os.Exit(1)
  }
  if n < 1 || m < n {
    options.PrintUsage(os.Stderr)
    os.Exit(1)
  }
  filenameFasta := ""
  filenameOut   := ""
  if len(options.Args()) >= 3 {
    filenameFasta = options.Args()[2]
  }
  if len(options.Args()) == 4 {
    filenameOut   = options.Args()[3]
  }
  kmerSearch(config, int(n), int(m), *optRegions, filenameFasta, filenameOut)
}
