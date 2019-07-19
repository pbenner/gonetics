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
import   "io"
import   "os"
import   "strconv"
import   "strings"

import   "github.com/pborman/getopt"

import . "github.com/pbenner/gonetics"
import   "github.com/pbenner/threadpool"

/* -------------------------------------------------------------------------- */

type Config struct {
  Alphabet       ComplementableAlphabet
  Binary         bool
  MaxAmbiguous []int
  Complement     bool
  Reverse        bool
  Revcomp        bool
  Human          bool
  Sparse         bool
  Header         bool
  Threads        int
  Verbose        int
}

/* i/o
 * -------------------------------------------------------------------------- */

func PrintStderr(config Config, level int, format string, args ...interface{}) {
  if config.Verbose >= level {
    fmt.Fprintf(os.Stderr, format, args...)
  }
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

func WriteResult(config Config, granges GRanges, countsList KmerCountsList, filenameOut string) {
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
    for i := 0; i < countsList.Kmers.Len(); i++ {
      fmt.Fprintf(writer, "# %s\n", countsList.Kmers[i].Name)
    }
  }
  // convert kmer counts to human readable string
  if config.Human {
    kmersMeta := make([]string, countsList.Len())
    for i := 0; i < countsList.Len(); i++ {
      counts := countsList.At(i)
      for it := counts.Iterate(); it.Ok(); it.Next() {
        if len(kmersMeta[i]) == 0 {
          kmersMeta[i] = fmt.Sprintf("%s=%d", it.GetName(), it.GetCount())
        } else {
          kmersMeta[i] = fmt.Sprintf("%s,%s=%d", kmersMeta[i], it.GetName(), it.GetCount())
        }
      }
    }
    granges.AddMeta("k-mers", kmersMeta)
  } else
  if config.Sparse {
    kmersMeta := make([]string, countsList.Len())
    for i := 0; i < countsList.Len(); i++ {
      counts := countsList.At(i)
      for it := counts.Iterate(); it.Ok(); it.Next() {
        if it.GetCount() == 0 {
          continue
        }
        if len(kmersMeta[i]) == 0 {
          kmersMeta[i] = fmt.Sprintf("%s=%d", it.GetName(), it.GetCount())
        } else {
          kmersMeta[i] = fmt.Sprintf("%s,%s=%d", kmersMeta[i], it.GetName(), it.GetCount())
        }
      }
    }
    granges.AddMeta("k-mers", kmersMeta)
  } else {
    kmersMeta := make([][]int, countsList.Len())
    for i := 0; i < countsList.Len(); i++ {
      counts := countsList.At(i)
      kmersMeta[i] = make([]int, counts.Len())
      for j, it := 0, counts.Iterate(); it.Ok(); it.Next() {
        kmersMeta[i][j] = it.GetCount()
        j++
      }
    }
    granges.AddMeta("k-mers", kmersMeta)
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

func scanSequence(config Config, kmersCounter *KmersCounter, sequence []byte) KmerCounts {
  if config.Binary {
    return kmersCounter.IdentifyKmers(sequence)
  } else {
    return kmersCounter.CountKmers(sequence)
  }
}

/* -------------------------------------------------------------------------- */

func countKmers(config Config, n, m int, filenameRegions, filenameFasta, filenameOut string) {
  pool := threadpool.New(config.Threads, 100*config.Threads)
  granges, sequences := ImportData(config, filenameRegions, filenameFasta)

  // create a KmersCounter for each thread
  kmersCounter := make([]*KmersCounter, pool.NumberOfThreads())
  for i := 0; i < pool.NumberOfThreads(); i++ {
    if r, err := NewKmersCounter(n, m, config.Complement, config.Reverse, config.Revcomp, config.MaxAmbiguous, config.Alphabet); err != nil {
      log.Fatal(err)
    } else {
      kmersCounter[i] = r
    }
  }
  counts := make([]KmerCounts, len(sequences))

  pool.RangeJob(0, len(sequences), func(i int, pool threadpool.ThreadPool, erf func() error) error {
    counts[i] = scanSequence(config, kmersCounter[pool.GetThreadId()], sequences[i])
    return nil
  })
  countsList := NewKmerCountsList(counts...)
  WriteResult(config, granges, countsList, filenameOut)
}

/* -------------------------------------------------------------------------- */

func main() {
  log.SetFlags(0)

  config  := Config{}
  options := getopt.New()

  optAlphabet     := options. StringLong("alphabet",      0 , "nucleotide", "nucleotide, gapped-nucleotide, or ambiguous-nucleotide")
  optBinary       := options.   BoolLong("binary",        0 ,               "count matrix only shows the presence or absence of k-mers")
  optMaxAmbiguous := options. StringLong("max-ambiguous", 0 , "-1",         "maxum number of ambiguous positions (either a scalar to set a global maximum or a comma separated list of length [MAX-K-MER-LENGTH]-[MIN-K-MER-LENGTH]+1)")
  optRegions      := options. StringLong("regions",       0 , "",           "bed with with regions")
  optHeader       := options.   BoolLong("header",        0 ,               "print k-mer header")
  optHuman        := options.   BoolLong("human",         0 ,               "print human readable k-mer statistics")
  optSparse       := options.   BoolLong("sparse",        0 ,               "print sparse k-mer statistics")
  optThreads      := options.    IntLong("threads",       0 ,  1,           "number of threads [default: 1]")
  optComplement   := options.   BoolLong("complement",    0 ,               "consider complement sequences")
  optReverse      := options.   BoolLong("reverse",       0 ,               "consider reverse sequences")
  optRevcomp      := options.   BoolLong("revcomp",       0 ,               "consider reverse complement sequences")
  optVerbose      := options.CounterLong("verbose",      'v',               "verbose level [-v or -vv]")
  optHelp         := options.   BoolLong("help",         'h',               "print help")

  options.SetParameters("<MIN-K-MER-LENGTH> <MAX-K-MER-LENGTH> [<INPUT.fasta> [OUTPUT.table]]")
  options.Parse(os.Args)

  if *optHelp {
    options.PrintUsage(os.Stdout)
    os.Exit(0)
  }
  if len(options.Args()) < 2 || len(options.Args()) > 4 {
    options.PrintUsage(os.Stderr)
    os.Exit(1)
  }
  switch strings.ToLower(*optAlphabet) {
  case "nucleotide"          : config.Alphabet =          NucleotideAlphabet{}
  case "gapped-nucleotide"   : config.Alphabet =    GappedNucleotideAlphabet{}
  case "ambiguous-nucleotide": config.Alphabet = AmbiguousNucleotideAlphabet{}
  default:
    options.PrintUsage(os.Stderr)
    os.Exit(1)
  }
  if *optHuman && *optSparse {
    log.Fatal("options `--human' and `--sparse' are incompatible")
  }
  if *optHeader && *optSparse {
    log.Fatal("options `--header' and `--sparse' are incompatible")
  }
  config.Binary       = *optBinary
  config.Complement   = *optComplement
  config.Reverse      = *optReverse
  config.Revcomp      = *optRevcomp
  config.Header       = *optHeader
  config.Human        = *optHuman
  config.Sparse       = *optSparse
  config.Threads      = *optThreads
  config.Verbose      = *optVerbose
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
  if fields := strings.Split(*optMaxAmbiguous, ","); len(fields) == 1 || len(fields) == int(m-n+1) {
    config.MaxAmbiguous = make([]int, len(fields))
    for i := 0; i < len(fields); i++ {
      if t, err := strconv.ParseInt(fields[i], 10, 64); err != nil {
        options.PrintUsage(os.Stderr)
        os.Exit(1)
      } else {
        config.MaxAmbiguous[i] = int(t)
      }
    }
  } else {
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

  countKmers(config, int(n), int(m), *optRegions, filenameFasta, filenameOut)
}
