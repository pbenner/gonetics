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
//import   "math"
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

func WriteResult(config Config, granges GRanges, similarities [][]float64, filenameOut string) {
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
  granges.AddMeta("similarities", similarities)

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

func computeSimilarity(config Config, countsSeq, countsRef KmerCounts) float64 {
  r := 0.0
  for i := 0; i < countsSeq.Len(); i++ {
    x := countsRef.At(i)
    y := countsSeq.At(i)
    r += float64(x*y)
  }
  return r
}

/* -------------------------------------------------------------------------- */

func sequenceSimilarity(config Config, n, m int, filenameRegions, filenameReference, filenameFasta, filenameOut string) {
  pool := threadpool.New(config.Threads, 100*config.Threads)

  kmersCounter, err := NewKmersCounter(n, m, config.Complement, config.Reverse, config.Revcomp, config.MaxAmbiguous, config.Alphabet); if err != nil {
    log.Fatal(err)
  }
  // import reference sequence
  _, reference := ImportData(config, "", filenameReference)
  // make sure the reference contains only a single sequence
  if len(reference) != 1 {
    log.Fatal("reference fasta file should contain a single sequence")
  }
  countsRef := scanSequence(config, kmersCounter, reference[0])
  // freeze k-mer set (i.e. do not consider any new k-mers)
  kmersCounter.Freeze()
  // import target sequences
  granges, sequences := ImportData(config, filenameRegions, filenameFasta)
  // allocate result
  result := make([][]float64, len(sequences))
  // count k-mers on target sequences
  pool.RangeJob(0, len(sequences), func(i int, pool threadpool.ThreadPool, erf func() error) error {
    n := len(reference[0])
    m := len(sequences[i]) - n
    if m <= 0 {
      return nil
    }    
    r := make([]float64, m)
    for j := 0; j < m; j++ {
      countsSeq := scanSequence(config, kmersCounter, sequences[i][j:j+n])
      if countsSeq.Len() != countsRef.Len() {
        panic("internal error")
      }
      r[j] = computeSimilarity(config, countsSeq, countsRef)
    }
    result[i] = r
    return nil
  })
  WriteResult(config, granges, result, filenameOut)
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
  optThreads      := options.    IntLong("threads",       0 ,  1,           "number of threads [default: 1]")
  optComplement   := options.   BoolLong("complement",    0 ,               "consider complement sequences")
  optReverse      := options.   BoolLong("reverse",       0 ,               "consider reverse sequences")
  optRevcomp      := options.   BoolLong("revcomp",       0 ,               "consider reverse complement sequences")
  optVerbose      := options.CounterLong("verbose",      'v',               "verbose level [-v or -vv]")
  optHelp         := options.   BoolLong("help",         'h',               "print help")

  options.SetParameters("<MIN-K-MER-LENGTH> <MAX-K-MER-LENGTH> <REFERENCE.fasta> [<INPUT.fasta> [OUTPUT.table]]")
  options.Parse(os.Args)

  if *optHelp {
    options.PrintUsage(os.Stdout)
    os.Exit(0)
  }
  if len(options.Args()) < 3 || len(options.Args()) > 5 {
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
  config.Binary       = *optBinary
  config.Complement   = *optComplement
  config.Reverse      = *optReverse
  config.Revcomp      = *optRevcomp
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
  filenameReference := options.Args()[2]
  filenameFasta     := ""
  filenameOut       := ""
  if len(options.Args()) >= 4 {
    filenameFasta = options.Args()[3]
  }
  if len(options.Args()) == 5 {
    filenameOut   = options.Args()[4]
  }
  sequenceSimilarity(config, int(n), int(m), *optRegions, filenameReference, filenameFasta, filenameOut)
}
