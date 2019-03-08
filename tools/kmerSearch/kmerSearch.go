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
  Human    bool
  Header   bool
  Threads  int
  Verbose  int
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

func WriteResult(config Config, kmersMeta KmersMeta, granges GRanges, filenameOut string) {
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
    for _, name := range kmersMeta.names {
      fmt.Fprintf(writer, "# %s\n", name)
    }
  }
  // convert kmer counts to human readable string
  if config.Human {
    kmers    := granges.GetMeta("kmers").([][]int)
    kmersNew := make([]string, len(kmers))
    for i, _ := range kmers {
      for j, _ := range kmers[i] {
        if len(kmersNew[i]) == 0 {
          kmersNew[i] = fmt.Sprintf("%s=%d", kmersMeta.names[j], kmers[i][j])
        } else {
          kmersNew[i] = fmt.Sprintf("%s,%s=%d", kmersNew[i], kmersMeta.names[j], kmers[i][j])
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

type KmersMeta struct {
  n, m, N     int
  indices [][]int
  names     []string
  p         []int
  al          NucleotideAlphabet
}

func NewKmersMeta(n, m int) KmersMeta {
  r := KmersMeta{n: n, m: m}
  p := make([]int, m+1)
  for k := 0; k <= m; k++ {
    p[k] = ipow(r.al.Length(), k)
  }
  names   := []string{}
  indices := make([][]int, m-n+1)
  idx     := 0
  for k := n; k <= m; k++ {
    kn := ipow(r.al.Length(), k)
    indices[k-n] = make([]int, kn)
    c1 := make([]byte, k)
    c2 := make([]byte, k)
    for i := 0; i < kn; i++ {
      // convert index to sequence
      for j, ix := 0, i; j < k; j++ {
        c1[k-j-1] = byte(ix % r.al.Length())
        ix        = ix / r.al.Length()
        if x, err := r.al.ComplementCoded(c1[k-j-1]); err != nil {
          log.Fatal(err)
        } else {
          c2[j] = x
        }
      }
      q := 0
      for j := 0; j < k; j++ {
        q += int(c2[k-j-1]) * p[j]
      }
      if i <= q {
        for j := 0; j < k; j++ {
          if x, err := r.al.Decode(c1[j]); err != nil {
            log.Fatal(err)
          } else {
            c1[j] = x
          }
          if x, err := r.al.Decode(c2[j]); err != nil {
            log.Fatal(err)
          } else {
            c2[j] = x
          }
        }
        names = append(names, fmt.Sprintf("%s|%s", string(c1), string(c2)))
        indices[k-n][i] = idx; idx += 1
      } else {
        indices[k-n][i] = -1
      }
    }
  }
  r.N       = idx
  r.p       = p
  r.names   = names
  r.indices = indices
  return r
}

func (obj KmersMeta) Index(k, i int) int {
  return obj.indices[k-obj.n][i]
}

/* -------------------------------------------------------------------------- */

func scanSequence(config Config, kmersMeta KmersMeta, sequence []byte) []int {
  c1 := make([]int, len(sequence))
  c2 := make([]int, len(sequence))
  for i := 0; i < len(sequence); i++ {
    if sequence[i] == 'n' || sequence[i] == 'N' {
      c1[i] = -1
      c2[i] = -1
      continue
    }
    if r, err := kmersMeta.al.Code(sequence[i]); err != nil {
      log.Fatal(err)
    } else {
      c1[i] = int(r)
    }
    if r, err := kmersMeta.al.ComplementCoded(byte(c1[i])); err != nil {
      log.Fatal(err)
    } else {
      c2[i] = int(r)
    }
  }
  // allocate results vector
  result := make([]int, kmersMeta.N)
  // loop over sequence
  for i := 0; i < len(sequence); i++ {
    // loop over all k-mers
kLoop:
    for k := kmersMeta.n; k <= kmersMeta.m && i+k-1 < len(sequence); k++ {
      // eval k-mer
      s := 0
      r := 0
      for j := 0; j < k; j++ {
        if c1[i+j] == -1 {
          break kLoop
        }
        s += c1[i+j] * kmersMeta.p[k-j-1]
        r += c2[i+j] * kmersMeta.p[j]
      }
      if s < r {
        result[kmersMeta.Index(k, s)] += 1
      } else {
        result[kmersMeta.Index(k, r)] += 1
      }
    }
  }
  return result
}

/* -------------------------------------------------------------------------- */

func kmerSearch(config Config, n, m int, filenameRegions, filenameFasta, filenameOut string) {
  pool := threadpool.New(config.Threads, 100*config.Threads)
  jg   := pool.NewJobGroup()
  granges, sequences := ImportData(config, filenameRegions, filenameFasta)

  result    := make([][]int, len(sequences))
  kmersMeta := NewKmersMeta(n, m)

  pool.AddRangeJob(0, len(sequences), jg, func(i int, pool threadpool.ThreadPool, erf func() error) error {
    result[i] = scanSequence(config, kmersMeta, sequences[i])
    return nil
  })

  granges.AddMeta("kmers", result)
  WriteResult(config, kmersMeta, granges, filenameOut)
}

/* -------------------------------------------------------------------------- */

func main() {

  config  := Config{}
  options := getopt.New()

  optRegions := options. StringLong("regions",     0 , "", "bed with with regions")
  optHeader  := options.   BoolLong("header",      0 ,     "print kmer header")
  optHuman   := options.   BoolLong("human",       0 ,     "print human readable kmer statistics")
  optThreads := options.    IntLong("threads",     0 ,  1, "number of threads [default: 1]")
  optVerbose := options.CounterLong("verbose",    'v',     "verbose level [-v or -vv]")
  optHelp    := options.   BoolLong("help",       'h',     "print help")

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
  config.Header  = *optHeader
  config.Human   = *optHuman
  config.Threads = *optThreads
  config.Verbose = *optVerbose
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
