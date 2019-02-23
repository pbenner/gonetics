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

import   "github.com/pborman/getopt"

import . "github.com/pbenner/gonetics"
import   "github.com/pbenner/threadpool"

/* -------------------------------------------------------------------------- */

type Config struct {
  Alphabet NucleotideAlphabet
  Pretty   bool
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

/* -------------------------------------------------------------------------- */

func prettyPrintKmer(config Config, result []int, k int, p []int) {
  c1 := make([]byte, k)
  c2 := make([]byte, k)
  for i := 0; i < len(result); i++ {
    // convert index to sequence
    for j, ix := 0, i; j < k; j++ {
      c1[k-j-1] = byte(ix % config.Alphabet.Length())
      ix        = ix / config.Alphabet.Length()
      if x, err := config.Alphabet.ComplementCoded(c1[k-j-1]); err != nil {
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
        if x, err := config.Alphabet.Decode(c1[j]); err != nil {
          log.Fatal(err)
        } else {
          c1[j] = x
        }
        if x, err := config.Alphabet.Decode(c2[j]); err != nil {
          log.Fatal(err)
        } else {
          c2[j] = x
        }
      }
      fmt.Printf("%s|%s = %d\n", string(c1), string(c2), result[i])
    }
  }
}

func prettyPrint(config Config, result [][]int, n, m int, p []int) {
  for k := n; k <= m; k++ {
    prettyPrintKmer(config, result[k-n], k, p)
  }
}

/* -------------------------------------------------------------------------- */

func scanSequence(config Config, sequence []byte, n, m int, p []int) [][]int {
  c1 := make([]int, len(sequence))
  c2 := make([]int, len(sequence))
  for i := 0; i < len(sequence); i++ {
    if r, err := config.Alphabet.Code(sequence[i]); err != nil {
      log.Fatal(err)
    } else {
      c1[i] = int(r)
    }
    if r, err := config.Alphabet.ComplementCoded(byte(c1[i])); err != nil {
      log.Fatal(err)
    } else {
      c2[i] = int(r)
    }
  }
  // allocate results matrix
  result := make([][]int, m-n+1)
  for k := n; k <= m; k++ {
    result[k-n] = make([]int, ipow(config.Alphabet.Length(), k))
  }
  // loop over sequence
  for i := 0; i < len(sequence); i++ {
    // loop over all k-mers
    for k := n; k <= m && i+k-1 < len(sequence); k++ {
      // eval k-mer
      s := 0
      r := 0
      for j := 0; j < k; j++ {
        s += c1[i+j] * p[k-j-1]
        r += c2[i+j] * p[j]
      }
      if s < r {
        result[k-n][s] += 1
      } else {
        result[k-n][r] += 1
      }
    }
  }
  return result
}

/* -------------------------------------------------------------------------- */

func kmerSearch(config Config, n, m int, filenameFasta, filenameOut string) {
  pool := threadpool.New(config.Threads, 100*config.Threads)
  jg   := pool.NewJobGroup()
  ss   := ImportFasta(config, filenameFasta)

  // evaluate powers 4^k
  p := make([]int, m+1)
  for k := 0; k <= m; k++ {
    p[k] = ipow(config.Alphabet.Length(), k)
  }
  pool.AddRangeJob(0, len(ss.Seqnames), jg, func(i int, pool threadpool.ThreadPool, erf func() error) error {
    name     := ss.Seqnames[i]
    sequence := ss.Sequences[name]
    result   := scanSequence(config, sequence, n, m, p)
    prettyPrint(config, result, n, m, p)
    return nil
  })
}

/* -------------------------------------------------------------------------- */

func main() {

  config  := Config{}
  options := getopt.New()

  optPretty  := options.   BoolLong("pretty",      0 ,         "print pretty")
  optThreads := options.    IntLong("threads",     0 ,      1, "number of threads [default: 1]")
  optVerbose := options.CounterLong("verbose",    'v',         "verbose level [-v or -vv]")
  optHelp    := options.   BoolLong("help",       'h',         "print help")

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
  config.Pretty  = *optPretty
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
  if len(options.Args()) == 3 {
    filenameFasta = options.Args()[2]
  }
  if len(options.Args()) == 4 {
    filenameOut   = options.Args()[3]
  }
  kmerSearch(config, int(n), int(m), filenameFasta, filenameOut)
}
