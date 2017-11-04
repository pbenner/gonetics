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
import   "bufio"
import   "log"
import   "io"
import   "math"
import   "os"
import   "strconv"
import   "strings"

import   "github.com/pborman/getopt"
import . "github.com/pbenner/gonetics"

/* i/o
 * -------------------------------------------------------------------------- */

func PrintStderr(verbose int, level int, format string, args ...interface{}) {
  if verbose >= level {
    fmt.Fprintf(os.Stderr, format, args...)
  }
}

/* -------------------------------------------------------------------------- */

func importGenome(filename string, verbose int) Genome {
  genome := Genome{}
  PrintStderr(verbose, 1, "Reading genome `%s'... ", filename)
  if err := genome.Import(filename); err != nil {
    PrintStderr(verbose, 1, "failed\n")
    log.Fatal(err)
  }
  PrintStderr(verbose, 1, "done\n")
  return genome
}

func exportTrack(track Track, trackFilename string, verbose int) {
  PrintStderr(verbose, 1, "Writing track `%s'... ", trackFilename)
  if err := (GenericTrack{track}).ExportBigWig(trackFilename); err != nil {
    PrintStderr(verbose, 1, "failed\n")
    log.Fatal(err)
  } else {
    PrintStderr(verbose, 1, "done\n")
  }
}

/* -------------------------------------------------------------------------- */

func importTable(r io.Reader, track SimpleTrack, columnName string) error {
  scanner := bufio.NewScanner(r)
  colIdx  := -1
  seq     := TrackMutableSequence{}

  // read sequence name and get track sequence
  if !scanner.Scan() {
    return fmt.Errorf("invalid table (no sequence information)")
  } else {
    fields := strings.Fields(scanner.Text())
    if len(fields) != 2 {
      return fmt.Errorf("invalid table (first line must have two entries")
    }
    if s, err := track.GetMutableSequence(fields[1]); err != nil {
      return err
    } else {
      seq = s
    }
  }
  // read header
  if !scanner.Scan() {
    return fmt.Errorf("invalid table (no header information)")
  } else {
    fields := strings.Fields(scanner.Text())
    for i, str := range fields {
      if str == columnName {
        colIdx = i; break
      }
    }
    if colIdx == -1 {
      return fmt.Errorf("column `%s' not found in table", columnName)
    }
  }
  for i := 0; scanner.Scan(); i++ {
    fields := strings.Fields(scanner.Text())
    if colIdx >= len(fields) {
      return fmt.Errorf("table has invalid number of columns at line", i+3)
    }
    if i >= seq.NBins() {
      return fmt.Errorf("table has more entries than expected (is the binsize correct?)")
    }
    if v, err := strconv.ParseFloat(fields[colIdx], 64); err != nil {
      return fmt.Errorf("parsing floating point number at line `%d' failed: %v", i+3, err)
    } else {
      seq.SetBin(i, v)
    }
  }
  return nil
}

/* -------------------------------------------------------------------------- */

func chromHmmTablesToBigWig(filenameGenome string, binSize int, columnName string, trackInit float64, filenameOut string, filenames []string, verbose int) {

  track := AllocSimpleTrack("", importGenome(filenameGenome, verbose), binSize)

  // initialize track to NaN
  if trackInit != 0.0 {
    GenericMutableTrack{track}.Map(track, func(seqname string, pos int, value float64) float64 {
      return math.NaN()
    })
  }

  for _, filename := range filenames {
    PrintStderr(verbose, 1, "Reading table `%s'... ", filename)
    // open file
    f, err := os.Open(filename)
    if err != nil {
      PrintStderr(verbose, 1, "failed\n")
      log.Fatal(err)
    }
    if err := importTable(f, track, columnName); err != nil {
      PrintStderr(verbose, 1, "failed\n")
      log.Fatal(err)
    }      
    PrintStderr(verbose, 1, "done\n")

    f.Close()
  }
  exportTrack(track, filenameOut, verbose)
}

/* -------------------------------------------------------------------------- */

func main() {

  trackInit := 0.0
  verbose   := 0

  options := getopt.New()

  optTrackInit := options. StringLong("initial-value",  0 , "", "track initial value [default: 0]")
  optVerbose   := options.CounterLong("verbose",       'v',     "verbose level [-v or -vv]")
  optHelp      := options.   BoolLong("help",          'h',     "print help")

  options.SetParameters("<GENOME_FILE> <BINSIZE> <COLUMN_NAME> <RESULT_FILE> <FILES...>")
  options.Parse(os.Args)

  // parse options
  //////////////////////////////////////////////////////////////////////////////
  if *optHelp {
    options.PrintUsage(os.Stdout)
    os.Exit(0)
  }
  if *optVerbose != 0 {
    verbose = *optVerbose
  }
  if len(options.Args()) < 5 {
    options.PrintUsage(os.Stderr)
    os.Exit(1)
  }
  if *optTrackInit != "" {
    v, err := strconv.ParseFloat(*optTrackInit, 64)
    if err != nil {
      log.Fatalf("parsing initial value failed: %v", err)
    }
    trackInit = v
  }

  filenameGenome := options.Args()[0]
  binSize        := 0
  columnName     := options.Args()[2]
  filenameOut    := options.Args()[3]
  filenames      := options.Args()[4:]

  if t, err := strconv.ParseInt(options.Args()[1], 10, 64); err != nil {
    log.Fatal(err)
  } else {
    binSize = int(t)
  }
  if binSize < 1 {
    log.Fatalf("invalid binsize `%d' specified", binSize)
  }
  chromHmmTablesToBigWig(filenameGenome, binSize, columnName, trackInit, filenameOut, filenames, verbose)
}
