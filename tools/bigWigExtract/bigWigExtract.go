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
import   "math"
import   "os"
import   "path/filepath"
import   "strconv"

import   "github.com/pborman/getopt"

import . "github.com/pbenner/gonetics"

/* -------------------------------------------------------------------------- */

type Config struct {
  Verbose    int
  Scientific bool
  BinSize    int
  BinOver    int
  TrackInit  float64
  BinStat    BinSummaryStatistics
}

/* -------------------------------------------------------------------------- */

func PrintStderr(config Config, level int, format string, args ...interface{}) {
  if config.Verbose >= level {
    fmt.Fprintf(os.Stderr, format, args...)
  }
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

/* -------------------------------------------------------------------------- */

func importBed3(config Config, filename string) GRanges {
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

func exportTable(config Config, granges GRanges, filename string) {
  if filename == "" {
    if err := granges.WriteTable(os.Stdout, true, true, false, OptionPrintScientific{config.Scientific}); err != nil {
      PrintStderr(config, 1, "failed\n")
      log.Fatal(err)
    }
  } else {
    PrintStderr(config, 1, "Writing table `%s'... ", filename)
    if err := granges.ExportTable(filename, true, true, false, OptionPrintScientific{config.Scientific}); err != nil {
      PrintStderr(config, 1, "failed\n")
      log.Fatal(err)
    } else {
      PrintStderr(config, 1, "done\n")
    }
  }
}

/* -------------------------------------------------------------------------- */

func extractTable(config Config, granges GRanges, filenameBw, filenameOut string) {
  if err := granges.ImportBigWig(filenameBw, "counts", config.BinStat, config.BinSize, config.BinOver, config.TrackInit, false); err != nil {
    log.Fatal(err)
  }
  exportTable(config, granges, filenameOut)
}

func extractBigWig(config Config, granges GRanges, filenameBw, filenameOut string) {
  f, err := os.Open(filenameBw)
  if err != nil {
    log.Fatal(err)
  }
  r, err := NewBigWigReader(f)
  if err != nil {
    log.Fatal(err)
  }
  var track *SimpleTrack

  for i := 0; i < granges.Length(); i++ {
    if src, binSize, err := r.QuerySlice(granges.Seqnames[i], granges.Ranges[i].From, granges.Ranges[i].To, config.BinStat, config.BinSize, config.BinOver, config.TrackInit); err != nil {
      log.Fatal(err)
    } else {
      config.BinSize = binSize
      if track == nil {
        t := AllocSimpleTrack("", r.Genome, binSize)
        // initialize track
        GenericMutableTrack{t}.Map(t, func(string, int, float64) float64 {
          return math.NaN()
        })
        track = &t
      }
      if dst, err := track.GetMutableSequence(granges.Seqnames[i]); err != nil {
        log.Fatal(err)
      } else {
        // copy result to track
        for j := granges.Ranges[i].From; j < granges.Ranges[i].To; j += binSize {
          dst.Set(j, src[(j-granges.Ranges[i].From)/binSize])
        }
      }
    }
  }
  if err := track.ExportBigWig(filenameOut); err != nil {
    log.Fatal(err)
  }
}

func extract(config Config, filenameBed, filenameBw, filenameOut string) {
  // import bed file first to check if it exists
  granges := importBed3(config, filenameBed)

  switch filepath.Ext(filenameOut) {
  case ".bw"    : fallthrough
  case ".bigWig":
    extractBigWig(config, granges, filenameBw, filenameOut)
  case ""       : fallthrough
  case ".table" :
    extractTable(config, granges, filenameBw, filenameOut)
  default:
    log.Fatal("output file has invalid file extension")
  }
}

/* -------------------------------------------------------------------------- */

func main() {

  config  := Config{}

  options := getopt.New()

  optBinSize    := options.    IntLong("bin-size",       0 ,      0, "bin size")
  optBinStat    := options. StringLong("bin-summary",    0 , "mean", "bin summary statistic [mean (default), max, min, discrete mean]")
  optBinOver    := options.    IntLong("bin-overlap",    0 ,      0, "number of overlapping bins when computing the summary")
  optTrackInit  := options. StringLong("initial-value",  0 ,     "", "track initial value [default: 0]")
  optScientific := options.   BoolLong("scientific",     0 ,         "use scientific format to print values")
  optHelp       := options.   BoolLong("help",          'h',         "print help")
  optVerbose    := options.CounterLong("verbose",       'v',         "be verbose")

  options.SetParameters(
    "<INPUT.bw> <INPUT.bed> [<OUTPUT.(table|bw)]>\n\n" +
    "  The extracted data is saved as a table if the output file has extension\n"  +
    "  `.table'. If the extension is `.bw' or '.bigWig' then the output format\n"  +
    "  is bigWig. The data is printed as a table to stdout if no output file is\n" +
    "  given.\n")
  options.Parse(os.Args)

  if *optHelp {
    options.PrintUsage(os.Stdout)
    os.Exit(0)
  }
  if len(options.Args()) != 2 && len(options.Args()) != 3 {
    options.PrintUsage(os.Stderr)
    os.Exit(1)
  }
  if *optTrackInit != "" {
    v, err := strconv.ParseFloat(*optTrackInit, 64)
    if err != nil {
      log.Fatalf("parsing initial value failed: %v", err)
    }
    config.TrackInit = v
  } else {
    config.TrackInit = 0.0
  }
  config.Verbose    = *optVerbose
  config.Scientific = *optScientific
  config.BinSize    = *optBinSize
  config.BinOver    = *optBinOver
  config.BinStat    = getBinSummaryStatistics(*optBinStat)

  filenameBw  := options.Args()[0]
  filenameBed := options.Args()[1]
  filenameOut := ""

  if len(options.Args()) == 3 {
    filenameOut = options.Args()[2]
  }
  extract(config, filenameBed, filenameBw, filenameOut)
}
