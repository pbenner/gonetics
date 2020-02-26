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
import   "os"
import   "strings"

import   "github.com/pborman/getopt"

import . "github.com/pbenner/gonetics"

/* -------------------------------------------------------------------------- */

type Config struct {
  InputFormat  string
  Verbose      int
}

/* -------------------------------------------------------------------------- */

func PrintStderr(config Config, level int, format string, args ...interface{}) {
  if config.Verbose >= level {
    fmt.Fprintf(os.Stderr, format, args...)
  }
}

/* -------------------------------------------------------------------------- */

func importBed3(config Config, filename string) GRanges {
  granges := GRanges{}
  if filename == "" {
    if err := granges.ReadBed3(os.Stdin); err != nil {
      log.Fatal(err)
    }
  } else {
    PrintStderr(config, 1, "Reading bed file `%s'... ", filename)
    if err := granges.ImportBed3(filename); err != nil {
      PrintStderr(config, 1, "failed\n")
      log.Fatal(err)
    } else {
      PrintStderr(config, 1, "done\n")
    }
  }
  return granges
}

func importBed6(config Config, filename string) GRanges {
  granges := GRanges{}
  if filename == "" {
    if err := granges.ReadBed6(os.Stdin); err != nil {
      log.Fatal(err)
    }
  } else {
    PrintStderr(config, 1, "Reading bed file `%s'... ", filename)
    if err := granges.ImportBed6(filename); err != nil {
      PrintStderr(config, 1, "failed\n")
      log.Fatal(err)
    } else {
      PrintStderr(config, 1, "done\n")
    }
  }
  return granges
}

func importBed9(config Config, filename string) GRanges {
  granges := GRanges{}
  if filename == "" {
    if err := granges.ReadBed9(os.Stdin); err != nil {
      log.Fatal(err)
    }
  } else {
    PrintStderr(config, 1, "Reading bed file `%s'... ", filename)
    if err := granges.ImportBed9(filename); err != nil {
      PrintStderr(config, 1, "failed\n")
      log.Fatal(err)
    } else {
      PrintStderr(config, 1, "done\n")
    }
  }
  return granges
}

func importTable(config Config, filename string) GRanges {
  granges := GRanges{}
  if filename == "" {
    if err := granges.ReadTableAll(os.Stdin); err != nil {
      log.Fatal(err)
    }
  } else {
    PrintStderr(config, 1, "Reading table file `%s'... ", filename)
    if err := granges.ImportTableAll(filename); err != nil {
      PrintStderr(config, 1, "failed\n")
      log.Fatal(err)
    } else {
      PrintStderr(config, 1, "done\n")
    }
  }
  return granges
}

func importInput(config Config, filename string) GRanges {
  switch strings.ToLower(config.InputFormat) {
  case "bed3" : return importBed3 (config, filename)
  case "bed6" : return importBed6 (config, filename)
  case "bed9" : return importBed9 (config, filename)
  case "table": return importTable(config, filename)
  default:
    log.Fatalf("invalid input format: %s", config.InputFormat)
    panic("internal error")
  }
}

/* -------------------------------------------------------------------------- */

func exportBed3(config Config, filename string, granges GRanges) {
  if filename == "" {
    if err := granges.WriteBed3(os.Stdout); err != nil {
      log.Fatal(err)
    }
  } else {
    PrintStderr(config, 1, "Exporting bed file to `%s'... ", filename)
    if err := granges.ExportBed3(filename, false); err != nil {
      PrintStderr(config, 1, "failed\n")
      log.Fatal(err)
    } else {
      PrintStderr(config, 1, "done\n")
    }
  }
}

func exportBed6(config Config, filename string, granges GRanges) {
  if filename == "" {
    if err := granges.WriteBed6(os.Stdout); err != nil {
      log.Fatal(err)
    }
  } else {
    PrintStderr(config, 1, "Exporting bed file to `%s'... ", filename)
    if err := granges.ExportBed6(filename, false); err != nil {
      PrintStderr(config, 1, "failed\n")
      log.Fatal(err)
    } else {
      PrintStderr(config, 1, "done\n")
    }
  }
}

func exportBed9(config Config, filename string, granges GRanges) {
  if filename == "" {
    if err := granges.WriteBed3(os.Stdout); err != nil {
      log.Fatal(err)
    }
  } else {
    PrintStderr(config, 1, "Exporting bed file to `%s'... ", filename)
    if err := granges.ExportBed9(filename, false); err != nil {
      PrintStderr(config, 1, "failed\n")
      log.Fatal(err)
    } else {
      PrintStderr(config, 1, "done\n")
    }
  }
}

func exportTable(config Config, filename string, granges GRanges) {
  if filename == "" {
    if err := granges.WriteTable(os.Stdout, true, true, false); err != nil {
      log.Fatal(err)
    }
  } else {
    PrintStderr(config, 1, "Exporting bed file to `%s'... ", filename)
    if err := granges.ExportTable(filename, true, true, false); err != nil {
      PrintStderr(config, 1, "failed\n")
      log.Fatal(err)
    } else {
      PrintStderr(config, 1, "done\n")
    }
  }
}

func exportOutput(config Config, filename string, granges GRanges) {
  switch strings.ToLower(config.InputFormat) {
  case "bed3" : exportBed3 (config, filename, granges)
  case "bed6" : exportBed6 (config, filename, granges)
  case "bed9" : exportBed9 (config, filename, granges)
  case "table": exportTable(config, filename, granges)
  default:
    log.Fatalf("invalid input format: %s", config.InputFormat)
  }
}

/* -------------------------------------------------------------------------- */

func removeOverlaps(config Config, filenameRm, filenameIn, filenameOut string) {
  regionsRm := importBed3 (config, filenameRm)
  regionsIn := importInput(config, filenameIn)
  regionsIn  = regionsIn.RemoveOverlapsWith(regionsRm)

  exportOutput(config, filenameOut, regionsIn)
}

/* -------------------------------------------------------------------------- */

func main() {

  config  := Config{}

  options := getopt.New()

  optInFormat  := options. StringLong("input-format",  0 , "bed3", "bed3 (default), bed6, or bed9")
  optHelp      := options.   BoolLong("help",         'h',         "print help")
  optVerbose   := options.CounterLong("verbose",      'v',         "be verbose")

  options.SetParameters("<INADMISSIBLE_REGIONS.bed> [INPUT.bed [OUTPUT.bed]]\n")
  options.Parse(os.Args)

  if *optHelp {
    options.PrintUsage(os.Stdout)
    os.Exit(0)
  }
  if len(options.Args()) != 1 && len(options.Args()) != 2 && len(options.Args()) != 3 {
    options.PrintUsage(os.Stderr)
    os.Exit(1)
  }
  config.Verbose     = *optVerbose
  config.InputFormat = *optInFormat

  filenameRm  := ""
  filenameIn  := ""
  filenameOut := ""
  if len(options.Args()) > 0 {
    filenameRm  = options.Args()[0]
  }
  if len(options.Args()) > 1 {
    filenameIn  = options.Args()[1]
  }
  if len(options.Args()) > 2 {
    filenameOut = options.Args()[2]
  }
  removeOverlaps(config, filenameRm, filenameIn, filenameOut)
}
