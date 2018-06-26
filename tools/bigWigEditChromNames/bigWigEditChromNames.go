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
import   "encoding/binary"
import   "log"
import   "os"
import   "regexp"
import   "strings"

import   "github.com/pborman/getopt"

import . "github.com/pbenner/gonetics"

/* -------------------------------------------------------------------------- */

func editChromNames(filename, regex, repl string, dryRun, verbose bool) {
  bwf := new(BbiFile)
  {
    f, err := os.Open(filename); if err != nil {
      log.Fatal(err)
    }
    // read chromosome list
    if err := bwf.Open(f); err != nil {
      log.Fatal(err)
    }
    f.Close()
  }

  // file pointer for writing
  fptr, err := os.OpenFile(filename, os.O_RDWR, 0); if err != nil {
    log.Fatal(err)
  }

  r, err := regexp.Compile(regex); if err != nil {
    log.Fatal("invalid regular expression:", err)
  }
  for i := 0; i < len(bwf.ChromData.Keys); i++ {
    if len(bwf.ChromData.Values[i]) != 8 {
      log.Fatalf("reading `%s' failed: invalid chromosome list", filename)
    }
    // read seqname
    seqname_old := strings.TrimRight(string(bwf.ChromData.Keys[i]), "\x00")
    // apply regular expression
    seqname_new := r.ReplaceAllString(seqname_old, repl)
    if dryRun {
      fmt.Printf("`%s' -> `%s'\n", seqname_old, seqname_new)
    } else {
      // update ChromData
      copy(bwf.ChromData.Keys[i], seqname_new)
      // write new seqname
      if _, err := fptr.Seek(bwf.ChromData.PtrKeys[i], 0); err != nil {
        log.Fatal(err)
      }
      if err := binary.Write(fptr, bwf.Order, bwf.ChromData.Keys[i]); err != nil {
        log.Fatal("writing new seqname failed:", err)
      }
    }
  }
}

func main() {
  options := getopt.New()
  options.SetProgram(fmt.Sprintf("%s", os.Args[0]))

  optDryRun     := options.  BoolLong("dry-run",      0 ,     "just print changes and do not edit the file")
  optHelp       := options.  BoolLong("help",        'h',     "print help")
  optVerbose    := options.  BoolLong("verbose",     'v',     "be verbose")

  options.SetParameters("<input.bw> <regex> <replacement>")
  options.Parse(os.Args)

  if *optHelp {
    options.PrintUsage(os.Stdout)
    os.Exit(0)
  }
  if len(options.Args()) != 3 {
    options.PrintUsage(os.Stderr)
    os.Exit(1)
  }
  filename := options.Args()[0]
  regex    := options.Args()[1]
  repl     := options.Args()[2]

  editChromNames(filename, regex, repl, *optDryRun, *optVerbose)
}
