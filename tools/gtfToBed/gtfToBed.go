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
import   "log"
import   "os"
//import   "strconv"

import   "github.com/pborman/getopt"

import . "github.com/pbenner/gonetics"

/* -------------------------------------------------------------------------- */

func printMsg(verbose bool, format string, args... interface{}) {
  if verbose {
    fmt.Fprintf(os.Stderr, format, args...)
  }
}

/* -------------------------------------------------------------------------- */

type gtfEntry struct {
  seqname string
  from    int
  to      int
  strand  byte
}

func (obj *gtfEntry) merge(r gtfEntry) error {
  if obj.seqname != r.seqname {
    return fmt.Errorf("cannot merge entries on different sequences")
  }
  if r.from < obj.from {
    obj.from = r.from
  }
  if r.to > obj.to {
    obj.to = r.to
  }
  // if strands do not match, set to `*'
  if r.strand != obj.strand {
    obj.strand = '*'
  }
  return nil
}

/* -------------------------------------------------------------------------- */

func mergeRows(g GRanges, mergeBy string) GRanges {
  m := make(map[string]gtfEntry)
  a := g.GetMetaStr(mergeBy)

  if a == nil {
    log.Fatalf("attribute `%s' missing", mergeBy)
  }

  for i := 0; i < g.Length(); i++ {
    // skip rows with no attribute
    if a[i] == ""{
      continue
    }
    r := gtfEntry{
      seqname: g.Seqnames[i],
      from   : g.Ranges[i].From,
      to     : g.Ranges[i].To,
      strand : g.Strand[i] }
    if entry, ok := m[a[i]]; ok {
      // attribute already present, merge both
      if err := entry.merge(r); err != nil {
        log.Fatal(err)
      }
    } else {
      // add new attribute
      m[a[i]] = r
    }
  }
  // convert map to new granges object
  seqnames := []string{}
  from     := []int{}
  to       := []int{}
  strand   := []byte{}
  names    := []string{}
  for name, entry := range m {
    names    = append(names,    name)
    seqnames = append(seqnames, entry.seqname)
    from     = append(from,     entry.from)
    to       = append(to,       entry.to)
    strand   = append(strand,   entry.strand)
  }
  r := NewGRanges(seqnames, from, to, strand)
  r.AddMeta("name", names)
  return r
}

/* -------------------------------------------------------------------------- */

func gtfToBed(filenameIn, filenameOut, mergeBy string, verbose bool) {
  attrNames := []string{}
  attrTypes := []string{}
  attrDef   := []interface{}{}
  // add mergeBy attribute if it is not listed
  if mergeBy != "" {
    attrNames = append(attrNames, mergeBy)
    attrTypes = append(attrTypes, "[]string")
    attrDef   = append(attrDef,   "")
  }

  g := GRanges{}
  printMsg(verbose, "Reading gtf from file `%s'... ", filenameIn)
  if err := g.ImportGTF(filenameIn, attrNames, attrTypes, attrDef); err != nil {
    printMsg(verbose, "failed\n")
    log.Fatal(err)
  }
  printMsg(verbose, "done\n")

  // rename score since it has the wrong type for bed
  g.RenameMeta("score", "gtfScore")

  if mergeBy == "" {
    // use feature as name colum
    g.RenameMeta("feature", "name")
  } else {
    // merge rows by optional attributes and use these
    // as name column
    g = mergeRows(g, mergeBy)
  }

  printMsg(verbose, "Writing result to file `%s'... ", filenameOut)
  if err := g.ExportBed6(filenameOut, false); err != nil {
    printMsg(verbose, "failed\n")
    log.Fatal(err)
  }
  printMsg(verbose, "done\n")
}

/* -------------------------------------------------------------------------- */

func main() {
  options := getopt.New()
  options.SetProgram(fmt.Sprintf("%s", os.Args[0]))

  optHelp       := options.  BoolLong("help",        'h',     "print help")
  optVerbose    := options.  BoolLong("verbose",     'v',     "be verbose")
  optMergeBy    := options.StringLong("merge-by",     0 , "", "merge rows by optional attribut [e.g. transcript_id, gene_id, gene_name]")

  options.SetParameters("<input.gtf> <output.bed>")
  options.Parse(os.Args)

  if *optHelp {
    options.PrintUsage(os.Stdout)
    os.Exit(0)
  }
  if len(options.Args()) != 2 {
    options.PrintUsage(os.Stderr)
    os.Exit(1)
  }
  filenameIn  := options.Args()[0]
  filenameOut := options.Args()[1]

  gtfToBed(filenameIn, filenameOut, *optMergeBy, *optVerbose)
}