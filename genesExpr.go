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

package gonetics

/* -------------------------------------------------------------------------- */

import "fmt"
import "bufio"
import "compress/gzip"
import "os"
import "strconv"
import "strings"
import "unicode"

import . "github.com/pbenner/pshape/Utility"

/* i/o
 * -------------------------------------------------------------------------- */

type gtfOptional struct {
  GeneId       string
  TranscriptId string
  FPKM         float64
  RPKM         float64
}

func readGTFParseOptional(fields []string) gtfOptional {
  if len(fields) % 2 == 1 {
    panic("ReadGTF(): invalid file format!")
  }
  gtfOpt := gtfOptional{}
  // loop through list
  for i := 0; i < len(fields); i += 2 {
    if fields[i] == "gene_id" {
      gtfOpt.GeneId = fields[i+1]
    }
    if fields[i] == "transcript_id" {
      gtfOpt.TranscriptId = fields[i+1]
    }
    if fields[i] == "FPKM" {
      t, err := strconv.ParseFloat(fields[i+1], 64); Check(err)
      gtfOpt.FPKM = t
    }
    if fields[i] == "RPKM" {
      t, err := strconv.ParseFloat(fields[i+1], 64); Check(err)
      gtfOpt.RPKM = t
    }
  }
  return gtfOpt
}

func readGTFParseLine(line string) []string {
  // if quoted
  q := false
  f := func(r rune) bool {
    if r == '"' {
      q = !q
    }
    // A quote is treated as a white space so that it is removed from the
    // line. Otherwise a white space is removed only if q (quote) is false.
    return r == '"' || ((unicode.IsSpace(r) || r == ';') && q == false)
  }
  return strings.FieldsFunc(line, f)
}

// Parse expression data from a GTF file (gene transfer format). The data
// is added as a meta column named "expr" to the gene list. Parameters:
//  geneIdName: Name of the optional field containing the gene id
//  exprIdName: Name of the optional field containing the expression data
//  genes: List of query genes
func (genes *Genes) ReadGTF(filename, geneIdName, exprIdName string, verbose bool) {
  var scanner *bufio.Scanner
  // open file
  f, err := os.Open(filename)
  Check(err)
  defer f.Close()
  // check if file is gzipped
  if IsGzip(filename) {
    g, err := gzip.NewReader(f)
    Check(err)
    defer g.Close()
    scanner = bufio.NewScanner(g)
  } else {
    scanner = bufio.NewScanner(f)
  }
  // it seems that buffering the data does not increase
  // performance
  expr := make([]float64, genes.Length())

  for scanner.Scan() {
    Check(scanner.Err())
    geneTmp := ""
    exprTmp := 0.0
    fields := readGTFParseLine(scanner.Text())
    if len(fields) == 0 {
      continue
    }
    if len(fields) <= 8 {
      panic("File must have more than eight columns!")
    }
    if fields[2] != "transcript" && fields[2] != "TSS" {
      continue
    }
    // parse optional fields to get the gene/transcript id and
    // expression level
    gtfOpt := readGTFParseOptional(fields[8:len(fields)])

    switch geneIdName {
    case "gene_id"      : geneTmp = gtfOpt.GeneId
    case "transcript_id": geneTmp = gtfOpt.TranscriptId
    default: panic("Invalid gene_id/trascript_id!")
    }
    switch exprIdName {
    case "FPKM": exprTmp = gtfOpt.FPKM
    case "RPKM": exprTmp = gtfOpt.RPKM
    default: panic("Invalid expression type!")
    }
    // parse gene/transcript list
    fields = strings.FieldsFunc(geneTmp, func(r rune) bool { return r == ',' })
    for _, gene := range fields {
      if i, ok := genes.FindGene(gene); ok {
        expr[i] += exprTmp
      } else {
        if verbose {
          fmt.Fprintf(os.Stderr, "`%s' not present in gene list!\n", gene)
        }
      }
    }
  }
  genes.AddMeta("expr", expr)
}

// Import expression data from cufflinks. The data is added to the gene
// list as a meta column named "expr".
func (genes *Genes) ReadCufflinksFPKMTracking(filename string, verbose bool) {
  var scanner *bufio.Scanner
  // open file
  f, err := os.Open(filename)
  Check(err)
  defer f.Close()
  // check if file is gzipped
  if IsGzip(filename) {
    g, err := gzip.NewReader(f)
    Check(err)
    defer g.Close()
    scanner = bufio.NewScanner(g)
  } else {
    scanner = bufio.NewScanner(f)
  }
  expr := make([]float64, genes.Length())

  // parse header
  scanner.Scan();
  Check(scanner.Err())
  header := strings.Fields(scanner.Text())
  if len(header) != 13 || header[0] != "tracking_id" ||
    (header[9] != "FPKM" && header[9] != "RPKM") {
    panic("Invalid header!")
  }
  // parse data
  for scanner.Scan() {
    Check(scanner.Err())

    fields := strings.Fields(scanner.Text())
    if len(fields) == 0 {
      continue
    }
    if len(fields) != 13 {
      panic("File must have 13 columns!")
    }
    if fields[12] == "OK" {
      geneStr := fields[0]
      exprStr := fields[9]
      if i, ok := genes.FindGene(geneStr); ok {
        t, err := strconv.ParseFloat(exprStr, 64); Check(err)
        expr[i] += t
      } else {
        if verbose {
          fmt.Fprintf(os.Stderr, "`%s' not present in gene list!\n", geneStr)
        }
      }
    }
  }
  genes.AddMeta("expr", expr)
}
