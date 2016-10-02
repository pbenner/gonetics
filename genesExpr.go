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

/* i/o
 * -------------------------------------------------------------------------- */

// Parse expression data from a GTF file (gene transfer format). The data
// is added as a meta column named "expr" to the gene list. Parameters:
//  geneIdName: Name of the optional field containing the gene id
//  exprIdName: Name of the optional field containing the expression data
//  genes: List of query genes
func (genes *Genes) ReadGTFExpr(filename, geneIdName, exprIdName string) error {
  granges := GRanges{}
  granges.ReadGTF(filename, []string{geneIdName, exprIdName}, []string{"[]string", "[]float64"})

  // slice containing expression values
  expr := make([]float64, genes.Length())

  geneIds := granges.GetMetaStr(geneIdName)
  exprVal := granges.GetMetaFloat(exprIdName)

  if len(geneIds) == 0 {
    return fmt.Errorf("invalid geneIdName `%s'", geneIdName)
  }
  if len(exprVal) == 0 {
    return fmt.Errorf("invalid exprIdName `%s'", exprIdName)
  }
  for i := 0; i < granges.Length(); i++ {
    if j, ok := genes.FindGene(geneIds[i]); ok {
      expr[j] += exprVal[i]
    }
  }
  genes.AddMeta("expr", expr)

  return nil
}

// Import expression data from cufflinks. The data is added to the gene
// list as a meta column named "expr".
func (genes *Genes) ReadCufflinksFPKMTracking(filename string, verbose bool) error {
  var scanner *bufio.Scanner
  // open file
  f, err := os.Open(filename)
  if err != nil {
    return err
  }
  defer f.Close()
  // check if file is gzipped
  if isGzip(filename) {
    g, err := gzip.NewReader(f)
    if err != nil {
      return err
    }
    defer g.Close()
    scanner = bufio.NewScanner(g)
  } else {
    scanner = bufio.NewScanner(f)
  }
  expr := make([]float64, genes.Length())

  // parse header
  scanner.Scan();
  if err := scanner.Err(); err != nil {
    return err
  }
  header := strings.Fields(scanner.Text())
  if len(header) != 13 || header[0] != "tracking_id" ||
    (header[9] != "FPKM" && header[9] != "RPKM") {
    return fmt.Errorf("invalid header")
  }
  // parse data
  for scanner.Scan() {
    if err := scanner.Err(); err != nil {
      return err
    }
    fields := strings.Fields(scanner.Text())
    if len(fields) == 0 {
      continue
    }
    if len(fields) != 13 {
      return fmt.Errorf("file must have 13 columns")
    }
    if fields[12] == "OK" {
      geneStr := fields[0]
      exprStr := fields[9]
      if i, ok := genes.FindGene(geneStr); ok {
        t, err := strconv.ParseFloat(exprStr, 64)
        if err != nil {
          return err
        }
        expr[i] += t
      } else {
        if verbose {
          fmt.Fprintf(os.Stderr, "`%s' not present in gene list!\n", geneStr)
        }
      }
    }
  }
  genes.AddMeta("expr", expr)

  return nil
}
