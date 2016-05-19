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

import "bufio"
import "bytes"
import "fmt"
import "log"
import "compress/gzip"
import "database/sql"
import "io"
import "io/ioutil"
import "os"
import "sort"
import "strconv"
import "strings"

import _ "github.com/go-sql-driver/mysql"
import . "github.com/pbenner/pshape/Utility"

/* -------------------------------------------------------------------------- */

// Container for genes. Tx contains the transcription start and end
// positions. Cds specifies the coding regions.
type Genes struct {
  Names      []string
  Seqnames   []string
  Tx         []Range
  Cds        []Range
  Strand     []byte
  index      map[string]int
  Meta
}

/* constructors
 * -------------------------------------------------------------------------- */

func newGenes(names, seqnames []string, tx, cds []Range, strand []byte) Genes {
  n := len(names)
  if len(tx)       != n || len(cds)    != n  ||
    (len(seqnames) != n || len(strand) != n) {
    panic("NewGenes(): invalid arguments!")
  }
  index := map[string]int{}
  for i := 0; i < n; i++ {
    // check if strand is valid
    if strand[i] != '+' && strand[i] != '-' {
      panic("NewGenes(): Invalid strand!")
    }
    index[names[i]] = i
  }
  return Genes{names, seqnames, tx, cds, strand, index, Meta{}}
}

func NewGenes(names, seqnames []string, txFrom, txTo, cdsFrom, cdsTo []int, strand []byte) Genes {
  n := len(names)
  if len(txFrom)   != n || len(txTo)   != n  ||
    (len(cdsFrom)  != n || len(cdsTo)  != n) ||
    (len(seqnames) != n || len(strand) != n) {
    panic("NewGenes(): invalid arguments!")
  }
  tx  := make([]Range, n)
  cds := make([]Range, n)
  for i := 0; i < n; i++ {
    // create tx and cds ranges
    tx [i] = NewRange( txFrom[i],  txTo[i])
    cds[i] = NewRange(cdsFrom[i], cdsTo[i])
  }
  return newGenes(names, seqnames, tx, cds, strand)
}

func (g *Genes) Clone() Genes {
  n := g.Length()
  names    := make([]string, n)
  seqnames := make([]string, n)
  cds      := make([]Range, n)
  tx       := make([]Range, n)
  strand   := make([]byte, n)
  copy(names,    g.Names)
  copy(seqnames, g.Seqnames)
  copy(cds,      g.Cds)
  copy(tx,       g.Tx)
  copy(strand,   g.Strand)
  result := newGenes(names, seqnames, cds, tx, strand)
  result.Meta = g.Meta.Clone()
  return result
}

/* -------------------------------------------------------------------------- */

// Number of genes in this object.
func (g Genes) Length() int {
  return len(g.Names)
}

// Returns the index of a gene.
func (g Genes) FindGene(name string) (int, bool) {
  i, ok := g.index[name]
  return i, ok
}

func (g1 *Genes) Append(g2 Genes) Genes {
  names    := append(g1.Names,    g2.Names...)
  seqnames := append(g1.Seqnames, g2.Seqnames...)
  cds      := append(g1.Cds,      g2.Cds...)
  tx       := append(g1.Tx,       g2.Tx...)
  strand   := append(g1.Strand,   g2.Strand...)

  result := newGenes(names, seqnames, cds, tx, strand)
  result.Meta = g1.Meta.Append(g2.Meta)

  return result
}

func (g *Genes) Remove(indices []int) Genes {
  if len(indices) == 0 {
    return g.Clone()
  }
  indices = RemoveDuplicatesInt(indices)
  sort.Ints(indices)

  n := g.Length()
  m := n - len(indices)
  // convert indices to subset indices
  idx := make([]int, m)
  for i, j, k := 0, 0, 0; i < g.Length(); i++ {
    for k < len(indices)-1 && i > indices[k] {
      k++
    }
    if i != indices[k] {
      idx[j] = i
      j++
    }
  }
  result := g.Subset(idx)
  result.Meta = g.Meta.Subset(idx)

  return result
}

// Returns a subset of the genes given by indices. This function may
// also be used to alter the order of genes.
func (g *Genes) Subset(indices []int) Genes {
  n := len(indices)
  names    := make([]string, n)
  seqnames := make([]string, n)
  txFrom   := make([]int, n)
  txTo     := make([]int, n)
  cdsFrom  := make([]int, n)
  cdsTo    := make([]int, n)
  strand   := make([]byte, n)

  for i := 0; i < n; i++ {
    names   [i] = g.Names   [indices[i]]
    seqnames[i] = g.Seqnames[indices[i]]
    txFrom  [i] = g.Tx      [indices[i]].From
    txTo    [i] = g.Tx      [indices[i]].To
    cdsFrom [i] = g.Cds     [indices[i]].From
    cdsTo   [i] = g.Cds     [indices[i]].To
    strand  [i] = g.Strand  [indices[i]]
  }
  result := NewGenes(names, seqnames, txFrom, txTo, cdsFrom, cdsTo, strand)
  result.Meta = g.Meta.Subset(indices)

  return result
}

// Returns the set of genes starting from row ifrom till ito (row ito
// is not included).
func (g *Genes) Slice(ifrom, ito int) Genes {
  n := ito-ifrom
  names    := make([]string, n)
  seqnames := make([]string, n)
  txFrom   := make([]int, n)
  txTo     := make([]int, n)
  cdsFrom  := make([]int, n)
  cdsTo    := make([]int, n)
  strand   := make([]byte, n)

  for i := ifrom; i < ito; i++ {
    names   [i-ifrom] = g.Names   [i]
    seqnames[i-ifrom] = g.Seqnames[i]
    txFrom  [i-ifrom] = g.Tx      [i].From
    txTo    [i-ifrom] = g.Tx      [i].To
    cdsFrom [i-ifrom] = g.Cds     [i].From
    cdsTo   [i-ifrom] = g.Cds     [i].To
    strand  [i-ifrom] = g.Strand  [i]
  }
  result := NewGenes(names, seqnames, txFrom, txTo, cdsFrom, cdsTo, strand)
  result.Meta = g.Meta.Slice(ifrom, ito)

  return result
}

// Sort rows using data from a meta column. The parameter name gives
// the name of the columns. If reverse is true, rows are sorted in
// descending order.
func (g *Genes) Sort(name string, reverse bool) (Genes, error) {
  j, err := g.sortedIndices(name, reverse)
  if err != nil {
    return Genes{}, err
  }
  return g.Subset(j), nil
}

/* convert to GRanges
 * -------------------------------------------------------------------------- */

// Export Genes object to GRanges. The transcription region is exported
// as the standard ranges field. Ranges of coding sequences and gene
// names are exported as meta data.
func (g *Genes) GRanges() *GRanges {
  n := g.Length()
  names   := make([]string, n)
  txFrom  := make([]int, n)
  txTo    := make([]int, n)
  cdsFrom := make([]int, n)
  cdsTo   := make([]int, n)
  for i := 0; i < n; i++ {
    names[i]   = g.Names[i]
    txFrom[i]  = g.Tx[i].From
    txTo[i]    = g.Tx[i].To
    cdsFrom[i] = g.Cds[i].From
    cdsTo[i]   = g.Cds[i].To
  }
  result := NewGRanges(g.Seqnames, txFrom, txTo, g.Strand)
  result.Meta = g.Meta.Clone()
  result.AddMeta("names",   names)
  result.AddMeta("cdsFrom", cdsFrom)
  result.AddMeta("cdsTo",   cdsTo)
  return &result
}

/* convert to string
 * -------------------------------------------------------------------------- */

func (genes Genes) PrettyPrint(n int) string {
  var buffer bytes.Buffer
  writer := bufio.NewWriter(&buffer)

  // pretty print meta data and create a scanner reading
  // the resulting string
  metaStr     := genes.Meta.PrettyPrint(n)
  metaReader  := strings.NewReader(metaStr)
  metaScanner := bufio.NewScanner(metaReader)

  // compute the width of a single cell
  updateMaxWidth := func(format string, widths []int, j int, args ...interface{}) {
    width, _ := fmt.Fprintf(ioutil.Discard, format, args...)
    if width > widths[j] {
      widths[j] = width
    }
  }
  // compute widths of all cells in row i
  updateMaxWidths := func(i int, widths []int) {
    updateMaxWidth("%d", widths, 0, i+1)
    updateMaxWidth("%s", widths, 1, genes.Names[i])
    updateMaxWidth("%s", widths, 2, genes.Seqnames[i])
    updateMaxWidth("%d", widths, 3, genes.Tx[i].From)
    updateMaxWidth("%d", widths, 4, genes.Tx[i].To)
    updateMaxWidth("%d", widths, 5, genes.Cds[i].From)
    updateMaxWidth("%d", widths, 6, genes.Cds[i].To)
    updateMaxWidth("%c", widths, 7, genes.Strand[i])
  }
  printMetaRow := func(writer io.Writer) {
    if genes.MetaLength() != 0 {
      fmt.Fprintf(writer, " | ")
      metaScanner.Scan()
      fmt.Fprintf(writer, "%s", metaScanner.Text())
    }
  }
  printHeader := func(writer io.Writer, format string) {
    fmt.Fprintf(writer, format,
      "", "names", "seqnames", "transcripts", "cds", "strand")
    printMetaRow(writer)
    fmt.Fprintf(writer, "\n")
  }
  printRow := func(writer io.Writer, format string, i int) {
    if i != 0 {
      fmt.Fprintf(writer, "\n")
    }
    fmt.Fprintf(writer, format,
      i+1,
      genes.Names[i],
      genes.Seqnames[i],
      genes.Tx[i].From,
      genes.Tx[i].To,
      genes.Cds[i].From,
      genes.Cds[i].To,
      genes.Strand[i])
    printMetaRow(writer)
  }
  applyRows := func(f1 func(i int), f2 func()) {
    if genes.Length() <= n+1 {
      // apply to all entries
      for i := 0; i < genes.Length(); i++ { f1(i) }
    } else {
      // apply to first n/2 rows
      for i := 0; i < n/2; i++ { f1(i) }
      // between first and last n/2 rows
      f2()
      // apply to last n/2 rows
      for i := genes.Length() - n/2; i < genes.Length(); i++ { f1(i) }
    }
  }
  // maximum column widths
  widths := []int{1, 5, 8, 1, 1, 1, 1, 6}
  // determine column widths
  applyRows(func(i int) { updateMaxWidths(i, widths) }, func() {})
  // generate format strings
  formatRow    := fmt.Sprintf("%%%dd %%%ds %%%ds [%%%dd, %%%dd) [%%%dd, %%%dd) %%%dc",
    widths[0], widths[1], widths[2], widths[3], widths[4], widths[5], widths[6], widths[7])
  formatHeader := fmt.Sprintf("%%%ds %%%ds %%%ds %%%ds %%%ds %%%ds",
    widths[0], widths[1], widths[2], widths[3]+widths[4]+4, widths[5]+widths[6]+4, widths[7])
  // pring header
  printHeader(writer, formatHeader)
  // print rows
  applyRows(
    func(i int) {
      printRow(writer, formatRow, i)
    },
    func() {
      fmt.Fprintf(writer, "\n")
      fmt.Fprintf(writer, formatHeader, "", "...", "...", "...", "...", "")
      printMetaRow(writer)
    })
  writer.Flush()

  return buffer.String()
}

func (genes Genes) String() string {
  return genes.PrettyPrint(10)
}

/* i/o
 * -------------------------------------------------------------------------- */

// Import genes from UCSC text files. The format is a whitespace separated
// table with columns: Name, Seqname, TranscriptStart, TranscriptEnd,
// CodingStart, CodingEnd, and Strand.
func ReadUCSCGenes(filename string) Genes {
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
  names    := []string{}
  seqnames := []string{}
  txFrom   := []int{}
  txTo     := []int{}
  cdsFrom  := []int{}
  cdsTo    := []int{}
  strand   := []byte{}

  for scanner.Scan() {
    err = scanner.Err()
    if err != nil {
      panic(err)
    }
    fields := strings.Fields(scanner.Text())
    if len(fields) == 0 {
      continue
    }
    if len(fields) != 7 {
      panic("File must have seven columns!")
    }
    t1, e := strconv.ParseInt(fields[3], 10, 64); Check(e)
    t2, e := strconv.ParseInt(fields[4], 10, 64); Check(e)
    t3, e := strconv.ParseInt(fields[5], 10, 64); Check(e)
    t4, e := strconv.ParseInt(fields[6], 10, 64); Check(e)
    names    = append(names,    fields[0])
    seqnames = append(seqnames, fields[1])
    txFrom   = append(txFrom,   int(t1))
    txTo     = append(txTo,     int(t2))
    cdsFrom  = append(cdsFrom,  int(t3))
    cdsTo    = append(cdsTo,    int(t4))
    strand   = append(strand,   fields[2][0])
  }
  return NewGenes(names, seqnames, txFrom, txTo, cdsFrom, cdsTo, strand)
}

func (genes Genes) WriteTable(filename string, header, compress bool) {
  var buffer bytes.Buffer

  w := bufio.NewWriter(&buffer)

  // print header
  if header {
    fmt.Fprintf(w, "%16s %10s %6s %10s %10s %10s %10s",
      "names", "seqnames", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd")
    genes.Meta.WriteTableRow(w, -1)
    fmt.Fprintf(w, "\n")
  }
  // print data
  for i := 0; i < genes.Length(); i++ {
    fmt.Fprintf(w,  "%16s", genes.Names[i])
    fmt.Fprintf(w, " %10s", genes.Seqnames[i])
    fmt.Fprintf(w, " %6c",  genes.Strand[i])
    fmt.Fprintf(w, " %10d", genes.Tx[i].From)
    fmt.Fprintf(w, " %10d", genes.Tx[i].To)
    fmt.Fprintf(w, " %10d", genes.Cds[i].From)
    fmt.Fprintf(w, " %10d", genes.Cds[i].To)
    genes.Meta.WriteTableRow(w, i)
    fmt.Fprintf(w, "\n")
  }
  w.Flush()
  writeFile(filename, &buffer, compress)
}


/* import genes from ucsc
 * -------------------------------------------------------------------------- */

func ImportGenesFromUCSC(genome, table string) Genes {
  /* variables for storing a single database row */
  var i_name, i_seqname, i_strand string
  var i_txFrom, i_txTo, i_cdsFrom, i_cdsTo int

  names    := []string{}
  seqnames := []string{}
  txFrom   := []int{}
  txTo     := []int{}
  cdsFrom  := []int{}
  cdsTo    := []int{}
  strand   := []byte{}

  /* open connection */
  db, err := sql.Open("mysql",
		fmt.Sprintf("genome@tcp(genome-mysql.cse.ucsc.edu:3306)/%s", genome))
	if err != nil {
		log.Fatal(err)
	}
	defer db.Close()

  err = db.Ping()
  if err != nil {
    log.Fatal(err)
  }

  /* receive data */
  rows, err := db.Query(
    fmt.Sprintf("SELECT name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd FROM %s", table))
  if err != nil {
    log.Fatal(err)
  }
  defer rows.Close()
  for rows.Next() {
    err := rows.Scan(&i_name, &i_seqname, &i_strand, &i_txFrom, &i_txTo, &i_cdsFrom, &i_cdsTo)
    if err != nil {
      log.Fatal(err)
    }
    names    = append(names,    i_name)
    seqnames = append(seqnames, i_seqname)
    txFrom   = append(txFrom,   i_txFrom)
    txTo     = append(txTo,     i_txTo)
    cdsFrom  = append(cdsFrom,  i_cdsFrom)
    cdsTo    = append(cdsTo,    i_cdsTo)
    strand   = append(strand,   i_strand[0])
  }
  return NewGenes(names, seqnames, txFrom, txTo, cdsFrom, cdsTo, strand)
}
