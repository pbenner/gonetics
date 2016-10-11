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
import "fmt"
import "compress/gzip"
import "database/sql"
import "os"
import "strconv"
import "strings"

import _ "github.com/go-sql-driver/mysql"

/* import genes from ucsc
 * -------------------------------------------------------------------------- */

// Import genes from UCSC text files. The format is a whitespace separated
// table with columns: Name, Seqname, TranscriptStart, TranscriptEnd,
// CodingStart, CodingEnd, and Strand.
func ReadUCSCGenes(filename string) (Genes, error) {
  var genes Genes
  var scanner *bufio.Scanner
  // open file
  f, err := os.Open(filename)
  if err != nil {
    return genes, err
  }
  defer f.Close()
  // check if file is gzipped
  if isGzip(filename) {
    g, err := gzip.NewReader(f)
    if err != nil {
      return genes, err
    }
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
      return genes, err
    }
    fields := strings.Fields(scanner.Text())
    if len(fields) == 0 {
      continue
    }
    if len(fields) != 7 {
      return genes, fmt.Errorf("file must have seven columns")
    }
    t1, e := strconv.ParseInt(fields[3], 10, 64)
    if e != nil {
      return genes, e
    }
    t2, e := strconv.ParseInt(fields[4], 10, 64)
    if e != nil {
      return genes, e
    }
    t3, e := strconv.ParseInt(fields[5], 10, 64)
    if e != nil {
      return genes, e
    }
    t4, e := strconv.ParseInt(fields[6], 10, 64)
    if e != nil {
      return genes, e
    }
    names    = append(names,    fields[0])
    seqnames = append(seqnames, fields[1])
    txFrom   = append(txFrom,   int(t1))
    txTo     = append(txTo,     int(t2))
    cdsFrom  = append(cdsFrom,  int(t3))
    cdsTo    = append(cdsTo,    int(t4))
    strand   = append(strand,   fields[2][0])
  }
  return NewGenes(names, seqnames, txFrom, txTo, cdsFrom, cdsTo, strand), nil
}

func ImportGenesFromUCSC(genome, table string) (Genes, error) {
  genes := Genes{}
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
		return genes, err
	}
	defer db.Close()

  err = db.Ping()
  if err != nil {
    return genes,err
  }

  /* receive data */
  rows, err := db.Query(
    fmt.Sprintf("SELECT name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd FROM %s", table))
  if err != nil {
    return genes, err
  }
  defer rows.Close()
  for rows.Next() {
    err := rows.Scan(&i_name, &i_seqname, &i_strand, &i_txFrom, &i_txTo, &i_cdsFrom, &i_cdsTo)
    if err != nil {
      return genes, err
    }
    names    = append(names,    i_name)
    seqnames = append(seqnames, i_seqname)
    txFrom   = append(txFrom,   i_txFrom)
    txTo     = append(txTo,     i_txTo)
    cdsFrom  = append(cdsFrom,  i_cdsFrom)
    cdsTo    = append(cdsTo,    i_cdsTo)
    strand   = append(strand,   i_strand[0])
  }
  return NewGenes(names, seqnames, txFrom, txTo, cdsFrom, cdsTo, strand), nil
}
