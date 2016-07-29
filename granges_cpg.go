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
import "log"
import "database/sql"

import _ "github.com/go-sql-driver/mysql"

/* import data from ucsc
 * -------------------------------------------------------------------------- */

func ImportCpGIslandsFromUCSC(genome string) GRanges {
  /* variables for storing a single database row */
  var i_seqname string
  var i_from, i_to, i_length, i_cpgNum, i_gcNum int
  var i_perCpg, i_perGc, i_obsExp float64

  seqnames := []string{}
  from     := []int{}
  to       := []int{}
  length   := []int{}
  cpgNum   := []int{}
  gcNum    := []int{}
  perCpg   := []float64{}
  perGc    := []float64{}
  obsExp   := []float64{}

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
    fmt.Sprintf("SELECT chrom, chromStart, chromEnd, length, cpgNum, gcNum, perCpg, perGc, obsExp FROM %s", "cpgIslandExt"))
  if err != nil {
    log.Fatal(err)
  }
  defer rows.Close()
  for rows.Next() {
    err := rows.Scan(&i_seqname, &i_from, &i_to, &i_length, &i_cpgNum, &i_gcNum, &i_perCpg, &i_perGc, &i_obsExp)
    if err != nil {
      log.Fatal(err)
    }
    seqnames = append(seqnames, i_seqname)
    from     = append(from,     i_from)
    to       = append(to,       i_to)
    length   = append(length,   i_length)
    cpgNum   = append(cpgNum,   i_cpgNum)
    gcNum    = append(gcNum,    i_gcNum)
    perCpg   = append(perCpg,   i_perCpg)
    perGc    = append(perGc,    i_perGc)
    obsExp   = append(obsExp,   i_obsExp)
  }
  r := NewGRanges(seqnames, from, to, []byte{})
  r.AddMeta("length", length)
  r.AddMeta("cpgNum", cpgNum)
  r.AddMeta("gcNum",  gcNum)
  r.AddMeta("perCpg", perCpg)
  r.AddMeta("perGc",  perGc)
  r.AddMeta("obsExp", obsExp)

  return r
}

func ReadCpGIslandsFromTable(filename string) GRanges {
  cpg := GRanges{}
  cpg.ReadTable(filename,
    []string{"length", "cpgNum", "gcNum", "perCpg", "perGc", "obsExp"},
    []string{"[]int", "[]int", "[]int", "[]float64", "[]float64", "[]float64"})
  return cpg
}
