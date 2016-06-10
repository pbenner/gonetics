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
import "errors"
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

type GRanges struct {
  Seqnames   []string
  Ranges     []Range
  Strand     []byte
  Meta
}

/* constructors
 * -------------------------------------------------------------------------- */

func NewGRanges(seqnames []string, from, to []int, strand []byte) GRanges {
  n := len(seqnames)
  if len(  from) != n || len(    to) != n ||
    (len(strand) != 0 && len(strand) != n) {
    panic("NewGRanges(): invalid arguments!")
  }
  if len(strand) == 0 {
    strand = make([]byte, n)
    for i := 0; i < n; i++ {
      strand[i] = '*'
    }
  }
  ranges := make([]Range, n)
  for i := 0; i < n; i++ {
    // create range
    ranges[i] = NewRange(from[i], to[i])
    // check if strand is valid
    if strand[i] != '+' && strand[i] != '-' && strand[i] != '*' {
      panic("NewGRanges(): Invalid strand!")
    }
  }
  return GRanges{seqnames, ranges, strand, Meta{}}
}

func NewEmptyGRanges(n int) GRanges {
  seqnames := make([]string, n)
  ranges   := make([]Range, n)
  strand   := make([]byte, n)
  for i := 0; i < n; i++ {
    strand[i] = '*'
  }
  return GRanges{seqnames, ranges, strand, Meta{}}
}

func (r *GRanges) Clone() GRanges {
  result := GRanges{}
  n := r.Length()
  result.Seqnames = make([]string, n)
  result.Ranges   = make([]Range, n)
  result.Strand   = make([]byte, n)
  copy(result.Seqnames, r.Seqnames)
  copy(result.Ranges,   r.Ranges)
  copy(result.Strand,   r.Strand)
  result.Meta = r.Meta.Clone()
  return result
}

/* -------------------------------------------------------------------------- */

func (r GRanges) Length() int {
  return len(r.Ranges)
}

func (r1 GRanges) Append(r2 GRanges) GRanges {
  result := GRanges{}

  result.Seqnames = append(r1.Seqnames, r2.Seqnames...)
  result.Ranges   = append(r1.Ranges,   r2.Ranges...)
  result.Strand   = append(r1.Strand,   r2.Strand...)

  result.Meta = r1.Meta.Append(r2.Meta)

  return result
}

func (r GRanges) Remove(indices []int) GRanges {
  if len(indices) == 0 {
    return r.Clone()
  }
  indices = RemoveDuplicatesInt(indices)
  sort.Ints(indices)

  n := r.Length()
  m := n - len(indices)
  // convert indices to subset indices
  idx := make([]int, m)
  for i, j, k := 0, 0, 0; i < r.Length(); i++ {
    for k < len(indices)-1 && i > indices[k] {
      k++
    }
    if i != indices[k] {
      idx[j] = i
      j++
    }
  }
  result := r.Subset(idx)
  result.Meta = r.Meta.Subset(idx)

  return result
}

func (r GRanges) RemoveOverlapsWith(subject GRanges) GRanges {
  queryHits, _ := FindOverlaps(r, subject)
  return r.Remove(queryHits)
}

func (r GRanges) Subset(indices []int) GRanges {
  n := len(indices)
  seqnames := make([]string, n)
  from     := make([]int, n)
  to       := make([]int, n)
  strand   := make([]byte, n)

  for i := 0; i < n; i++ {
    seqnames[i] = r.Seqnames[indices[i]]
    from    [i] = r.Ranges  [indices[i]].From
    to      [i] = r.Ranges  [indices[i]].To
    strand  [i] = r.Strand  [indices[i]]
  }
  result := NewGRanges(seqnames, from, to, strand)
  result.Meta = r.Meta.Subset(indices)

  return result
}

func (r GRanges) Slice(ifrom, ito int) GRanges {
  n := ito-ifrom
  seqnames := make([]string, n)
  from     := make([]int, n)
  to       := make([]int, n)
  strand   := make([]byte, n)

  for i := ifrom; i < ito; i++ {
    seqnames[i-ifrom] = r.Seqnames[i]
    from    [i-ifrom] = r.Ranges  [i].From
    to      [i-ifrom] = r.Ranges  [i].To
    strand  [i-ifrom] = r.Strand  [i]
  }
  result := NewGRanges(seqnames, from, to, strand)
  result.Meta = r.Meta.Slice(ifrom, ito)

  return result
}

func (r GRanges) Sort(name string, reverse bool) (GRanges, error) {
  j, err := r.sortedIndices(name, reverse)
  if err != nil {
    return GRanges{}, err
  }
  return r.Subset(j), nil
}

// Remove all entries that are not in the given genome.
func (r GRanges) FilterGenome(genome Genome) GRanges {
  idx      := []int{}
  seqnames := make(map[string]int)
  for i := 0; i < genome.Length(); i++ {
    seqnames[genome.Seqnames[i]] = genome.Lengths[i]
  }
  for i := 0; i < r.Length(); i++ {
    length, ok := seqnames[r.Seqnames[i]]
    if ok && r.Ranges[i].To <= length {
      continue
    }
    idx = append(idx, i)
  }
  return r.Remove(idx)
}

// Add data from a track to the GRanges object. The data will be
// contained in a meta-data column with the same name as the track.
// It is required that each range has the same length.
func (r *GRanges) ImportTrack(track Track, revNegStrand bool) *GRanges {
  n := r.Length()
  m := -1
  data    := make([][]float64, n)
  binsize := make([]int, n)
  // fill matrix
  for i := 0; i < n; i++ {
    from := r.Ranges[i].From
    to   := r.Ranges[i].To
    seq  := r.Seqnames[i]
    if m == -1 {
      m = DivIntUp(to - from, track.Binsize)
    } else if m != DivIntUp(to - from, track.Binsize) {
      panic("varying window sizes are not allowed")
    }
    // all rows are using the same binsize
    binsize[i] = track.Binsize
    // loop over window
    data[i] = make([]float64, m)
    if r.Strand[i] == '+' || revNegStrand == false {
      for j, k := 0, from; k < to; k, j = k+track.Binsize, j+1 {
        value, err := track.At(seq, k)
        if err == nil {
          data[i][j] = value
        }
      }
    } else if r.Strand[i] == '-' {
      for j, k := 0, to-1; k >= from; k, j = k-track.Binsize, j+1 {
        value, err := track.At(seq, k)
        if err == nil {
          data[i][j] = value
        }
      }
    } else {
      panic("range has no strand information")
    }
  }
  r.AddMeta("binsize", binsize)
  r.AddMeta(track.Name, data)
  return r
}

/* convert to string
 * -------------------------------------------------------------------------- */

func (granges GRanges) PrettyPrint(n int) string {
  var buffer bytes.Buffer
  writer := bufio.NewWriter(&buffer)

  // pretty print meta data and create a scanner reading
  // the resulting string
  metaStr     := granges.Meta.PrettyPrint(n)
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
    updateMaxWidth("%s", widths, 1, granges.Seqnames[i])
    updateMaxWidth("%d", widths, 2, granges.Ranges[i].From)
    updateMaxWidth("%d", widths, 3, granges.Ranges[i].To)
    updateMaxWidth("%c", widths, 4, granges.Strand[i])
  }
  printMetaRow := func(writer io.Writer) {
    if granges.MetaLength() != 0 {
      fmt.Fprintf(writer, " | ")
      metaScanner.Scan()
      fmt.Fprintf(writer, "%s", metaScanner.Text())
    }
  }
  printHeader := func(writer io.Writer, format string) {
    fmt.Fprintf(writer, format,
      "", "seqnames", "ranges", "strand")
    printMetaRow(writer)
  }
  printRow := func(writer io.Writer, format string, i int) {
    fmt.Fprintf(writer, "\n")
    fmt.Fprintf(writer, format,
      i+1,
      granges.Seqnames[i],
      granges.Ranges[i].From,
      granges.Ranges[i].To,
      granges.Strand[i])
    printMetaRow(writer)
  }
  applyRows := func(f1 func(i int), f2 func()) {
    if granges.Length() <= n+1 {
      // apply to all entries
      for i := 0; i < granges.Length(); i++ { f1(i) }
    } else {
      // apply to first n/2 rows
      for i := 0; i < n/2; i++ { f1(i) }
      // between first and last n/2 rows
      f2()
      // apply to last n/2 rows
      for i := granges.Length() - n/2; i < granges.Length(); i++ { f1(i) }
    }
  }
  // maximum column widths
  widths := []int{1, 8, 1, 1, 6}
  // determine column widths
  applyRows(func(i int) { updateMaxWidths(i, widths) }, func() {})
  // generate format strings
  formatRow    := fmt.Sprintf("%%%dd %%%ds [%%%dd, %%%dd) %%%dc",
    widths[0], widths[1], widths[2], widths[3], widths[4])
  formatHeader := fmt.Sprintf("%%%ds %%%ds %%%ds %%%ds",
    widths[0], widths[1], widths[2]+widths[3]+4, widths[4])
  // pring header
  printHeader(writer, formatHeader)
  // print rows
  applyRows(
    func(i int) {
      printRow(writer, formatRow, i)
    },
    func() {
      fmt.Fprintf(writer, "\n")
      fmt.Fprintf(writer, formatHeader, "", "...", "...", "")
      printMetaRow(writer)
    })
  writer.Flush()

  return buffer.String()
}

func (granges GRanges) String() string {
  return granges.PrettyPrint(10)
}

/* i/o
 * -------------------------------------------------------------------------- */

// Export GRanges as GTF file. Required GTF fields should be provided
// as meta columns named sources, features, scores, and frames. All other
// meta columns are exported as optional fields.
func (granges GRanges) WriteGTF(filename string) {
  f, err := os.Create(filename); Check(err)
  defer f.Close()

  w := bufio.NewWriter(f)
  defer w.Flush()

  sources  := granges.GetMetaStr("sources")
  features := granges.GetMetaStr("features")
  scores   := granges.GetMetaInt("scores")
  frames   := granges.GetMetaInt("frames")

  for i := 0; i < granges.Length(); i++ {
    fmt.Fprintf(w, "%s", granges.Seqnames[i])
    if len(sources) > 0 {
      fmt.Fprintf(w, "\t%s", sources[i])
    } else {
      fmt.Fprintf(w, "\t%s", ".")
    }
    if len(features) > 0 {
      fmt.Fprintf(w, "\t%s", features[i])
    } else {
      fmt.Fprintf(w, "\t%s", ".")
    }
    fmt.Fprintf(w, "\t%d", granges.Ranges[i].From)
    fmt.Fprintf(w, "\t%d", granges.Ranges[i].To)
    if len(scores) > 0 {
      fmt.Fprintf(w, "\t%d", scores[i])
    } else {
      fmt.Fprintf(w, "\t%s", ".")
    }
    if len(granges.Strand) > 0 && granges.Strand[i] != '*' {
      fmt.Fprintf(w, "\t%c", granges.Strand[i])
    } else {
      fmt.Fprintf(w, "\t%c", '.')
    }
    if len(frames) > 0 {
      fmt.Fprintf(w, "\t%d", frames[i])
    } else {
      fmt.Fprintf(w, "\t%s", ".")
    }

    if granges.MetaLength() != 0 {
      printedTab := false
      for k := 0; k < granges.MetaLength(); k++ {
        if granges.MetaName[k] == "sources" {
          continue
        }
        if granges.MetaName[k] == "features" {
          continue
        }
        if granges.MetaName[k] == "scores" {
          continue
        }
        if granges.MetaName[k] == "frames" {
          continue
        }
        if printedTab {
          w.WriteString(" ")
        } else {
          w.WriteString("\t")
          printedTab = true
        }
        // print name of the meta data
        fmt.Fprintf(w, "%s ", granges.MetaName[k])
        // print data
        switch v := granges.MetaData[k].(type) {
        case []string : fmt.Fprintf(w, "\"%s\"", v[i])
        case []float64: fmt.Fprintf(w, "\"%f\"", v[i])
        case []int    : fmt.Fprintf(w, "\"%d\"", v[i])
        case [][]string:
          w.WriteString("\"")
          for j := 0; j < len(v[i]); j++ {
            if j != 0 {
              w.WriteString(" ")
            }
            fmt.Fprintf(w, "%s", v[i][j])
          }
          w.WriteString("\"")
        case [][]float64:
          w.WriteString("\"")
          for j := 0; j < len(v[i]); j++ {
            if j != 0 {
              w.WriteString(" ")
            }
            fmt.Fprintf(w, "%f", v[i][j])
          }
          w.WriteString("\"")
        case [][]int:
          w.WriteString("\"")
          for j := 0; j < len(v[i]); j++ {
            if j != 0 {
              w.WriteString(" ")
            }
            fmt.Fprintf(w, "%d", v[i][j])
          }
          w.WriteString("\"")
        }
        w.WriteString(";")
      }
    }
    w.WriteString("\n")
  }
}

// Export GRanges as a table. The first line contains the header
// of the table.
func (granges GRanges) WriteTable(filename string, header, strand, compress bool) {
  var buffer bytes.Buffer

  w := bufio.NewWriter(&buffer)

  // print header
  if header {
    if strand {
      fmt.Fprintf(w, "%10s %10s %10s %6s", "seqnames", "from", "to", "strand")
    } else {
      fmt.Fprintf(w, "%10s %10s %10s", "seqnames", "from", "to")
    }
    granges.Meta.WriteTableRow(w, -1)
    fmt.Fprintf(w, "\n")
  }
  // print data
  for i := 0; i < granges.Length(); i++ {
    fmt.Fprintf(w,  "%10s", granges.Seqnames[i])
    fmt.Fprintf(w, " %10d", granges.Ranges[i].From)
    fmt.Fprintf(w, " %10d", granges.Ranges[i].To)
    if strand {
      if len(granges.Strand) > 0 {
        fmt.Fprintf(w, " %6c", granges.Strand[i])
      } else {
        fmt.Fprintf(w, " %6c", '*')
      }
    }
    granges.Meta.WriteTableRow(w, i)
    fmt.Fprintf(w, "\n")
  }
  w.Flush()
  writeFile(filename, &buffer, compress)
}

// Export GRanges object as bed file with three columns.
func (granges GRanges) WriteBed3(filename string, compress bool) {
  var buffer bytes.Buffer

  w := bufio.NewWriter(&buffer)

  // print data
  for i := 0; i < granges.Length(); i++ {
    fmt.Fprintf(w,   "%s", granges.Seqnames[i])
    fmt.Fprintf(w, "\t%d", granges.Ranges[i].From)
    fmt.Fprintf(w, "\t%d", granges.Ranges[i].To)
    fmt.Fprintf(w, "\n")
  }
  w.Flush()
  writeFile(filename, &buffer, compress)
}

func (granges GRanges) WriteBed6(filename string, compress bool) {
  var buffer bytes.Buffer

  w := bufio.NewWriter(&buffer)

  name  := granges.GetMetaStr  ("name")
  score := granges.GetMetaFloat("score")

  for i := 0; i < granges.Length(); i++ {
    fmt.Fprintf(w,   "%s", granges.Seqnames[i])
    fmt.Fprintf(w, "\t%d", granges.Ranges[i].From)
    fmt.Fprintf(w, "\t%d", granges.Ranges[i].To)
    if len(name) > 0 {
      fmt.Fprintf(w, "\t%s", name[i])
    } else {
      fmt.Fprintf(w, "\t%s", ".")
    }
    if len(score) > 0 {
      fmt.Fprintf(w, "\t%f", score[i])
    } else {
      fmt.Fprintf(w, "\t%s", ".")
    }
    if len(granges.Strand) > 0 {
      fmt.Fprintf(w, "\t%c", granges.Strand[i])
    } else {
      fmt.Fprintf(w, "\t%c", '*')
    }
    fmt.Fprintf(w, "\n")
  }
  w.Flush()
  writeFile(filename, &buffer, compress)
}

func ReadGRangesFromTable(filename string, names, types []string) (GRanges, error) {
  result    := GRanges{}
  hasStrand := false

  var scanner *bufio.Scanner
  // open file
  f, err := os.Open(filename)
  if err != nil {
    return result, err
  }
  defer f.Close()

  // check if file is gzipped
  if IsGzip(filename) {
    g, err := gzip.NewReader(f)
    if err != nil {
      return result, err
    }
    defer g.Close()
    scanner = bufio.NewScanner(g)
  } else {
    scanner = bufio.NewScanner(f)
  }

  // scan header
  if scanner.Scan() {
    fields := strings.Fields(scanner.Text())
    // get number of columns
    if len(fields) < 4 {
      return result, errors.New("invalid table")
    }
    if fields[0] != "seqnames" || fields[1] != "from" || fields[2] != "to" {
      return result, errors.New("invalid table")
    }
    if fields[3] == "strand" {
      hasStrand = true
    }
  }
  for scanner.Scan() {
    fields := strings.Fields(scanner.Text())
    if len(fields) == 0 {
      continue
    }
    if len(fields) < 4 {
      return result, errors.New("invalid table")
    }
    v1, err := strconv.ParseInt(fields[1], 10, 64)
    if err != nil {
      return result, err
    }
    v2, err := strconv.ParseInt(fields[2], 10, 64)
    if err != nil {
      return result, err
    }
    result.Seqnames = append(result.Seqnames, fields[0])
    result.Ranges   = append(result.Ranges,   NewRange(int(v1), int(v2)))
    if hasStrand {
      result.Strand = append(result.Strand,   fields[3][0])
    } else {
      result.Strand = append(result.Strand,   '*')
    }
  }
  result.Meta, err = ReadMetaFromTable(filename, names, types)

  return result, err
}

// Import reads from a Bed file. The file is required to have at least
// six columns.
func ReadBedReads(filename string) GRanges {
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
  seqnames := []string{}
  from     := []int{}
  to       := []int{}
  strand   := []byte{}

  for scanner.Scan() {
    fields := strings.Fields(scanner.Text())
    if len(fields) == 0 {
      continue
    }
    if len(fields) < 6 {
      panic("Bed file must have at least six columns!")
    }
    t1, e1 := strconv.ParseInt(fields[1], 10, 64)
    t2, e2 := strconv.ParseInt(fields[2], 10, 64)
    Check(e1)
    Check(e2)
    seqnames = append(seqnames, fields[0])
    from     = append(from,     int(t1))
    to       = append(to,       int(t2))
    strand   = append(strand,   fields[5][0])
  }
  return NewGRanges(seqnames, from, to, strand)
}

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
  cpg, _ := ReadGRangesFromTable(filename,
    []string{"length", "cpgNum", "gcNum", "perCpg", "perGc", "obsExp"},
    []string{"[]int", "[]int", "[]int", "[]float64", "[]float64", "[]float64"})
  return cpg
}
