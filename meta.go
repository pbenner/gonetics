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
import "bytes"
import "bufio"
import "compress/gzip"
import "errors"
import "io"
import "os"
import "sort"
import "strconv"
import "strings"

import . "github.com/pbenner/pshape/Utility"

/* -------------------------------------------------------------------------- */

type Meta struct {
  MetaName []string
  MetaData []interface{}
  rows int
}

/* constructors
 * -------------------------------------------------------------------------- */

func NewMeta(names []string, data []interface{}) Meta {
  meta := Meta{}
  if len(names) != len(data) {
    panic("NewMeta(): invalid parameters!")
  }
  for i := 0; i < len(names); i++ {
    meta.AddMeta(names[i], data[i])
  }
  return meta
}

// Deep copy the Meta object.
func (m *Meta) Clone() Meta {
  // empty meta data
  result := NewMeta([]string{}, []interface{}{})
  // loop over meta data and clone elements
  for i := 0; i < m.MetaLength(); i++ {
    switch v := m.MetaData[i].(type) {
    case [][]string:
      r := make([][]string, len(v))
      for j := 0; j < len(v); j++ {
        r[j] = make([]string, len(v[j]))
        copy(r[j], v[j])
      }
      result.AddMeta(m.MetaName[i], r)
    case []string:
      r := make([]string, len(v))
      copy(r, v)
      result.AddMeta(m.MetaName[i], r)
    case [][]float64:
      r := make([][]float64, len(v))
      for j := 0; j < len(v); j++ {
        r[j] = make([]float64, len(v[j]))
        copy(r[j], v[j])
      }
      result.AddMeta(m.MetaName[i], r)
    case []float64:
      r := make([]float64, len(v))
      copy(r, v)
      result.AddMeta(m.MetaName[i], r)
    case [][]int:
      r := make([][]int, len(v))
      for j := 0; j < len(v); j++ {
        r[j] = make([]int, len(v[j]))
        copy(r[j], v[j])
      }
      result.AddMeta(m.MetaName[i], r)
    case []int:
      r := make([]int, len(v))
      copy(r, v)
      result.AddMeta(m.MetaName[i], r)
    default: panic("AddMeta(): invalid type!")
    }
  }
  return result
}

/* -------------------------------------------------------------------------- */

// Returns the number of rows.
func (m *Meta) Length() int {
  return m.rows
}

// Returns the number of columns.
func (m *Meta) MetaLength() int {
  return len(m.MetaName)
}

func (m *Meta) AddMeta(name string, meta interface{}) {
  if m.MetaLength() > 0 {
    // this is not the first column added; check length
    switch v := meta.(type) {
    case [][]string:  if len(v) != m.rows { goto lengthErr }
    case   []string:  if len(v) != m.rows { goto lengthErr }
    case [][]float64: if len(v) != m.rows { goto lengthErr }
    case   []float64: if len(v) != m.rows { goto lengthErr }
    case [][]int:     if len(v) != m.rows { goto lengthErr }
    case   []int:     if len(v) != m.rows { goto lengthErr }
    default: panic("AddMeta(): invalid type!")
    }
  } else {
    // this is the first column, determine length
    switch v := meta.(type) {
    case [][]string:  m.rows = len(v)
    case   []string:  m.rows = len(v)
    case [][]float64: m.rows = len(v)
    case   []float64: m.rows = len(v)
    case [][]int:     m.rows = len(v)
    case   []int:     m.rows = len(v)
    default: panic("AddMeta(): invalid type!")
    }
  }
  m.MetaData = append(m.MetaData, meta)
  m.MetaName = append(m.MetaName, name)
  return
lengthErr:
  panic("AddMeta(): column has invalid length!")
}

func (m *Meta) DeleteMeta(name string) {
  for i := 0; i < m.MetaLength(); i++ {
    if m.MetaName[i] == name {
      m.MetaName = append(m.MetaName[:i], m.MetaName[i+1:]...)
      m.MetaData = append(m.MetaData[:i], m.MetaData[i+1:]...)
    }
  }
}

func (m *Meta) RenameMeta(nameOld, nameNew string) {
  for i := 0; i < m.MetaLength(); i++ {
    if m.MetaName[i] == nameOld {
      m.MetaName[i] = nameNew
    }
  }
}

func (m *Meta) GetMeta(name string) interface{} {
  for i := 0; i < m.MetaLength(); i++ {
    if m.MetaName[i] == name {
      return m.MetaData[i]
    }
  }
  return nil
}

func (m *Meta) GetMetaStr(name string) []string {
  r := m.GetMeta(name)
  if r != nil {
    return r.([]string)
  }
  return []string{}
}

func (m *Meta) GetMetaFloat(name string) []float64 {
  r := m.GetMeta(name)
  if r != nil {
    return r.([]float64)
  }
  return []float64{}
}

func (m *Meta) GetMetaInt(name string) []int {
  r := m.GetMeta(name)
  if r != nil {
    return r.([]int)
  }
  return []int{}
}

func (meta1 *Meta) Append(meta2 Meta) Meta {
  result := Meta{}

  // clone data so we do not have to deep copy
  // two-dimensional slices
  m1 := meta1.Clone()
  m2 := meta2.Clone()

  for j := 0; j < m1.MetaLength(); j++ {
    var t interface{}

    name := m1.MetaName[j]
    dat1 := m1.MetaData[j]
    dat2 := m2.GetMeta(name)

    switch v := dat1.(type) {
    case [][]string:  t = append(v, dat2.([][]string)...)
    case [][]int:     t = append(v, dat2.([][]int)...)
    case [][]float64: t = append(v, dat2.([][]float64)...)
    case   []string:  t = append(v, dat2.(  []string)...)
    case   []int:     t = append(v, dat2.(  []int)...)
    case   []float64: t = append(v, dat2.(  []float64)...)
    }
    result.AddMeta(name, t)
  }
  return result
}

func (meta *Meta) Remove(indices []int) Meta {
  if len(indices) == 0 {
    return meta.Clone()
  }
  indices = RemoveDuplicatesInt(indices)
  sort.Ints(indices)

  n := meta.Length()
  m := n - len(indices)
  // convert indices to subset indices
  idx := make([]int, m)
  for i, j, k := 0, 0, 0; i < meta.Length(); i++ {
    for k < len(indices)-1 && i > indices[k] {
      k++
    }
    if i != indices[k] {
      idx[j] = i
      j++
    }
  }
  return meta.Subset(idx)
}

// Return a new Meta object with a subset of the rows from
// this object.
func (meta *Meta) Subset(indices []int) Meta {
  n := len(indices)
  m := meta.MetaLength()
  data := []interface{}{}

  for j := 0; j < m; j++ {
    switch v := meta.MetaData[j].(type) {
    case [][]string:
      l := make([][]string, n)
      for i := 0; i < n; i++ {
        l[i] = v[indices[i]]
      }
      data = append(data, l)
    case []string :
      l := make([]string, n)
      for i := 0; i < n; i++ {
        l[i] = v[indices[i]]
      }
      data = append(data, l)
    case []float64:
      l := make([]float64, n)
      for i := 0; i < n; i++ {
        l[i] = v[indices[i]]
      }
      data = append(data, l)
    case [][]float64:
      l := make([][]float64, n)
      for i := 0; i < n; i++ {
        l[i] = v[indices[i]]
      }
      data = append(data, l)
    case []int    :
      l := make([]int, n)
      for i := 0; i < n; i++ {
        l[i] = v[indices[i]]
      }
      data = append(data, l)
    case [][]int    :
      l := make([][]int, n)
      for i := 0; i < n; i++ {
        l[i] = v[indices[i]]
      }
      data = append(data, l)
    }
  }
  return NewMeta(meta.MetaName, data)
}

// Return a new Meta object containing rows given by the range
// [ifrom, ito).
func (meta *Meta) Slice(ifrom, ito int) Meta {
  n := ito-ifrom
  m := meta.MetaLength()
  data := []interface{}{}

  for j := 0; j < m; j++ {
    switch v := meta.MetaData[j].(type) {
    case [][]string:
      l := make([][]string, n)
      for i := ifrom; i < ito; i++ {
        l[i-ifrom] = v[i]
      }
      data = append(data, l)
    case []string :
      l := make([]string, n)
      for i := ifrom; i < ito; i++ {
        l[i-ifrom] = v[i]
      }
      data = append(data, l)
    case []float64:
      l := make([]float64, n)
      for i := ifrom; i < ito; i++ {
        l[i-ifrom] = v[i]
      }
      data = append(data, l)
    case [][]float64:
      l := make([][]float64, n)
      for i := ifrom; i < ito; i++ {
        l[i-ifrom] = v[i]
      }
      data = append(data, l)
    case []int    :
      l := make([]int, n)
      for i := ifrom; i < ito; i++ {
        l[i-ifrom] = v[i]
      }
      data = append(data, l)
    case [][]int    :
      l := make([][]int, n)
      for i := ifrom; i < ito; i++ {
        l[i-ifrom] = v[i]
      }
      data = append(data, l)
    }
  }
  return NewMeta(meta.MetaName, data)
}

// Return a new Meta object where a given set of rows has been merged. The argument indices
// should assign the same target index to all rows that should be merged, i.e.
// indices := []int{0, 1, 1, 2, 3}
// would merge rows 1 and 2 but leave all rows 0, 3 and 4 as they are. Rows are merged by
// replacing one-dimensional slices by two-dimensional slices. The Reduce{String,Float64,Int}
// methods may be used afterwards to apply a function to the merged data. A Meta object that
// already contains two-dimensional slices cannot be merged.
func (meta *Meta) Merge(indices []int) Meta {
  sliceMax := func(s []int) int {
    max := 0
    for _, v := range s {
      if v > max {
        max = v
      }
    }
    return max
  }
  n := sliceMax(indices)+1
  m := meta.MetaLength()
  data := []interface{}{}

  for j := 0; j < m; j++ {
    switch v := meta.MetaData[j].(type) {
    case []string :
      l := make([][]string, n)
      for i := 0; i < len(v); i++ {
        l[indices[i]] = append(l[indices[i]], v[i])
      }
      data = append(data, l)
    case []float64:
      l := make([][]float64, n)
      for i := 0; i < len(v); i++ {
        l[indices[i]] = append(l[indices[i]], v[i])
      }
      data = append(data, l)
    case []int    :
      l := make([][]int, n)
      for i := 0; i < len(v); i++ {
        l[indices[i]] = append(l[indices[i]], v[i])
      }
      data = append(data, l)
    case [][]string : panic("cannot merge [][]string")
    case [][]int    : panic("cannot merge [][]int")
    case [][]float64: panic("cannot merge [][]float64")
    }
  }
  return NewMeta(meta.MetaName, data)
}

// Reduce a column with a two-dimensional string slice by applying
// the given function f. If nameNew != "" the old meta column
// is kept.
func (meta *Meta) ReduceString(name, nameNew string, f func([]string) string) {
  t := meta.GetMeta(name).([][]string)
  r := make([]string, len(t))
  for i := 0; i < len(t); i++ {
    r[i] = f(t[i])
  }
  if nameNew != "" && name != nameNew {
    meta.AddMeta(nameNew, r)
  } else {
    meta.DeleteMeta(name)
    meta.AddMeta(name, r)
  }
}

// Reduce a column with a two-dimensional float64 slice by applying
// the given function f. If nameNew != "" the old meta column
// is kept.
func (meta *Meta) ReduceFloat(name, nameNew string, f func([]float64) float64) {
  t := meta.GetMeta(name).([][]float64)
  r := make([]float64, len(t))
  for i := 0; i < len(t); i++ {
    r[i] = f(t[i])
  }
  if nameNew != "" && name != nameNew {
    meta.AddMeta(nameNew, r)
  } else {
    meta.DeleteMeta(name)
    meta.AddMeta(name, r)
  }
}

// Reduce a column with a two-dimensional int slice by applying
// the given function f. If nameNew != "" the old meta column
// is kept.
func (meta *Meta) ReduceInt(name, nameNew string, f func([]int) int) {
  t := meta.GetMeta(name).([][]int)
  r := make([]int, len(t))
  for i := 0; i < len(t); i++ {
    r[i] = f(t[i])
  }
  if nameNew != "" && name != nameNew {
    meta.AddMeta(nameNew, r)
  } else {
    meta.DeleteMeta(name)
    meta.AddMeta(name, r)
  }
}

/* sorting
 * -------------------------------------------------------------------------- */

type metaPair struct {
  Key   int
  Value interface{}
}
type metaPairList []metaPair

func (p metaPairList) Len() int {
  return len(p)
}
func (p metaPairList) Less(i, j int) bool {
  switch p[i].Value.(type) {
  case float64: return p[i].Value.(float64) < p[j].Value.(float64)
  case int    : return p[i].Value.(int)     < p[j].Value.(int)
  case string : return p[i].Value.(string)  < p[j].Value.(string)
  default:
    panic("Invalid type for sorting!")
  }
}

func (p metaPairList) Swap(i, j int) {
  p[i], p[j] = p[j], p[i]
}

func (meta *Meta) sortedIndices(name string, reverse bool) ([]int, error) {
  var l metaPairList
  if t := meta.GetMeta(name); t != nil {
    switch s := t.(type) {
    case []float64: for i, v := range s { l = append(l, metaPair{i, v}) }
    case []int    : for i, v := range s { l = append(l, metaPair{i, v}) }
    case []string : for i, v := range s { l = append(l, metaPair{i, v}) }
    default:
      panic("Invalid type for sorting!")
    }
  } else {
    return []int{}, errors.New("Meta column not found!")
  }
  if reverse {
    sort.Sort(sort.Reverse(l))
  } else {
    sort.Sort(l)
  }
  // slice of indices
  j := make([]int, len(l))
  // extract indices
  for i := 0; i < len(l); i++ {
    j[i] = l[i].Key
  }
  return j, nil
}

func (meta *Meta) Sort(name string, reverse bool) (Meta, error) {
  j, err := meta.sortedIndices(name, reverse)
  if err != nil {
    return Meta{}, err
  }
  return meta.Subset(j), nil
}

/* i/o
 * -------------------------------------------------------------------------- */

func (meta *Meta) PrettyPrint(n int) string {
  var buffer bytes.Buffer
  var bufferHeader bytes.Buffer
  writer       := bufio.NewWriter(&buffer)
  writerHeader := bufio.NewWriter(&bufferHeader)
  m := meta.MetaLength()
  maxLengths := make([]int, m)

  printRow := func(writer io.Writer, i int) {
    if i != 0 {
      fmt.Fprintf(writer, "\n")
    }
    for j := 0; j < m; j++ {
      length := 0
      switch v := meta.MetaData[j].(type) {
      case []string :
        length, _ = fmt.Fprintf(writer, " %12s", v[i])
      case []float64:
        length, _ = fmt.Fprintf(writer, " %12f", v[i])
      case []int    :
        length, _ = fmt.Fprintf(writer, " %12d", v[i])
      case [][]string:
        for k := 0; k < len(v[i]); k++ {
          l, _ := fmt.Fprintf(writer, " %12s", v[i][k])
          length += l
        }
      case [][]float64:
        for k := 0; k < len(v[i]); k++ {
          l, _ := fmt.Fprintf(writer, " %12f", v[i][k])
          length += l
        }
      case [][]int:
        for k := 0; j < len(v[i]); k++ {
          l, _ := fmt.Fprintf(writer, " %12d", v[i][k])
          length += l
        }
      }
      if maxLengths[j] > length {
        format := fmt.Sprintf("%%%ds", maxLengths[j]-length)
        fmt.Fprintf(writer, format, "")
      }
      if maxLengths[j] < length {
        maxLengths[j] = length
      }
    }
  }
  // select rows to print
  printData := func(writer io.Writer, ) {
    if meta.Length() <= n+1 {
      // print all entries
      for i := 0; i < meta.Length(); i++ {
        printRow(writer, i)
      }
    } else {
      // print first n/2 rows
      for i := 0; i < n/2; i++ {
        printRow(writer, i)
      }
      fmt.Fprintf(writer, "\n %12s", "...")
      // print last n/2 rows
      for i := meta.Length() - n/2; i < meta.Length(); i++ {
        printRow(writer, i)
      }
    }
  }
  // determine column widths
  printData(writer)
  writer.Flush()
  buffer.Reset()
  // actually print data
  printData(writer)
  writer.Flush()
  // print header
  for j := 0; j < m; j++ {
    format := fmt.Sprintf("%%%ds", maxLengths[j])
    fmt.Fprintf(writerHeader, format, meta.MetaName[j])
  }
  writerHeader.WriteString("\n")
  writerHeader.Flush()

  return bufferHeader.String() + buffer.String()
}

func (meta *Meta) String() string {
  return meta.PrettyPrint(10)
}

func (meta *Meta) WriteTableRow(w io.Writer, i int) {
  if i == -1 {
    // write header
    for k := 0; k < meta.MetaLength(); k++ {
      fmt.Fprintf(w, " %10s", meta.MetaName[k])
    }
  } else {
    for k := 0; k < meta.MetaLength(); k++ {
      switch v := meta.MetaData[k].(type) {
      case []string : fmt.Fprintf(w, " %10s", v[i])
      case []float64: fmt.Fprintf(w, " %10f", v[i])
      case []int    : fmt.Fprintf(w, " %10d", v[i])
      case [][]string:
        for j := 0; j < len(v[i]); j++ {
          if j == 0 {
            fmt.Fprintf(w, " %s", v[i][j])
          } else {
            fmt.Fprintf(w, ",%s", v[i][j])
          }
        }
      case [][]float64:
        for j := 0; j < len(v[i]); j++ {
          if j == 0 {
            fmt.Fprintf(w, " %f", v[i][j])
          } else {
            fmt.Fprintf(w, ",%f", v[i][j])
          }
        }
      case [][]int:
        for j := 0; j < len(v[i]); j++ {
          if j == 0 {
            fmt.Fprintf(w, " %d", v[i][j])
          } else {
            fmt.Fprintf(w, ",%d", v[i][j])
          }
        }
      }
    }
  }
}

func ReadMetaFromTable(filename string, names, types []string) (Meta, error) {
  if len(names) != len(types) {
    panic("invalid arguments")
  }
  result := Meta{}
  // header information
   idxMap := make(map[string]int)
  metaMap := make(map[string]interface{})
  for i := 0; i < len(names); i++ {
    idxMap[names[i]] = -1
    switch types[i] {
    case "[]string":  metaMap[names[i]] = []string{}
    case "[]int":     metaMap[names[i]] = []int{}
    case "[]float64": metaMap[names[i]] = []float64{}
    case "[][]string":  metaMap[names[i]] = [][]string{}
    case "[][]int":     metaMap[names[i]] = [][]int{}
    case "[][]float64": metaMap[names[i]] = [][]float64{}
    default: panic("invalid types argument")
    }
  }

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
    for i := 0; i < len(fields); i++ {
      if _, ok := idxMap[fields[i]]; ok {
        idxMap[fields[i]] = i
      }
    }
  }
  for scanner.Scan() {
    fields := strings.Fields(scanner.Text())
    if len(fields) == 0 {
      continue
    }
    for name, idx := range idxMap {
      if idx == -1 {
        // column not found, skip
        continue
      }
      if idx >= len(fields) {
        return result, errors.New("invalid table")
      }
      switch entry := metaMap[name].(type) {
      case []string:
        entry = append(entry, fields[idx])
        metaMap[name] = entry
      case []int:
        v, err := strconv.ParseInt(fields[idx], 10, 64)
        if err != nil {
          return result, err
        }
        entry = append(entry, int(v))
        metaMap[name] = entry
      case []float64:
        v, err := strconv.ParseFloat(fields[idx], 64)
        if err != nil {
          return result, err
        }
        entry = append(entry, v)
        metaMap[name] = entry
      case [][]int:
        data := strings.FieldsFunc(fields[idx], func(x rune) bool { return x == ',' })
        // parse counts
        entry = append(entry, make([]int, len(data)))
        // loop over count vector
        for i := 0; i < len(data); i++ {
          v, err := strconv.ParseInt(data[i], 10, 64)
          if err != nil {
            return result, err
          }
          entry[len(entry)-1][i] = int(v)
        }
        metaMap[name] = entry
      case [][]float64:
        data := strings.FieldsFunc(fields[idx], func(x rune) bool { return x == ',' })
        // parse counts
        entry = append(entry, make([]float64, len(data)))
        // loop over count vector
        for i := 0; i < len(data); i++ {
          v, err := strconv.ParseFloat(data[i], 64)
          if err != nil {
            return result, err
          }
          entry[len(entry)-1][i] = v
        }
        metaMap[name] = entry
      case [][]string:
        data := strings.FieldsFunc(fields[idx], func(x rune) bool { return x == ',' })
        // parse counts
        entry = append(entry, make([]string, len(data)))
        // loop over count vector
        for i := 0; i < len(data); i++ {
          entry[len(entry)-1][i] = data[i]
        }
        metaMap[name] = entry
      }
    }
  }
  for name, idx := range idxMap {
    if idx != -1 {
      result.AddMeta(name, metaMap[name])
    }
  }
  return result, nil
}
