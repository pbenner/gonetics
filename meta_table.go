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
import "io"
import "os"
import "strconv"
import "strings"

/* -------------------------------------------------------------------------- */

func (meta Meta) WriteTableRow(w io.Writer, i int) {
  if i == -1 {
    // write header
    for k := 0; k < meta.MetaLength(); k++ {
      fmt.Fprintf(w, " %s", meta.MetaName[k])
    }
  } else {
    for k := 0; k < meta.MetaLength(); k++ {
      switch v := meta.MetaData[k].(type) {
      case []string : fmt.Fprintf(w, " %s", v[i])
      case []float64: fmt.Fprintf(w, " %f", v[i])
      case []int    : fmt.Fprintf(w, " %d", v[i])
      case [][]string:
        if len(v[i]) == 0 {
          fmt.Fprintf(w, " nil")
        }
        for j := 0; j < len(v[i]); j++ {
          if j == 0 {
            fmt.Fprintf(w, " %s", v[i][j])
          } else {
            fmt.Fprintf(w, ",%s", v[i][j])
          }
        }
      case [][]float64:
        if len(v[i]) == 0 {
          fmt.Fprintf(w, " nil")
        }
        for j := 0; j < len(v[i]); j++ {
          if j == 0 {
            fmt.Fprintf(w, " %f", v[i][j])
          } else {
            fmt.Fprintf(w, ",%f", v[i][j])
          }
        }
      case [][]int:
        if len(v[i]) == 0 {
          fmt.Fprintf(w, " nil")
        }
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

func (meta *Meta) ReadTable(filename string, names, types []string) error {
  if len(names) != len(types) {
    panic("invalid arguments")
  }
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

  // scan header
  if scanner.Scan() {
    fields := strings.Fields(scanner.Text())
    // get number of columns
    if len(fields) < 4 {
      return fmt.Errorf("ReadTable(): invalid table `%s'", filename)
    }
    for i := 0; i < len(fields); i++ {
      if _, ok := idxMap[fields[i]]; ok {
        idxMap[fields[i]] = i
      }
    }
  }
  for i := 2; scanner.Scan(); i++ {
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
        return fmt.Errorf("ReadTable(): invalid table `%s'", filename)
      }
      switch entry := metaMap[name].(type) {
      case []string:
        entry = append(entry, fields[idx])
        metaMap[name] = entry
      case []int:
        v, err := strconv.ParseInt(fields[idx], 10, 64)
        if err != nil {
          return fmt.Errorf("parsing meta information failed at line `%d': %v", i, err)
        }
        entry = append(entry, int(v))
        metaMap[name] = entry
      case []float64:
        v, err := strconv.ParseFloat(fields[idx], 64)
        if err != nil {
          return fmt.Errorf("parsing meta information failed at line `%d': %v", i, err)
        }
        entry = append(entry, v)
        metaMap[name] = entry
      case [][]int:
        data := strings.FieldsFunc(fields[idx], func(x rune) bool { return x == ',' })
        // parse counts
        if len(data) == 1 && data[0] == "nil" {
          entry = append(entry, []int{})
        } else {
          entry = append(entry, make([]int, len(data)))
          // loop over count vector
          for i := 0; i < len(data); i++ {
            v, err := strconv.ParseInt(data[i], 10, 64)
            if err != nil {
              return fmt.Errorf("parsing meta information failed at line `%d': %v", i, err)
            }
            entry[len(entry)-1][i] = int(v)
          }
        }
        metaMap[name] = entry
      case [][]float64:
        data := strings.FieldsFunc(fields[idx], func(x rune) bool { return x == ',' })
        // parse counts
        if len(data) == 1 && data[0] == "nil" {
          entry = append(entry, []float64{})
        } else {
          entry = append(entry, make([]float64, len(data)))
          // loop over count vector
          for i := 0; i < len(data); i++ {
            v, err := strconv.ParseFloat(data[i], 64)
            if err != nil {
              return fmt.Errorf("parsing meta information failed at line `%d': %v", i, err)
            }
            entry[len(entry)-1][i] = v
          }
        }
        metaMap[name] = entry
      case [][]string:
        data := strings.FieldsFunc(fields[idx], func(x rune) bool { return x == ',' })
        // parse counts
        if len(data) == 1 && data[0] == "nil" {
          entry = append(entry, []string{})
        } else {
          entry = append(entry, make([]string, len(data)))
          // loop over count vector
          for i := 0; i < len(data); i++ {
            entry[len(entry)-1][i] = data[i]
          }
        }
        metaMap[name] = entry
      }
    }
  }
  for name, idx := range idxMap {
    if idx != -1 {
      meta.AddMeta(name, metaMap[name])
    }
  }
  return nil
}
