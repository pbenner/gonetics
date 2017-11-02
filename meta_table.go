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
import "math"
import "io"
import "io/ioutil"
import "strconv"
import "strings"

/* -------------------------------------------------------------------------- */

type OptionGRangesScientific struct {
  Value bool
}

/* -------------------------------------------------------------------------- */

func (meta Meta) WriteTable(writer io.Writer, header bool, args ...interface{}) error {
  useScientific := false
  for _, arg := range args {
    switch a := arg.(type) {
    case OptionGRangesScientific:
      useScientific = a.Value
    default:
    }
  }
  printCellSlice := func(writer io.Writer, widths []int, i, j int, data interface{}) (int, error) {
    var tmpBuffer bytes.Buffer
    tmpWriter := bufio.NewWriter(&tmpBuffer)
    switch v := data.(type) {
    case [][]string:
      if len(v[i]) == 0 {
        if _, err := fmt.Fprintf(tmpWriter, "nil"); err != nil {
          return 0, err
        }
      }
      for k := 0; k < len(v[i]); k++ {
        if k != 0 {
          fmt.Fprintf(tmpWriter, ",")
        }
        if _, err := fmt.Fprintf(tmpWriter, "%s", v[i][k]); err != nil {
          return 0, err
        }
      }
    case [][]float64:
      if len(v[i]) == 0 {
        if _, err := fmt.Fprintf(tmpWriter, "nil"); err != nil {
          return 0, err
        }
      }
      for k := 0; k < len(v[i]); k++ {
        if k != 0 {
          if _, err := fmt.Fprintf(tmpWriter, ","); err != nil {
            return 0, err
          }
        }
        if useScientific {
          if _, err := fmt.Fprintf(tmpWriter, "%e", v[i][k]); err != nil {
            return 0, err
          }
        } else {
          if _, err := fmt.Fprintf(tmpWriter, "%f", v[i][k]); err != nil {
            return 0, err
          }
        }
      }
    case [][]int:
      if len(v[i]) == 0 {
        if _, err := fmt.Fprintf(tmpWriter, "nil"); err != nil {
          return 0, err
        }
      }
      for k := 0; k < len(v[i]); k++ {
        if k != 0 {
          fmt.Fprintf(tmpWriter, ",")
        }
        if _, err := fmt.Fprintf(tmpWriter, "%d", v[i][k]); err != nil {
          return 0, err
        }
      }
    default:
      panic("invalid meta data")
    }
    tmpWriter.Flush()
    format := fmt.Sprintf(" %%%ds", widths[j]-1)
    l, _ := fmt.Fprintf(writer, format, tmpBuffer.String())
    return l, nil
  }
  printCell := func(writer io.Writer, widths []int, i, j int) (int, error) {
    switch v := meta.MetaData[j].(type) {
    case []string:
      format := fmt.Sprintf(" %%%ds", widths[j]-1)
      return fmt.Fprintf(writer, format, v[i])
    case []float64:
      if useScientific {
        format := fmt.Sprintf(" %%%de", widths[j]-1)
        return fmt.Fprintf(writer, format, v[i])
      } else {
        format := fmt.Sprintf(" %%%df", widths[j]-1)
        return fmt.Fprintf(writer, format, v[i])
      }
    case []int:
      format := fmt.Sprintf(" %%%dd", widths[j]-1)
      return fmt.Fprintf(writer, format, v[i])
    default:
      return printCellSlice(writer, widths, i, j, v)
    }
  }
  printRow := func(writer io.Writer, widths []int, i int) error {
    if i != 0 {
      if _, err := fmt.Fprintf(writer, "\n"); err != nil {
        return err
      }
    }
    for j := 0; j < meta.MetaLength(); j++ {
      if _, err := printCell(writer, widths, i, j); err != nil {
        return err
      }
    }
    return nil
  }
  // compute widths of all cells in row i
  updateMaxWidths := func(i int, widths []int) error {
    for j := 0; j < meta.MetaLength(); j++ {
      if width, err := printCell(ioutil.Discard, widths, i, j); err != nil {
        return err
      } else {
        if width > widths[j] {
          widths[j] = width
        }
      }
    }
    return nil
  }
  printHeader := func(writer io.Writer, widths []int) error {
    for j := 0; j < meta.MetaLength(); j++ {
      format := fmt.Sprintf(" %%%ds", widths[j]-1)
      if _, err := fmt.Fprintf(writer, format, meta.MetaName[j]); err != nil {
        return err
      }
    }
    if _, err := fmt.Fprintf(writer, "\n"); err != nil {
      return err
    }
    return nil
  }
  applyRows := func(f1 func(i int) error) error {
    for i := 0; i < meta.Length(); i++ {
      if err := f1(i); err != nil {
        return err
      }
    }
    return nil
  }
  // maximum column widths
  widths := make([]int, meta.MetaLength())
  for j := 0; j < meta.MetaLength(); j++ {
    if width, err := fmt.Fprintf(ioutil.Discard, " %s", meta.MetaName[j]); err != nil {
      return err
    } else {
      widths[j] = width
    }
  }
  // determine column widths
  if err := applyRows(func(i int) error { return updateMaxWidths(i, widths) }); err != nil {
    return err
  }
  // pring header
  if header {
    if err := printHeader(writer, widths); err != nil {
      return err
    }
  }
  // print rows
  if err := applyRows(
    func(i int) error { return printRow(writer, widths, i) }); err != nil {
    return err
  }
  return nil
}

/* -------------------------------------------------------------------------- */

func (meta *Meta) PrintTable(header bool, args ...interface{}) string {
  var buffer bytes.Buffer
  writer := bufio.NewWriter(&buffer)

  if err := meta.WriteTable(writer, header, args...); err != nil {
    return ""
  }
  writer.Flush()

  return buffer.String()
}

/* -------------------------------------------------------------------------- */

func (meta *Meta) ReadTable(r io.Reader, names, types []string) error {
  scanner := bufio.NewScanner(r)

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
  // scan header
  if scanner.Scan() {
    fields := strings.Fields(scanner.Text())
    // get number of columns
    if len(fields) < 4 {
      return fmt.Errorf("invalid table")
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
        return fmt.Errorf("invalid table")
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
        if fields[idx] == "NA" || fields[idx] == "NaN" {
          entry = append(entry, math.NaN())
          metaMap[name] = entry
        } else {
          v, err := strconv.ParseFloat(fields[idx], 64)
          if err != nil {
            return fmt.Errorf("parsing meta information failed at line `%d': %v", i, err)
          }
          entry = append(entry, v)
          metaMap[name] = entry
        }
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
