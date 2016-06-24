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
import "compress/gzip"
import "fmt"
import "os"
import "strconv"
import "strings"

/* -------------------------------------------------------------------------- */

type TFMatrix struct {
  Values [][]float64
}

/* -------------------------------------------------------------------------- */

func EmptyTFMatrix() TFMatrix {
  return TFMatrix{}
}

/* -------------------------------------------------------------------------- */

func (t TFMatrix) Length() int {
  if len(t.Values) == 0 {
    return -1
  }
  return len(t.Values[0])
}

func (t TFMatrix) Get(c byte, i int) float64 {
  switch c {
  case 'A': fallthrough
  case 'a':
    return t.Values[0][i]
  case 'C': fallthrough
  case 'c':
    return t.Values[1][i]
  case 'G': fallthrough
  case 'g':
    return t.Values[2][i]
  case 'T': fallthrough
  case 't':
    return t.Values[3][i]
  default:
    panic("Get(): invalid parameters")
  }
}

func (t TFMatrix) GetRow(c byte) []float64 {
  switch c {
  case 'A': fallthrough
  case 'a':
    return t.Values[0]
  case 'C': fallthrough
  case 'c':
    return t.Values[1]
  case 'G': fallthrough
  case 'g':
    return t.Values[2]
  case 'T': fallthrough
  case 't':
    return t.Values[3]
  default:
    panic("GetRow(): invalid parameters")
  }
}

/* -------------------------------------------------------------------------- */

func (t *TFMatrix) ReadMatrix(filename string) error {

  ncols := -1
  // allocate memory
  t.Values = make([][]float64, 4)

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

  for scanner.Scan() {
    fields := strings.Fields(scanner.Text())
    // if empty line, continue scanning
    if len(fields) == 0 {
      continue
    }
    if len(fields) <= 1 {
      return fmt.Errorf("ReadMatrix(): invalid tf matrix")
    }
    // if first line, set number of columns
    if ncols == -1 {
      ncols = len(fields)-1
    }
    if len(fields) != ncols+1 {
      return fmt.Errorf("ReadMatrix(): invalid tf matrix")
    }
    data := []float64{}
    // read one row of the matrix
    for i := 1; i < len(fields); i++ {
      v, err := strconv.ParseFloat(fields[i], 64)
      if err != nil {
        return err
      }
      data = append(data, v)
    }
    switch fields[0][0] {
    case 'A': fallthrough
    case 'a':
      t.Values[0] = data
    case 'C': fallthrough
    case 'c':
      t.Values[1] = data
    case 'G': fallthrough
    case 'g':
      t.Values[2] = data
    case 'T': fallthrough
    case 't':
      t.Values[3] = data
    default:
      return fmt.Errorf("ReadMatrix(): invalid tf matrix")
    }
  }
  return nil
}

/* scanning
 * -------------------------------------------------------------------------- */

type PWM struct {
  TFMatrix
}

func (t PWM) Scan(sequence []byte) []float64 {
  // number of positions where the pwm could fit
  n := len(sequence)-t.Length()+1
  // allocate memory
  result := make([]float64, n)
  // loop over sequence
  for i := 0; i < n; i++ {
    // loop over pwm
    for j := 0; j < t.Length(); j++ {
      result[i] += t.Get(sequence[i+j], j)
    }
  }
  return result
}
