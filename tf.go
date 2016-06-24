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

type TFMatrix [][]float64

/* -------------------------------------------------------------------------- */

func EmptyTFMatrix() TFMatrix {
  return TFMatrix{}
}

/* -------------------------------------------------------------------------- */

func (t TFMatrix) Get(c byte, i int) float64 {
  switch c {
  case 'A': fallthrough
  case 'a':
    return t[0][i]
  case 'C': fallthrough
  case 'c':
    return t[1][i]
  case 'G': fallthrough
  case 'g':
    return t[2][i]
  case 'T': fallthrough
  case 't':
    return t[3][i]
  default:
    panic("Get(): invalid parameters")
  }
}

func (t TFMatrix) GetRow(c byte) []float64 {
  switch c {
  case 'A': fallthrough
  case 'a':
    return t[0]
  case 'C': fallthrough
  case 'c':
    return t[1]
  case 'G': fallthrough
  case 'g':
    return t[2]
  case 'T': fallthrough
  case 't':
    return t[3]
  default:
    panic("GetRow(): invalid parameters")
  }
}

/* -------------------------------------------------------------------------- */

func (tp *TFMatrix) ReadMatrix(filename string) error {

  ncols := -1
  // allocate memory
  *tp = make([][]float64, 4)
  // dereference pointer to slice
  t := *tp

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
      t[0] = data
    case 'C': fallthrough
    case 'c':
      t[1] = data
    case 'G': fallthrough
    case 'g':
      t[2] = data
    case 'T': fallthrough
    case 't':
      t[3] = data
    default:
      return fmt.Errorf("ReadMatrix(): invalid tf matrix")
    }
  }
  return nil
}
