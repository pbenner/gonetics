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
import "math"
import "io"
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

func (t TFMatrix) Get(c byte, j int) float64 {
  i, err := NucleotideAlphabet{}.Code(c)
  if err != nil {
    panic(err)
  }
  return t.Values[i][j]
}

func (t TFMatrix) GetRow(c byte) []float64 {
  i, err := NucleotideAlphabet{}.Code(c)
  if err != nil {
    panic(err)
  }
  return t.Values[i]
}

func (t TFMatrix) RevComp() TFMatrix {
  alphabet := NucleotideAlphabet{}
  s := make([][]float64, alphabet.Length())
  for i := 0; i < alphabet.Length(); i++ {
    j, _ := alphabet.ComplementCoded(byte(i))
    s[j] = reverseFloat64(t.Values[i])
  }
  return TFMatrix{s}
}

func (t TFMatrix) Scan(sequence []byte, x0 float64, revcomp bool, f func(float64, float64) float64) float64 {
  if len(sequence) != t.Length() {
    panic("Scan(): sequence has invalid length")
  }
  alphabet := NucleotideAlphabet{}
  // number of positions where the pwm could fit
  x := x0
  if revcomp {
    // loop over pwm
    for j := 0; j < t.Length(); j++ {
      if a := sequence[t.Length()-j-1]; a != 'N' && a != 'n' {
        c, _ := alphabet.Complement(a)
        x = f(x, t.Get(c, j))
      }
    }
  } else {
    // loop over pwm
    for j := 0; j < t.Length(); j++ {
      if a := sequence[j]; a != 'N' && a != 'n' {
        x = f(x, t.Get(a, j))
      }
    }
  }
  return x
}

/* -------------------------------------------------------------------------- */

func (t *TFMatrix) ReadMatrix(reader io.Reader) error {

  scanner := bufio.NewScanner(reader)

  ncols := -1
  // allocate memory
  t.Values = make([][]float64, NucleotideAlphabet{}.Length())

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
    i, err := NucleotideAlphabet{}.Code(fields[0][0])
    if err != nil {
      return err
    }
    t.Values[i] = data
  }
  return nil
}

func (t *TFMatrix) ImportMatrix(filename string) error {
  var reader io.Reader
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
    reader = g
  } else {
    reader = f
  }
  return t.ReadMatrix(reader)
}

/* scanning
 * -------------------------------------------------------------------------- */

type PWM struct {
  TFMatrix
}

func (t PWM) MaxScore(sequence []byte, revcomp bool) float64 {
  // number of positions where the pwm could fit
  n := len(sequence)-t.Length()+1
  // function for adding scanning results
  f := func(a, b float64) float64 { return a+b }
  // maximum score
  result := math.Inf(-1)
  // loop over sequence
  for i := 0; i < n; i++ {
    if tmp := t.Scan(sequence[i:i+t.Length()], 0.0, revcomp, f); tmp > result {
      result = tmp
    }
  }
  return result
}
