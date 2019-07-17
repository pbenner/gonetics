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
import "compress/gzip"
import "io"
import "io/ioutil"
import "math"
import "regexp"
import "strings"
import "sort"
import "os"
import "unicode"

/* -------------------------------------------------------------------------- */

func removeDuplicatesInt(s []int) []int {
  m := map[int]bool{}
  r := []int{}

  for _, v := range s {
    if m[v] != true {
	    m[v] = true
	    r    = append(r, v)
    }
  }
  return r
}

/* -------------------------------------------------------------------------- */

func iMin(a, b int) int {
  if a < b {
    return a
  } else {
    return b
  }
}

func iMax(a, b int) int {
  if a > b {
    return a
  } else {
    return b
  }
}

func iPow(x, k int) int {
  return int(math.Pow(float64(x), float64(k)))
}

// Divide a by b, the result is rounded down.
func divIntDown(a, b int) int {
  return a/b
}

// Divide a by b, the result is rounded up.
func divIntUp(a, b int) int {
  return (a+b-1)/b
}

/* -------------------------------------------------------------------------- */

func writeFile(filename string, r io.Reader, compress bool) error {
  var buffer bytes.Buffer

  if compress {
    w := gzip.NewWriter(&buffer)
    io.Copy(w, r)
    w.Close()
  } else {
    w := bufio.NewWriter(&buffer)
    io.Copy(w, r)
    w.Flush()
  }
  return ioutil.WriteFile(filename, buffer.Bytes(), 0666)
}

func isGzip(filename string) bool {

  f, err := os.Open(filename)
  if err != nil {
    return false
  }
  defer f.Close()

  b := make([]byte, 2)
  n, err := f.Read(b)
  if err != nil {
    return false
  }

  if n == 2 && b[0] == 31 && b[1] == 139 {
    return true
  }
  return false
}

/* -------------------------------------------------------------------------- */

func fieldsQuoted(line string) []string {
  // if quoted
  q := false
  f := func(r rune) bool {
    if r == '"' {
      q = !q
    }
    return unicode.IsSpace(r) && q == false
  }
  return strings.FieldsFunc(line, f)
}

func removeQuotes(str string) string {
  reg := regexp.MustCompile(`"([^"]*)"`)
  if reg.MatchString(str) {
    return reg.ReplaceAllString(str, "${1}")
  }
  return str
}

func reverseFloat64(x []float64) []float64 {
  y := make([]float64, len(x))
  for i := 0; i < len(x); i++ {
    y[len(x)-i-1] = x[i]
  }
  return y
}

/* -------------------------------------------------------------------------- */

func bufioReadLine(reader *bufio.Reader) (string, error) {
  l, err := reader.ReadString('\n')
  if err != nil {
    // ignore EOF errors if some bytes were read
    if len(l) > 0 && err == io.EOF {
      return l, nil
    }
    return l, err
  }
  // remove newline character
  return l[0:len(l)-1], err
}

/* -------------------------------------------------------------------------- */

type sortIntPairs struct {
  a []int
  b []int
}

func (obj sortIntPairs) Len() int {
  return len(obj.a)
}

func (obj sortIntPairs) Less(i, j int) bool {
  return obj.a[i] < obj.a[j]
}

func (obj sortIntPairs) Swap(i, j int) {
  obj.a[i], obj.a[j] = obj.a[j], obj.a[i]
  obj.b[i], obj.b[j] = obj.b[j], obj.b[i]
}

func (obj sortIntPairs) Sort() {
  sort.Sort(obj)
}

func (obj sortIntPairs) SortRev() {
  sort.Sort(sort.Reverse(obj))
}
