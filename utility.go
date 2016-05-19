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

/* -------------------------------------------------------------------------- */

func RemoveDuplicatesInt(s []int) []int {
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

// Divide a by b, the result is rounded down.
func divIntDown(a, b int) int {
  return a/b
}

// Divide a by b, the result is rounded up.
func divIntUp(a, b int) int {
  return (a-1)/b+1
}

/* -------------------------------------------------------------------------- */

func writeFile(filename string, r io.Reader, compress bool) {
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
  ioutil.WriteFile(filename, buffer.Bytes(), 0666)
}
