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

package progress

/* -------------------------------------------------------------------------- */

import "bytes"
import "bufio"
import "fmt"
import "os"

/* -------------------------------------------------------------------------- */

type Progress struct {
  N, K, LineWidth int
}

/* -------------------------------------------------------------------------- */

func New(n, k int) Progress {
  progress := Progress{n, n/k, 40}
  if k > n {
    progress.K = 1
  }
  return progress
}

/* -------------------------------------------------------------------------- */

const __line_del__ = "\033[2K\r"

func (progress Progress) Exec(i int) string {
  var buffer bytes.Buffer
  writer := bufio.NewWriter(&buffer)

  p := float64(i)/float64(progress.N)
  // carriage return
  fmt.Fprintf(writer, "%s|", __line_del__)

  for i := 1; i < progress.LineWidth-1; i++ {
    if float64(i)/float64(progress.LineWidth) < p {
      fmt.Fprintf(writer, ">")
    } else {
      fmt.Fprintf(writer, " ")
    }
  }
  fmt.Fprintf(writer, "| %6.2f%%", p*100)
  // add newline if finished
  if p == 1.0 {
    fmt.Fprintf(writer, "\n")
  }
  writer.Flush()

  return buffer.String()
}

func (progress Progress) PrintStdout(i int) {
  if i == 0 || i == progress.N || (i % progress.K == 0) {
    fmt.Fprint(os.Stdout, progress.Exec(i))
  }
}

func (progress Progress) PrintStderr(i int) {
  if i == 0 || i == progress.N || (i % progress.K == 0) {
    fmt.Fprint(os.Stderr, progress.Exec(i))
  }
}
