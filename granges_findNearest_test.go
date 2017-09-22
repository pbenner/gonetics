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

//import "fmt"
import "testing"

/* -------------------------------------------------------------------------- */

func TestNearest(t *testing.T) {

  rQuery := NewGRanges(
    []string{"chr4", "chr4"},
    []int{600, 850},
    []int{950, 950},
    []byte{})
  rSubjects := NewGRanges(
    []string{"chr4", "chr4", "chr4", "chr4"},
    []int{100, 200, 300, 400},
    []int{900, 300, 700, 600},
    []byte{})

  queryHits, subjectHits, distances := FindNearest(rQuery, rSubjects, 2)

  rq := []int{0, 0, 0, 0, 1, 1, 1}
  rs := []int{0, 2, 3, 1, 0, 2, 3}
  rd := []int{0, 0, 0, 300, 0, 150, 250}
  rn := len(rq)

  if len(queryHits) != rn || len(subjectHits) != rn || len(distances) != rn {
    t.Error("test failed")
  } else {
    for i := 0; i < rn; i++ {
      if rq[i] != queryHits[i] {
        t.Error("test failed")
      }
      if rs[i] != subjectHits[i] {
        t.Error("test failed")
      }
      if rd[i] != distances[i] {
        t.Error("test failed")
      }
    }
  }
}
