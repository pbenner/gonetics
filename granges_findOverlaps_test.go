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

func TestOverlapsList(t *testing.T) {

  s := NewEndPointList()

  r1 := endPoint{100, nil, nil, 1, true}
  r2 := endPoint{200, nil, nil, 1, true}
  r3 := endPoint{300, nil, nil, 1, true}
  r4 := endPoint{300, nil, nil, 1, true}

  s.Append(r1)
  s.Append(r2)
  s.Append(r3)

  s.Remove(r2)

  if s[0] != r1 {
    t.Error("TestOverlapsStack failed!")
  }
  if s[1] != r3 {
    t.Error("TestOverlapsStack failed!")
  }
  if s[1] != r4 {
    t.Error("TestOverlapsStack failed!")
  }
}

func TestOverlaps1(t *testing.T) {

  rSubjects := NewGRanges(
    []string{"chr4", "chr4", "chr4", "chr4"},
    []int{100, 200, 300, 400},
    []int{150, 250, 350, 450},
    []byte{})
  rQuery := NewGRanges(
    []string{"chr1", "chr4", "chr4", "chr4", "chr4", "chr4"},
    []int{100, 110, 190, 340, 390, 450},
    []int{150, 120, 220, 360, 400, 500},
    []byte{})

  queryHits, subjectHits := FindOverlaps(rQuery, rSubjects)

  if len(queryHits) != 3 {
    t.Error("TestOverlaps1 failed!")
  }
  if   queryHits[0] != 1 ||   queryHits[1] != 2 ||   queryHits[2] != 3 ||
    (subjectHits[0] != 0 || subjectHits[1] != 1 || subjectHits[2] != 2) {
    t.Error("TestOverlaps1 failed!")
  }
}

func TestOverlaps2(t *testing.T) {

  rSubjects := NewGRanges(
    []string{"chr4", "chr4", "chr4", "chr4"},
    []int{100, 200, 300, 400},
    []int{900, 800, 700, 600},
    []byte{})
  rQuery := NewGRanges(
    []string{"chr4", "chr4"},
    []int{600, 850},
    []int{950, 950},
    []byte{})

  queryHits, subjectHits := FindOverlaps(rQuery, rSubjects)

  if len(queryHits) != 4 {
    t.Error("TestOverlaps1 failed!")
  }
  if   queryHits[0] != 0 ||   queryHits[1] != 0 ||   queryHits[2] != 0 ||   queryHits[3] != 1 ||
    (subjectHits[0] != 0 || subjectHits[1] != 1 || subjectHits[2] != 2 || subjectHits[3] != 0) {
    t.Error("TestOverlaps2 failed!")
  }
}
