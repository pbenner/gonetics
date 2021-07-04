/* Copyright (C) 2021 Philipp Benner
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

/* -------------------------------------------------------------------------- */

type grangesSort struct {
  GRanges
  indices []int
}

/* -------------------------------------------------------------------------- */

func newGRangesSort(g GRanges) grangesSort {
  indices := make([]int, g.Length())
  for i := 0; i < len(indices); i++ {
    indices[i] = i
  }
  return grangesSort{g, indices}
}

/* -------------------------------------------------------------------------- */

func (r grangesSort) Len() int {
  return r.Length()
}

func (r grangesSort) Less(i, j int) bool {
  si := r.Seqnames[r.indices[i]]
  sj := r.Seqnames[r.indices[j]]
  li := len(si)
  lj := len(sj)
  if li < lj {
    return true
  } else
  if li > lj {
    return false
  } else {
    if si < sj {
      return true
    } else
    if si > sj {
      return false
    } else {
      fi := r.Ranges[r.indices[i]].From
      fj := r.Ranges[r.indices[j]].From
      if fi < fj {
        return true
      } else
      if fi > fj {
        return false
      } else {
        ti := r.Ranges[r.indices[i]].To
        tj := r.Ranges[r.indices[j]].To
        return ti < tj
      }
    }
  }
}

func (r grangesSort) Swap(i, j int) {
  r.indices[i], r.indices[j] = r.indices[j], r.indices[i]
}
