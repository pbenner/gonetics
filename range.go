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

/* -------------------------------------------------------------------------- */

type Range struct {
  From, To int
}

/* constructors
 * -------------------------------------------------------------------------- */

// Range object used to identify a genomic subsequence. By convention the first
// position in a sequence is numbered 0. The arguments from, to are interpreted
// as the interval [from, to).
func NewRange(from, to int) Range {
  if from > to {
    panic("NewRange(): from > to")
  }
  return Range{from, to}
}

/* -------------------------------------------------------------------------- */

func (r Range) Intersection(s Range) Range {
  from := iMax(r.From, s.From)
  to   := iMin(r.To,   s.To)
  // this shouldn't happen if r and s overlap
  if to < from {
    to = from
  }
  return NewRange(from, to)
}

/* -------------------------------------------------------------------------- */

func (r Range) String() string {
  return fmt.Sprintf("[%d %d)", r.From, r.To)
}
