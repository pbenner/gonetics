/* Copyright (C) 2019 Philipp Benner
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

type Kmer string

/* -------------------------------------------------------------------------- */

func (obj Kmer) matches(kmer1, kmer2 []byte, alphabet ComplementableAlphabet) bool {
  if len(kmer1) != len(kmer2) {
    return false
  }
  for i := 0; i < len(kmer1); i++ {
    if kmer1[i] == kmer2[i] {
      continue
    }
    b1, err1 := alphabet.Bases(kmer1[i])
    b2, err2 := alphabet.Bases(kmer2[i])
    if err1 != nil {
      panic(err1)
    }
    if err2 != nil {
      panic(err2)
    }
    // check if b1 is a superset of b2
    if !byteSuperset(b1, b2) {
      return false
    }
  }
  return true
}

func (obj Kmer) Matches(b Kmer, alphabet ComplementableAlphabet) bool {
  kmer1 := []byte(obj)
  kmer2 := []byte(b)
  if len(kmer1) > len(kmer2) {
    panic("kmer1 must be smaller or equal in length than kmer2")
  }
  for i := 0; i < len(kmer2)-len(kmer1)+1; i++ {
    if obj.matches(kmer1, kmer2[0:i+len(kmer1)], alphabet) {
      return true
    }
  }
  return false
}
