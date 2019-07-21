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

import "fmt"

/* -------------------------------------------------------------------------- */

type KmerEquivalenceRelation struct {
  n, m        int              // min and max kmer size
  complement  bool
  reverse     bool
  revcomp     bool
  ma        []int              // maximum number of ambiguous letters
  al          ComplementableAlphabet
  p         []int              // pre-evaluated powers
}

/* -------------------------------------------------------------------------- */

func NewKmerEquivalenceRelation(n, m int, comp, rev, rc bool, maxAmbiguous []int, al ComplementableAlphabet) (KmerEquivalenceRelation, error) {
  if len(maxAmbiguous) == 0 {
    maxAmbiguous = make([]int, m-n+1)
    for i := 0; i < m-n+1; i++ {
      maxAmbiguous[i] = -1
    }
  } else
  if len(maxAmbiguous) == 1 {
    x := maxAmbiguous[0]
    maxAmbiguous = make([]int, m-n+1)
    for i := 0; i < m-n+1; i++ {
      maxAmbiguous[i] = x
    }
  } else
  if len(maxAmbiguous) != m-n+1 {
    return KmerEquivalenceRelation{}, fmt.Errorf("parameter `maxAmbiguous' has invalid length")
  }
  p := make([]int, m+1)
  for k := 0; k <= m; k++ {
    p[k] = iPow(al.Length(), k)
  }
  r := KmerEquivalenceRelation{
    n         : n,
    m         : m,
    p         : p,
    complement: comp,
    reverse   : rev,
    revcomp   : rc,
    ma        : maxAmbiguous,
    al        : al }
  return r, nil
}

/* -------------------------------------------------------------------------- */

func (obj KmerEquivalenceRelation) MinKmerSize() int {
  return obj.n
}

func (obj KmerEquivalenceRelation) MaxKmerSize() int {
  return obj.m
}

func (obj KmerEquivalenceRelation) MaxAmbiguous() []int {
  return obj.ma
}

func (obj KmerEquivalenceRelation) Complement() bool {
  return obj.complement
}

func (obj KmerEquivalenceRelation) Reverse() bool {
  return obj.reverse
}

func (obj KmerEquivalenceRelation) Revcomp() bool {
  return obj.revcomp
}

func (obj KmerEquivalenceRelation) Alphabet() ComplementableAlphabet {
  return obj.al
}

/* -------------------------------------------------------------------------- */

func (obj KmerEquivalenceRelation) EquivalenceClass(kmer string) KmerClass {
  k  := len(kmer)
  c1 := []byte(kmer)
  c2 := make([]byte, k)
  c3 := make([]byte, k)
  c4 := make([]byte, k)
  // compute indices
  i     := 0 // id
  i_c   := 0 // id of complement
  i_r   := 0 // id of reverse
  i_rc  := 0 // id of reverse complement
  for j := 0; j < k; j++ {
    x1, _ := obj.al.Code(c1[    j])
    x2, _ := obj.al.Code(c1[k-j-1])
    y1, _ := obj.al.ComplementCoded(x1)
    y2, _ := obj.al.ComplementCoded(x2)
    i    += int(x2) * obj.p[j]
    i_c  += int(y2) * obj.p[j]
    i_r  += int(x1) * obj.p[j]
    i_rc += int(y1) * obj.p[j]
  }
  // find minimum
  if obj.Complement() && i > i_c {
    i = i_c
  } else {
    i_c = -1
  }
  if obj.Reverse() && i > i_r {
    i = i_r
  } else {
    i_r = -1
  }
  if obj.Revcomp() && i > i_rc {
    i = i_rc
  } else {
    i_rc = -1
  }
  // compute name
  b_c2 := false
  b_c3 := false
  b_c4 := false
  switch i {
  case i_c:
    obj.comp(c2, c1)
    c1, c2 = c2, c1
    b_c2 = true
  case i_r:
    obj.rev (c3, c1)
    c1, c3 = c3, c1
    b_c3 = true
  case i_rc:
    obj.rev (c2, c1)
    obj.comp(c4, c2)
    c1, c4 = c4, c1
    b_c2 = true
    b_c4 = true
  }
  if !b_c2 && obj.Complement() || obj.Revcomp() {
    obj.comp(c2, c1)
  }
  if !b_c3 && obj.Reverse() || obj.Revcomp() {
    obj.rev (c3, c1)
  }
  if !b_c4 && obj.Revcomp() {
    obj.rev (c4, c2)
  }
  name := string(c1)
  if obj.Complement() {
    name = fmt.Sprintf("%s|%s", name, string(c2))
  }
  if obj.Reverse() {
    name = fmt.Sprintf("%s|%s", name, string(c3))
  }
  if obj.Revcomp() {
    name = fmt.Sprintf("%s|%s", name, string(c4))
  }
  return KmerClass{K: k, I: i, Name: name}
}

/* -------------------------------------------------------------------------- */

func (obj KmerEquivalenceRelation) comp(dest, src []byte) error {
  for j := 0; j < len(src); j++ {
    if x, err := obj.al.Complement(src[j]); err != nil {
      return err
    } else {
      dest[j] = x
    }
  }
  return nil
}

func (obj KmerEquivalenceRelation) rev(dest, src []byte) {
  for i, j := 0, len(src)-1; i <= j; i, j = i+1, j-1 {
    dest[i], dest[j] = src[j], src[i]
  }
}
