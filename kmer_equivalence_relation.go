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

type KmerEquivalence struct {
  N, M           int              // min and max kmer size
  Complement     bool
  Reverse        bool
  Revcomp        bool
  MaxAmbiguous []int              // maximum number of ambiguous letters
  Alphabet       ComplementableAlphabet
}

/* -------------------------------------------------------------------------- */

func NewKmerEquivalence(n, m int, comp, rev, rc bool, maxAmbiguous []int, al ComplementableAlphabet) (KmerEquivalence, error) {
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
    return KmerEquivalence{}, fmt.Errorf("parameter `maxAmbiguous' has invalid length")
  }
  r := KmerEquivalence{}
  r.N            = n
  r.M            = m
  r.Complement   = comp
  r.Reverse      = rev
  r.Revcomp      = rc
  r.MaxAmbiguous = maxAmbiguous
  r.Alphabet     = al
  return r, nil
}

/* -------------------------------------------------------------------------- */

func (a KmerEquivalence) Equals(b KmerEquivalence) bool {
  if a.M != b.M {
    return false
  }
  if a.N != b.N {
    return false
  }
  if a.Complement != b.Complement {
    return false
  }
  if a.Reverse != b.Reverse {
    return false
  }
  if a.Revcomp != b.Revcomp {
    return false
  }
  if a.Alphabet.String() != b.Alphabet.String() {
    return false
  }
  if len(a.MaxAmbiguous) != len(b.MaxAmbiguous) {
    return false
  }
  for i := 0; i < len(a.MaxAmbiguous); i++ {
    if a.MaxAmbiguous[i] != b.MaxAmbiguous[i] {
      return false
    }
  }
  return true
}

/* -------------------------------------------------------------------------- */

type KmerEquivalenceRelation struct {
  KmerEquivalence
  p []int // pre-evaluated powers
}

/* -------------------------------------------------------------------------- */

func NewKmerEquivalenceRelation(n, m int, comp, rev, rc bool, maxAmbiguous []int, al ComplementableAlphabet) (KmerEquivalenceRelation, error) {
  r := KmerEquivalenceRelation{}
  p := make([]int, m+1)
  for k := 0; k <= m; k++ {
    p[k] = iPow(al.Length(), k)
  }
  t, err := NewKmerEquivalence(n, m, comp, rev, rc, maxAmbiguous, al)
  if err != nil {
    return r, err
  }
  r.KmerEquivalence = t
  r.p               = p
  return r, nil
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
    x1, _ := obj.Alphabet.Code(c1[    j])
    x2, _ := obj.Alphabet.Code(c1[k-j-1])
    y1, _ := obj.Alphabet.ComplementCoded(x1)
    y2, _ := obj.Alphabet.ComplementCoded(x2)
    i    += int(x2) * obj.p[j]
    i_c  += int(y2) * obj.p[j]
    i_r  += int(x1) * obj.p[j]
    i_rc += int(y1) * obj.p[j]
  }
  // find minimum
  if obj.Complement && i > i_c {
    i = i_c
  } else {
    i_c = -1
  }
  if obj.Reverse && i > i_r {
    i = i_r
  } else {
    i_r = -1
  }
  if obj.Revcomp && i > i_rc {
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
  if !b_c2 && obj.Complement || obj.Revcomp {
    obj.comp(c2, c1)
  }
  if !b_c3 && obj.Reverse    || obj.Revcomp {
    obj.rev (c3, c1)
  }
  if !b_c4 && obj.Revcomp {
    obj.rev (c4, c2)
  }
  elements := []string{string(c1)}
  if obj.Complement {
    elements = append(elements, string(c2))
  }
  if obj.Reverse {
    elements = append(elements, string(c3))
  }
  if obj.Revcomp {
    elements = append(elements, string(c4))
  }
  return NewKmerClass(k, i, elements)
}

/* -------------------------------------------------------------------------- */

func (obj KmerEquivalenceRelation) comp(dest, src []byte) error {
  for j := 0; j < len(src); j++ {
    if x, err := obj.Alphabet.Complement(src[j]); err != nil {
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
