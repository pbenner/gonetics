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

type KmerIterator struct {
  c  []byte
  al   ComplementableAlphabet
  ok   bool
  ma   int
  na   int
}

func NewKmerIterator(k int, maxAmbiguous int, alphabet ComplementableAlphabet) KmerIterator {
  c := make([]byte, k)
  for i := 0; i < k; i++ {
    c[i], _ = alphabet.Decode(0)
  }
  return KmerIterator{c: c, al: alphabet, ok: true, ma: maxAmbiguous}
}

func (obj KmerIterator) Get() string {
  return string(obj.c)
}

func (obj KmerIterator) Ok() bool {
  return obj.ok
}

func (obj *KmerIterator) Next() {
  k := len(obj.c)
  // increment d
  for i := 0; i < k; {
    ret  := false
    t, _ := obj.al.IsAmbiguous(obj.c[k-i-1])
    if obj.incrementPosition(k-i-1) {
      ret = true
    }
    // update ambiguous letter counter
    if s, _ := obj.al.IsAmbiguous(obj.c[k-i-1]); s {
      if !t { obj.na += 1 }
    } else {
      if  t { obj.na -= 1 }
    }
    if obj.ma >= 0 && obj.na > obj.ma {
      ret = false
    } else {
      i  += 1
    }
    if ret {
      return
    }
  }
  obj.ok = false
}

func (obj KmerIterator) incrementPosition(i int) bool {
  if c, _ := obj.al.Code(obj.c[i]); int(c+1) < obj.al.Length() {
    obj.c[i], _ = obj.al.Decode(c+1)
    return true
  } else {
    obj.c[i], _ = obj.al.Decode(0)
    return false
  }
}

/* -------------------------------------------------------------------------- */

type KmerCylinderIterator struct {
  c  []byte
  al   ComplementableAlphabet
  ok   bool
  ma   int
  na   int
  j    int
  m    int
}

func NewKmerCylinderIterator(k int, maxAmbiguous int, alphabet ComplementableAlphabet, j int, a_ string) KmerCylinderIterator {
  a := []byte(a_)
  m := j + len(a)
  c := make([]byte, k)
  if m > k {
    panic("invalid parameters")
  }
  for i := 0; i < k; i++ {
    c[i], _ = alphabet.Decode(0)
  }
  for i := j; i < m; i++ {
    c[i] = a[i-j]
  }
  return KmerCylinderIterator{c: c, al: alphabet, ok: true, ma: maxAmbiguous, j: j, m: m}
}

func (obj KmerCylinderIterator) Get() string {
  return string(obj.c)
}

func (obj KmerCylinderIterator) Ok() bool {
  return obj.ok
}

func (obj *KmerCylinderIterator) Next() {
  k := len(obj.c)
  i := k-1
  // skip fixed sub-k-mer
  if i >= obj.j && i < obj.m {
    i = obj.j-1
  }
  // increment d
  for i >= 0 {
    ret  := false
    t, _ := obj.al.IsAmbiguous(obj.c[i])
    if obj.incrementPosition(i) {
      ret = true
    }
    // update ambiguous letter counter
    if s, _ := obj.al.IsAmbiguous(obj.c[i]); s {
      if !t { obj.na += 1 }
    } else {
      if  t { obj.na -= 1 }
    }
    if obj.ma >= 0 && obj.na > obj.ma {
      ret = false
    } else {
      i  -= 1
      // skip fixed sub-k-mer
      if i >= obj.j && i < obj.m {
        i = obj.j-1
      }
    }
    if ret {
      return
    }
  }
  obj.ok = false
}

func (obj KmerCylinderIterator) incrementPosition(i int) bool {
  if c, _ := obj.al.Code(obj.c[i]); int(c+1) < obj.al.Length() {
    obj.c[i], _ = obj.al.Decode(c+1)
    return true
  } else {
    obj.c[i], _ = obj.al.Decode(0)
    return false
  }
}

/* -------------------------------------------------------------------------- */

type KmerInstantiationIterator struct {
  a     []byte
  c     []byte
  b   [][]byte
  j     []int
  ok      bool
  partial bool
  al      ComplementableAlphabet
}

func NewKmerInstantiationIterator(alphabet ComplementableAlphabet, a_ string, partial bool) KmerInstantiationIterator {
  a  := []byte(a_)
  k  := len(a)
  c  := make(  []byte, k)
  b  := make([][]byte, k)
  j  := make(  []int , k)
  ok := false
  copy(c, a)
  for i, _ := range a {
    // skip if character is not ambiguous
    if ok, err := alphabet.IsAmbiguous(a[i]); err != nil {
      panic(err)
    } else {
      if !ok {
        continue
      }
    }
    if r, err := alphabet.Bases(a[i]); err != nil {
      panic(err)
    } else {
      if partial {
        b[i] = r
      } else {
        b[i] = r
        c[i] = r[0]
        j[i] = 1
      }
    }
    ok = true
  }
  
  return KmerInstantiationIterator{a: a, c: c, b: b, j: j, ok: ok, partial: partial, al: alphabet}
}

func (obj KmerInstantiationIterator) Get() string {
  return string(obj.c)
}

func (obj KmerInstantiationIterator) Ok() bool {
  return obj.ok
}

func (obj *KmerInstantiationIterator) Next() {
  // find next position where an update is possible
  i := len(obj.b)-1
  for i >= 0 && (len(obj.b[i]) == 0 || len(obj.b[i]) == obj.j[i]) {
    i--
  }
  if i < 0 {
    obj.ok = false
    obj.c  = nil
  } else {
    obj.c[i] = obj.b[i][obj.j[i]]
    for i_ := len(obj.b)-1; i_ > i; i_-- {
      if len(obj.b[i_]) > 0 {
        if obj.partial {
          obj.j[i_] = 0
          obj.c[i_] = obj.a[i_]
        } else {
          obj.j[i_] = 1
          obj.c[i_] = obj.b[i_][0]
        }
      }
    }
    obj.j[i]++
  }
}
