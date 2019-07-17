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

type KmersIterator struct {
  c  []byte
  al   ComplementableAlphabet
  ok   bool
  ma   int
  na   int
}

func NewKmersIterator(k int, maxAmbiguous int, alphabet ComplementableAlphabet) KmersIterator {
  c := make([]byte, k)
  for i := 0; i < k; i++ {
    c[i], _ = alphabet.Decode(0)
  }
  return KmersIterator{c: c, al: alphabet, ok: true, ma: maxAmbiguous}
}

func (obj KmersIterator) Get() string {
  return string(obj.c)
}

func (obj KmersIterator) Ok() bool {
  return obj.ok
}

func (obj *KmersIterator) Next() {
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

func (obj KmersIterator) incrementPosition(i int) bool {
  if c, _ := obj.al.Code(obj.c[i]); int(c+1) < obj.al.Length() {
    obj.c[i], _ = obj.al.Decode(c+1)
    return true
  } else {
    obj.c[i], _ = obj.al.Decode(0)
    return false
  }
}

/* -------------------------------------------------------------------------- */

type KmersCylinderIterator struct {
  c  []byte
  al   ComplementableAlphabet
  ok   bool
  ma   int
  na   int
  j    int
  m    int
}

func NewKmersCylinderIterator(k int, maxAmbiguous int, alphabet ComplementableAlphabet, j int, a_ string) KmersCylinderIterator {
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
  return KmersCylinderIterator{c: c, al: alphabet, ok: true, ma: maxAmbiguous, j: j, m: m}
}

func (obj KmersCylinderIterator) Get() string {
  return string(obj.c)
}

func (obj KmersCylinderIterator) Ok() bool {
  return obj.ok
}

func (obj *KmersCylinderIterator) Next() {
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

func (obj KmersCylinderIterator) incrementPosition(i int) bool {
  if c, _ := obj.al.Code(obj.c[i]); int(c+1) < obj.al.Length() {
    obj.c[i], _ = obj.al.Decode(c+1)
    return true
  } else {
    obj.c[i], _ = obj.al.Decode(0)
    return false
  }
}

/* -------------------------------------------------------------------------- */

type KmersInstantiationIterator struct {
  c   []byte
  b [][]byte
  j     int
  ok    bool
  al    ComplementableAlphabet
}

func NewKmersInstantiationIterator(alphabet ComplementableAlphabet, a_ string) KmersInstantiationIterator {
  a  := []byte(a_)
  k  := len(a)
  c  := make(  []byte, k)
  b  := make([][]byte, k)
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
      c[i] = r[0]
      b[i] = r[1:]
    }
    ok = true
  }
  
  return KmersInstantiationIterator{c: c, b: b, j: 0, ok: ok, al: alphabet}
}

func (obj KmersInstantiationIterator) Get() string {
  return string(obj.c)
}

func (obj KmersInstantiationIterator) Ok() bool {
  return obj.ok
}

func (obj *KmersInstantiationIterator) Next() {
  for obj.j < len(obj.b) && len(obj.b[obj.j]) == 0 {
    obj.j++
  }
  if obj.j == len(obj.b) {
    obj.ok = false
    obj.c  = nil
  } else {
    obj.c[obj.j] = obj.b[obj.j][0]
    obj.b[obj.j] = obj.b[obj.j][1:]
  }
}
