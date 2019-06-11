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
import "sort"
import "strings"

/* -------------------------------------------------------------------------- */

type KmersCounter struct {
  n, m        int              // min and max kmer size
  length      int              // code length
  names     []string           // kmer names
  p         []int              // pre-evaluated powers
  kmap      []map[string][]int // map kmer strings to indices
  complement  bool
  reverse     bool
  revcomp     bool
  ma        []int              // maximum number of ambiguous letters
  na          int              // current number of ambiguous letters
  al          ComplementableAlphabet
}

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

func (obj KmersIterator) Get() []byte {
  return obj.c
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

func NewKmersCylinderIterator(k int, maxAmbiguous int, alphabet ComplementableAlphabet, j int, a []byte) KmersCylinderIterator {
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

func (obj KmersCylinderIterator) Get() []byte {
  return obj.c
}

func (obj KmersCylinderIterator) Ok() bool {
  return obj.ok
}

func (obj *KmersCylinderIterator) Next() {
  k := len(obj.c)
  // increment d
  for i := 0; i < k; {
    // skip fixed sub-kmer
    for i >= obj.j && i < obj.m {
      i++
    }
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

func NewKmersCounter(n, m int, comp, rev, rc bool, maxAmbiguous []int, al ComplementableAlphabet) (KmersCounter, error) {
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
    return KmersCounter{}, fmt.Errorf("parameter `maxAmbiguous' has invalid length")
  }
  r := KmersCounter{n: n, m: m, complement: comp, reverse: rev, revcomp: rc, ma: maxAmbiguous, al: al}
  p := make([]int, m+1)
  for k := 0; k <= m; k++ {
    p[k] = iPow(r.al.Length(), k)
  }
  r.kmap = make([]map[string][]int, m-n+1)
  names := []string{}
  idx   := 0
  for k := n; k <= m; k++ {
    r.kmap[k-n] = make(map[string][]int)
    c2 := make([]byte, k)
    c3 := make([]byte, k)
    c4 := make([]byte, k)
    cr := make([]byte, k)
    it := NewKmersIterator(k, maxAmbiguous[k-n], al)
    for c1 := it.Get(); it.Ok(); it.Next() {
      // do not allow gaps at the ends
      if ok, _ := r.al.IsWildcard(c1[  0]); ok {
        continue
      }
      if ok, _ := r.al.IsWildcard(c1[k-1]); ok {
        continue
      }
      // compute indices
      i     := 0 // index
      i_c   := 0 // index of complement
      i_r   := 0 // index of reverse
      i_rc  := 0 // index of reverse complement
      i_res := 0 // final index
      for j := 0; j < k; j++ {
        x1, _ := r.al.Code(c1[    j])
        x2, _ := r.al.Code(c1[k-j-1])
        y1, _ := r.al.ComplementCoded(x1)
        y2, _ := r.al.ComplementCoded(x2)
        i    += int(x2) * p[j]
        i_c  += int(y2) * p[j]
        i_r  += int(x1) * p[j]
        i_rc += int(y1) * p[j]
      }
      i_res = i
      // find minimum
      if comp && i_res > i_c {
        i_res = i_c
      }
      if rev && i_res > i_r {
        i_res = i_r
      }
      if rc && i_res > i_rc {
        i_res = i_rc
      }
      if i_res == i {
        // new kmer found
        r.addEquivalentKmers(cr, c1, idx)
        // compute strings
        if comp || rc {
          if err := r.comp(c2, c1); err != nil {
            return r, err
          }
        }
        name := string(c1)
        if comp {
          if string(c2) != string(c1) {
            r.addEquivalentKmers(cr, c2, idx)
          }
          name = fmt.Sprintf("%s|%s", name, string(c2))
        }
        if rev {
          r.rev(c3, c1)
          if comp {
            if string(c3) != string(c2) && string(c3) != string(c1) {
              r.addEquivalentKmers(cr, c3, idx)
            }
          } else {
            if string(c3) != string(c1) {
              r.addEquivalentKmers(cr, c3, idx)
            }
          }
          name = fmt.Sprintf("%s|%s", name, string(c3))
        }
        if rc {
          r.rev(c4, c2)
          if comp {
            if string(c4) != string(c3) && string(c4) != string(c2) && string(c4) != string(c1) {
              r.addEquivalentKmers(cr, c4, idx)
            }
          } else {
            if string(c4) != string(c3) && string(c4) != string(c1) {
              r.addEquivalentKmers(cr, c4, idx)
            }
          }
          name = fmt.Sprintf("%s|%s", name, string(c4))
        }
        names = append(names, name)
        idx  += 1
      }
    }
  }
  r.length  = idx
  r.p       = p
  r.names   = names
  return r, nil
}

/* -------------------------------------------------------------------------- */

func (obj KmersCounter) CountKmers(result []int, sequence []byte) error {
  if len(result) != obj.Length() {
    return fmt.Errorf("result slice has invalid length")
  }
  c := strings.ToLower(string(sequence))
  // loop over sequence
  for i := 0; i < len(c); i++ {
    // loop over all k-mers
    for k := obj.n; k <= obj.m && i+k-1 < len(c); k++ {
      for _, j := range obj.kmap[k-obj.n][c[i:i+k]] {
        result[j] += 1
      }
    }
  }
  return nil
}

func (obj KmersCounter) CountKmersSparse(sequence []byte) ([]int, []int) {
  c := strings.ToLower(string(sequence))
  r := make(map[int]int)
  // loop over sequence
  for i := 0; i < len(c); i++ {
    // loop over all k-mers
    for k := obj.n; k <= obj.m && i+k-1 < len(c); k++ {
      for _, j := range obj.kmap[k-obj.n][c[i:i+k]] {
        r[j] += 1
      }
    }
  }
  indices := []int{}
  counts  := []int{}
  for k, v := range r {
    indices = append(indices, k)
    counts  = append(counts,  v)
  }
  sort.Sort(sortIntPairs{indices, counts})
  return indices, counts
}

func (obj KmersCounter) IdentifyKmers(result []int, sequence []byte) error {
  if len(result) != obj.Length() {
    return fmt.Errorf("result slice has invalid length")
  }
  c := strings.ToLower(string(sequence))
  // loop over sequence
  for i := 0; i < len(c); i++ {
    // loop over all k-mers
    for k := obj.n; k <= obj.m && i+k-1 < len(c); k++ {
      for _, j := range obj.kmap[k-obj.n][c[i:i+k]] {
        result[j] = 1
      }
    }
  }
  return nil
}

func (obj KmersCounter) IdentifyKmersSparse(sequence []byte) []int {
  c := strings.ToLower(string(sequence))
  r := make(map[int]struct{})
  // loop over sequence
  for i := 0; i < len(c); i++ {
    // loop over all k-mers
    for k := obj.n; k <= obj.m && i+k-1 < len(c); k++ {
      for _, j := range obj.kmap[k-obj.n][c[i:i+k]] {
        r[j] = struct{}{}
      }
    }
  }
  indices := []int{}
  for k, _ := range r {
    indices = append(indices, k)
  }
  sort.Ints(indices)
  return indices
}

/* -------------------------------------------------------------------------- */

func (obj KmersCounter) Length() int {
  return obj.length
}

func (obj KmersCounter) KmerName(i int) string {
  return obj.names[i]
}

func (obj KmersCounter) MinKmerSize() int {
  return obj.n
}

func (obj KmersCounter) MaxKmerSize() int {
  return obj.m
}

func (obj KmersCounter) MaxAmbiguous() []int {
  return obj.ma
}

func (obj KmersCounter) Complement() bool {
  return obj.complement
}

func (obj KmersCounter) Reverse() bool {
  return obj.reverse
}

func (obj KmersCounter) Revcomp() bool {
  return obj.revcomp
}

func (obj KmersCounter) Alphabet() ComplementableAlphabet {
  return obj.al
}

/* -------------------------------------------------------------------------- */

func (obj KmersCounter) comp(dest, src []byte) error {
  for j := 0; j < len(src); j++ {
    if x, err := obj.al.Complement(src[j]); err != nil {
      return err
    } else {
      dest[j] = x
    }
  }
  return nil
}

func (obj KmersCounter) rev(dest, src []byte) {
  for i, j := 0, len(src)-1; i <= j; i, j = i+1, j-1 {
    dest[i], dest[j] = src[j], src[i]
  }
}

func (obj KmersCounter) addEquivalentKmersRec(dest, src []byte, i, idx int) {
  if i == len(src) {
    obj.kmap[len(dest)-obj.n][string(dest)] = append(obj.kmap[len(dest)-obj.n][string(dest)], idx)
  } else {
    x, _ := obj.al.Bases(src[i])
    for _, k := range x {
      dest[i] = k
      obj.addEquivalentKmersRec(dest, src, i+1, idx)
    }
  }
}

func (obj KmersCounter) addEquivalentKmers(dest, src []byte, idx int) {
  obj.addEquivalentKmersRec(dest, src, 0, idx)
}
