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
  al          ComplementableAlphabet
}

/* -------------------------------------------------------------------------- */

func NewKmersCounter(n, m int, comp, rev, rc bool, al ComplementableAlphabet) (KmersCounter, error) {
  r := KmersCounter{n: n, m: m, complement: comp, reverse: rev, revcomp: rc, al: al}
  p := make([]int, m+1)
  for k := 0; k <= m; k++ {
    p[k] = iPow(r.al.Length(), k)
  }
  r.kmap = make([]map[string][]int, m-n+1)
  names := []string{}
  idx   := 0
  for k := n; k <= m; k++ {
    r.kmap[k-n] = make(map[string][]int)
    kn := iPow(r.al.Length(), k)
    c1 := make([]byte, k)
    c2 := make([]byte, k)
    c3 := make([]byte, k)
    c4 := make([]byte, k)
    cr := make([]byte, k)
    for i := 0; i < kn; i++ {
      // convert index to sequence
      for j, ix := 0, i; j < k; j++ {
        c1[k-j-1] = byte(ix % r.al.Length())
        ix        = ix / r.al.Length()
        if x, err := r.al.ComplementCoded(c1[k-j-1]); err != nil {
          return r, err
        } else {
          c2[k-j-1] = x
        }
      }
      // do not allow gaps at the ends
      if x, _ := al.Decode(c1[  0]); r.al.IsWildcard(x) {
        continue
      }
      if x, _ := al.Decode(c1[k-1]); r.al.IsWildcard(x) {
        continue
      }
      // compute indices
      i_c   := 0 // index of complement
      i_r   := 0 // index of reverse
      i_rc  := 0 // index of reverse complement
      i_res := i // final index
      for j := 0; j < k; j++ {
        i_c  += int(c2[k-j-1]) * p[j]
        i_r  += int(c1[    j]) * p[j]
        i_rc += int(c2[    j]) * p[j]
      }
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
        if err := r.decode(c1, c1); err != nil {
          return r, err
        }
        r.addEquivalentKmers(cr, c1, idx)
        // compute strings
        if comp || rc {
          if err := r.decode(c2, c2); err != nil {
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

func (obj KmersCounter) decode(dest, src []byte) error {
  for j := 0; j < len(src); j++ {
    if x, err := obj.al.Decode(src[j]); err != nil {
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
