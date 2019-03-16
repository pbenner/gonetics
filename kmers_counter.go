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

type KmersCounter struct {
  n, m        int      // min and max kmer size
  length      int      // code length
  indices [][]int      // index map
  names     []string   // kmer names
  p         []int      // pre-evaluated powers
  complement  bool
  reverse     bool
  revcomp     bool
  al          NucleotideAlphabet
}

/* -------------------------------------------------------------------------- */

func NewKmersCounter(n, m int, comp, rev, rc bool) (KmersCounter, error) {
  r := KmersCounter{n: n, m: m, complement: comp, reverse: rev, revcomp: rc}
  p := make([]int, m+1)
  for k := 0; k <= m; k++ {
    p[k] = iPow(r.al.Length(), k)
  }
  names   := []string{}
  indices := make([][]int, m-n+1)
  idx     := 0
  for k := n; k <= m; k++ {
    kn := iPow(r.al.Length(), k)
    indices[k-n] = make([]int, kn)
    c1 := make([]byte, k)
    c2 := make([]byte, k)
    c3 := make([]byte, k)
    c4 := make([]byte, k)
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
      // compute indices
      i_c   := 0 // index of complement
      i_r   := 0 // index of reverse
      i_rc  := 0 // index of reverse complement
      i_res := i // final index
      for j := 0; j < k; j++ {
        i_c  += int(c2[    j]) * p[j]
        i_r  += int(c1[k-j-1]) * p[j]
        i_rc += int(c2[k-j-1]) * p[j]
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
      if i_res != i {
        indices[k-n][i] = indices[k-n][i_res]
      } else {
        // new kmer found
        indices[k-n][i] = idx; idx += 1
        // compute strings
        if err := r.decode(c1, c1); err != nil {
          return r, err
        }
        if comp || rc {
          if err := r.decode(c2, c2); err != nil {
            return r, err
          }
        }
        name := string(c1)
        if comp {
          name = fmt.Sprintf("%s|%s", name, string(c2))
        }
        if rev {
          r.rev(c3, c1)
          name = fmt.Sprintf("%s|%s", name, string(c3))
        }
        if rc {
          r.rev(c4, c2)
          name = fmt.Sprintf("%s|%s", name, string(c4))
        }
        names = append(names, name)
      }
    }
  }
  r.length  = idx
  r.p       = p
  r.names   = names
  r.indices = indices
  return r, nil
}

/* -------------------------------------------------------------------------- */

func (obj KmersCounter) CountKmers(result []int, sequence []byte) error {
  c := make([]int, len(sequence))
  for i := 0; i < len(sequence); i++ {
    if sequence[i] == 'n' || sequence[i] == 'N' {
      c[i] = -1
      continue
    }
    if r, err := obj.al.Code(sequence[i]); err != nil {
      return err
    } else {
      c[i] = int(r)
    }
  }
  // loop over sequence
  for i := 0; i < len(sequence); i++ {
    // loop over all k-mers
kLoop:
    for k := obj.n; k <= obj.m && i+k-1 < len(c); k++ {
      // eval k-mer
      s := 0
      for j := 0; j < k; j++ {
        if c[i+j] == -1 {
          break kLoop
        }
        s += c[i+j] * obj.p[k-j-1]
      }
      result[obj.index(k, s)] += 1
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

/* -------------------------------------------------------------------------- */

func (obj KmersCounter) index(k, i int) int {
  return obj.indices[k-obj.n][i]
}

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
  for i, j := 0, len(src)-1; i < j; i, j = i+1, j-1 {
    dest[i], dest[j] = src[j], src[i]
  }
}
