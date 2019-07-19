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

type KmersCatalogue struct {
  n, m        int              // min and max kmer size
  complement  bool
  reverse     bool
  revcomp     bool
  ma        []int              // maximum number of ambiguous letters
  al          ComplementableAlphabet
  p         []int              // pre-evaluated powers
  names     []map[int]string   // k-mer names (equivalent k-mers separated by pipe)
  idmap     []map[string]int   // unique k-mer IDs
}

/* -------------------------------------------------------------------------- */

func NewKmersCatalogue(n, m int, comp, rev, rc bool, maxAmbiguous []int, al ComplementableAlphabet) (*KmersCatalogue, error) {
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
    return nil, fmt.Errorf("parameter `maxAmbiguous' has invalid length")
  }
  p := make([]int, m+1)
  for k := 0; k <= m; k++ {
    p[k] = iPow(al.Length(), k)
  }
  idmap := make([]map[string]int, m-n+1)
  names := make([]map[int]string, m-n+1)
  for k := n; k <= m; k++ {
    idmap[k-n] = make(map[string]int)
    names[k-n] = make(map[int]string)
  }
  r := KmersCatalogue{
    n         : n,
    m         : m,
    p         : p,
    complement: comp,
    reverse   : rev,
    revcomp   : rc,
    ma        : maxAmbiguous,
    idmap     : idmap,
    names     : names,
    al        : al }
  return &r, nil
}

/* -------------------------------------------------------------------------- */

func (obj *KmersCatalogue) MinKmerSize() int {
  return obj.n
}

func (obj *KmersCatalogue) MaxKmerSize() int {
  return obj.m
}

func (obj *KmersCatalogue) MaxAmbiguous() []int {
  return obj.ma
}

func (obj *KmersCatalogue) Complement() bool {
  return obj.complement
}

func (obj *KmersCatalogue) Reverse() bool {
  return obj.reverse
}

func (obj *KmersCatalogue) Revcomp() bool {
  return obj.revcomp
}

func (obj *KmersCatalogue) Alphabet() ComplementableAlphabet {
  return obj.al
}

/* -------------------------------------------------------------------------- */

func (obj *KmersCatalogue) GetId(kmer string) int {
  k := len(kmer)
  if k < obj.n || k > obj.m {
    panic("k-mer has invalid length")
  }
  if i, ok := obj.idmap[k-obj.n][kmer]; ok {
    return i
  } else {
    i, _ = obj.computeId([]byte(kmer))
    return i
  }
}

func (obj *KmersCatalogue) IdToName(k, id int) string {
  if name, ok := obj.names[k-obj.n][id]; ok {
    return name
  }
  panic("k-mer name not found")
}

func (obj *KmersCatalogue) GetName(kmer string) string {
  k := len(kmer)
  if k < obj.n || k > obj.m {
    panic("k-mer has invalid length")
  }
  if i, ok := obj.idmap[k-obj.n][kmer]; ok {
    return obj.names[k-obj.n][i]
  } else {
    _, s := obj.computeId([]byte(kmer))
    return s
  }
}

func (obj *KmersCatalogue) GetNames(kmer string) []string {
  return strings.Split(obj.GetName(kmer), "|")
}

func (obj *KmersCatalogue) ObservedKmers() int {
  r := 0
  for i := 0; i < len(obj.names); i++ {
    r += len(obj.names[i])
  }
  return r
}

/* -------------------------------------------------------------------------- */

func (obj *KmersCatalogue) computeId(c1 []byte) (int, string) {
  k  := len(c1)
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
  obj.idmap[k-obj.n][name] = i
  if obj.Complement() {
    name = fmt.Sprintf("%s|%s", name, string(c2))
    obj.idmap[k-obj.n][string(c2)] = i
  }
  if obj.Reverse() {
    name = fmt.Sprintf("%s|%s", name, string(c3))
    obj.idmap[k-obj.n][string(c3)] = i
  }
  if obj.Revcomp() {
    name = fmt.Sprintf("%s|%s", name, string(c4))
    obj.idmap[k-obj.n][string(c4)] = i
  }
  obj.names[k-obj.n][i] = name
  return i, name
}

/* -------------------------------------------------------------------------- */

func (obj *KmersCatalogue) comp(dest, src []byte) error {
  for j := 0; j < len(src); j++ {
    if x, err := obj.al.Complement(src[j]); err != nil {
      return err
    } else {
      dest[j] = x
    }
  }
  return nil
}

func (obj *KmersCatalogue) rev(dest, src []byte) {
  for i, j := 0, len(src)-1; i <= j; i, j = i+1, j-1 {
    dest[i], dest[j] = src[j], src[i]
  }
}

/* -------------------------------------------------------------------------- */

func (obj *KmersCatalogue) scanSubKmers_(kmer []byte, k int) []int {
  idMap := make(map[int]struct{})
  // loop over sequence
  for i := 0; i < len(kmer); i++ {
    if i+k-1 >= len(kmer) {
      break
    }
    it := NewKmersInstantiationIterator(obj.al, string(kmer[i:i+k]), true)
    for ; it.Ok(); it.Next() {
      if id, ok := obj.idmap[k-obj.n][it.Get()]; ok {
        idMap[id] = struct{}{}
      }
    }
  }
  ids := []int{}
  for id, _ := range idMap {
    ids = append(ids, id)
  }
  sort.Ints(ids)
  return ids
}

func (obj *KmersCatalogue) scanSubKmers(kmer []byte) []string {
  names := []string{}
  for k := obj.n; k <= obj.m; k++ {
    ids := obj.scanSubKmers_(kmer, k)
    for _, id := range ids {
      names = append(names , obj.IdToName(k, id))
    }
  }
  return names
}

/* -------------------------------------------------------------------------- */

func (obj *KmersCatalogue) relatedKmers(s []byte, m, k int) []int {
  idMap := make(map[int]struct{})
  // loop over positions where the k-mer can be fixed
  for j := 0; j <= k - len(s); j++ {
    for it := NewKmersCylinderIterator(k, obj.ma[k-obj.n] - m, obj.al, j, string(s)); it.Ok(); it.Next() {
      if id, ok := obj.idmap[k-obj.n][it.Get()]; ok {
        idMap[id] = struct{}{}
      }
    }
  }
  ids := []int{}
  for id, _ := range idMap {
    ids = append(ids, id)
  }
  sort.Ints(ids)
  return ids
}

func (obj *KmersCatalogue) RelatedKmers(kmer string) []string {
  s := []byte(obj.GetNames(kmer)[0])
  m := obj.countAmbiguous(kmer)
  // scan k-mer for sub-k-mers
  names := obj.scanSubKmers(s)
  // loop over k-mer sizes
  for k := len(s)+1; k <= obj.m; k++ {
    ids := obj.relatedKmers(s, m, k)
    for _, id := range ids {
      names = append(names , obj.IdToName(k, id))
    }
  }
  return names
}

func (obj *KmersCatalogue) countAmbiguous(kmer string) int {
  m := 0
  s := []byte(kmer)
  for i := 0; i < len(s); i++ {
    if ok, err := obj.al.IsAmbiguous(s[i]); err != nil {
      panic("internal error")
    } else {
      if ok {
        m++
      }
    }
  }
  return m
}
