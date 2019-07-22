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
import "sort"

/* -------------------------------------------------------------------------- */

type KmerCatalogue struct {
  KmerEquivalenceRelation
  elements []map[int][]string   // map k-mer class ID (k,i) to class elements
  idmap    []map[string]int     // map k-mer (class instances) to class ID (k,i)
}

/* -------------------------------------------------------------------------- */

func NewKmerCatalogue(n, m int, comp, rev, rc bool, maxAmbiguous []int, al ComplementableAlphabet) (*KmerCatalogue, error) {
  r := KmerCatalogue{}
  if f, err := NewKmerEquivalenceRelation(n, m, comp, rev, rc, maxAmbiguous, al); err != nil {
    return nil, err
  } else {
    r.KmerEquivalenceRelation = f
  }
  idmap    := make([]map[string]int  , m-n+1)
  elements := make([]map[int][]string, m-n+1)
  for k := n; k <= m; k++ {
    idmap   [k-n] = make(map[string]int)
    elements[k-n] = make(map[int][]string)
  }
  r.elements = elements
  r.idmap    = idmap
  return &r, nil
}

/* -------------------------------------------------------------------------- */

func (obj *KmerCatalogue) GetKmerClassIfPresent(kmer string) (KmerClass, bool) {
  k := len(kmer)
  if k < obj.n || k > obj.m {
    panic("k-mer has invalid length")
  }
  if i, ok := obj.idmap[k-obj.n][kmer]; ok {
    return NewKmerClass(k, i, obj.elements[k-obj.n][i]), true
  } else {
    return KmerClass{}, false
  }
}

func (obj *KmerCatalogue) GetKmerClass(kmer string) KmerClass {
  k := len(kmer)
  if k < obj.n || k > obj.m {
    panic("k-mer has invalid length")
  }
  if i, ok := obj.idmap[k-obj.n][kmer]; ok {
    return NewKmerClass(k, i, obj.elements[k-obj.n][i])
  } else {
    r := obj.EquivalenceClass(kmer)
    for _, kmer := range r.Elements {
      obj.idmap[k-obj.n][kmer] = r.I
    }
    obj.elements[k-obj.n][r.I] = r.Elements
    return r
  }
}

func (obj *KmerCatalogue) GetKmerClassFromId(k, id int) KmerClass {
  if elements, ok := obj.elements[k-obj.n][id]; ok {
    return NewKmerClass(k, id, elements)
  }
  panic("k-mer name not found")
}

func (obj *KmerCatalogue) CatalogueSize() int {
  r := 0
  for i := 0; i < len(obj.elements); i++ {
    r += len(obj.elements[i])
  }
  return r
}

/* TODO: extend code below to k-mer equivalence graph
 * -------------------------------------------------------------------------- */

func (obj *KmerCatalogue) scanSubKmers_(kmer []byte, k int) []int {
  idMap := make(map[int]struct{})
  // loop over sequence
  for i := 0; i < len(kmer); i++ {
    if i+k-1 >= len(kmer) {
      break
    }
    it := NewKmerInstantiationIterator(obj.al, string(kmer[i:i+k]), true)
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

func (obj *KmerCatalogue) scanSubKmers(kmer []byte) []string {
  names := []string{}
  for k := obj.n; k <= obj.m; k++ {
    ids := obj.scanSubKmers_(kmer, k)
    for _, id := range ids {
      r := obj.GetKmerClassFromId(k, id)
      names = append(names , r.String())
    }
  }
  return names
}

/* -------------------------------------------------------------------------- */

func (obj *KmerCatalogue) relatedKmers(s []byte, m, k int) []int {
  idMap := make(map[int]struct{})
  // loop over positions where the k-mer can be fixed
  for j := 0; j <= k - len(s); j++ {
    for it := NewKmerCylinderIterator(k, obj.ma[k-obj.n] - m, obj.al, j, string(s)); it.Ok(); it.Next() {
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

func (obj *KmerCatalogue) RelatedKmers(kmer string) []string {
  s := []byte(obj.GetKmerClass(kmer).Elements[0])
  m := obj.countAmbiguous(kmer)
  // scan k-mer for sub-k-mers
  names := obj.scanSubKmers(s)
  // loop over k-mer sizes
  for k := len(s)+1; k <= obj.m; k++ {
    ids := obj.relatedKmers(s, m, k)
    for _, id := range ids {
      r := obj.GetKmerClassFromId(k, id)
      names = append(names , r.String())
    }
  }
  return names
}

func (obj *KmerCatalogue) countAmbiguous(kmer string) int {
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
