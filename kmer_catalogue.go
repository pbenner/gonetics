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

type KmerCatalogue struct {
  KmerEquivalenceRelation
  elements []map[int][]string   // map k-mer class ID (k,i) to class elements
  idmap    []map[string]int     // map k-mer (class instances) to class ID (k,i)
}

/* -------------------------------------------------------------------------- */

func NewKmerCatalogue(n, m int, comp, rev, rc bool, maxAmbiguous []int, al ComplementableAlphabet) (*KmerCatalogue, error) {
  if rel, err := NewKmerEquivalenceRelation(n, m, comp, rev, rc, maxAmbiguous, al); err != nil {
    return nil, err
  } else {
    return newKmerCatalogue(rel), nil
  }
}

func newKmerCatalogue(rel KmerEquivalenceRelation) *KmerCatalogue {
  idmap    := make([]map[string]int  , rel.M-rel.N+1)
  elements := make([]map[int][]string, rel.M-rel.N+1)
  for k := rel.N; k <= rel.M; k++ {
    idmap   [k-rel.N] = make(map[string]int)
    elements[k-rel.N] = make(map[int][]string)
  }
  r := KmerCatalogue{}
  r.KmerEquivalenceRelation = rel
  r.elements = elements
  r.idmap    = idmap
  return &r
}

/* -------------------------------------------------------------------------- */

func (obj *KmerCatalogue) Clone() *KmerCatalogue {
  r := KmerCatalogue{}
  r.KmerEquivalenceRelation = obj.KmerEquivalenceRelation
  r.idmap    = make([]map[string]int  , len(obj.idmap))
  r.elements = make([]map[int][]string, len(obj.elements))
  for i := 0; i < len(obj.idmap); i++ {
    r.idmap[i] = make(map[string]int)
    for k, v := range obj.idmap[i] {
      r.idmap[i][k] = v
    }
  }
  for i := 0; i < len(obj.elements); i++ {
    r.elements[i] = make(map[int][]string)
    for k, v := range obj.elements[i] {
      r.elements[i][k] = v
    }
  }
  return &r
}

/* -------------------------------------------------------------------------- */

func (obj *KmerCatalogue) GetKmerClassIfPresent(kmer string) (KmerClass, bool) {
  k := len(kmer)
  if k < obj.N || k > obj.M {
    panic("k-mer has invalid length")
  }
  if i, ok := obj.idmap[k-obj.N][kmer]; ok {
    return NewKmerClass(k, i, obj.elements[k-obj.N][i]), true
  } else {
    return KmerClass{}, false
  }
}

func (obj *KmerCatalogue) AddKmerClass(kmer KmerClass) {
  for _, s := range kmer.Elements {
    obj.idmap[kmer.K-obj.N][s] = kmer.I
  }
  obj.elements[kmer.K-obj.N][kmer.I] = kmer.Elements
}

func (obj *KmerCatalogue) GetKmerClass(kmer string) KmerClass {
  k := len(kmer)
  if k < obj.N || k > obj.M {
    panic("k-mer has invalid length")
  }
  if i, ok := obj.idmap[k-obj.N][kmer]; ok {
    return NewKmerClass(k, i, obj.elements[k-obj.N][i])
  } else {
    r := obj.EquivalenceClass(kmer)
    obj.AddKmerClass(r)
    return r
  }
}

func (obj *KmerCatalogue) GetKmerClassFromId(k, id int) KmerClass {
  if elements, ok := obj.elements[k-obj.N][id]; ok {
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
