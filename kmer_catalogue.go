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

func (obj *KmerCatalogue) Clone() *KmerCatalogue {
  r := KmerCatalogue{}
  r.KmerEquivalenceRelation = obj.KmerEquivalenceRelation
  r.idmap    = make([]map[string]int  , len(obj.idmap))
  r.elements = make([]map[int][]string, len(obj.elements))
  for i := 0; i < len(obj.idmap); i++ {
    for k, v := range obj.idmap[i] {
      r.idmap[i][k] = v
    }
  }
  for i := 0; i < len(obj.elements); i++ {
    for k, v := range obj.elements[i] {
      r.elements[i][k] = v
    }
  }
  return &r
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
