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

type KmerCatalogueIterator struct {
  kmers KmerClassList
  ids   []int
  i       int
}

/* -------------------------------------------------------------------------- */

func NewKmerCatalogueIterator(catalogue KmerCatalogue) KmerCatalogueIterator {
  kmers := KmerClassList{}
  for k := 0; k < len(catalogue.elements); k++ {
    for i, elements := range catalogue.elements[k] {
      kmers = append(kmers, NewKmerClass(k, i, elements))
    }
  }
  kmers.Sort()
  return KmerCatalogueIterator{kmers: kmers, i: 0}
}

/* -------------------------------------------------------------------------- */

func (obj KmerCatalogueIterator) Ok() bool {
  return obj.i < len(obj.kmers)
}

func (obj KmerCatalogueIterator) Get() KmerClass {
  return obj.kmers[obj.i]
}

func (obj *KmerCatalogueIterator) Next() {
  obj.i++
}
