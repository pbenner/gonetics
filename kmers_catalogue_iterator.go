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

type KmersCatalogueIterator struct {
  names []string
  ids   []int
  i       int
}

/* -------------------------------------------------------------------------- */

func NewKmersCatalogueIterator(kmersSet KmersCatalogue) KmersCatalogueIterator {
  names := []string{}
  ids   := []int{}
  for k := 0; k < len(kmersSet.names); k++ {
    n := []string{}
    i := []int{}
    for id, name := range kmersSet.names[k] {
      n = append(n, name)
      i = append(i, id)
    }
    sortIntStringPairs{i, n}.Sort()
    names = append(names, n...)
    ids   = append(ids,   i...)
  }
  return KmersCatalogueIterator{names: names, ids: ids, i: 0}
}

/* -------------------------------------------------------------------------- */

func (obj KmersCatalogueIterator) Ok() bool {
  return obj.i < len(obj.names)
}

func (obj KmersCatalogueIterator) Get() string {
  return obj.names[obj.i]
}

func (obj KmersCatalogueIterator) GetId() int {
  return obj.ids[obj.i]
}

func (obj *KmersCatalogueIterator) Next() {
  obj.i++
}
