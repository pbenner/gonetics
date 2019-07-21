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

type KmerCountsList struct {
  Kmers    KmerList
  Counts []map[KmerClass]int
}

func NewKmerCountsList(counts ...KmerCounts) KmerCountsList {
  r := KmerCountsList{}
  return r.Append(counts...)
}

func (obj KmerCountsList) Append(args ...KmerCounts) KmerCountsList {
  if len(args) == 0 {
    return obj
  }
  idLists := make([]KmerList, len(args))
  counts  := obj.Counts
  for i, c := range args {
    idLists[i] = c.Kmers
    counts = append(counts, c.Counts)
  }
  ids := obj.Kmers.Union(idLists...)
  return KmerCountsList{Kmers: ids, Counts: counts}
}

func (obj KmerCountsList) Len() int {
  return len(obj.Counts)
}

func (obj KmerCountsList) At(i int) KmerCounts {
  return KmerCounts{Kmers: obj.Kmers, Counts: obj.Counts[i]}
}

func (obj *KmerCountsList) Slice(i, j int) KmerCountsList {
  return KmerCountsList{Kmers: obj.Kmers, Counts: obj.Counts[i:j]}
}

/* -------------------------------------------------------------------------- */

type KmerCountsIterator struct {
  KmerCounts
  i int
}

func (obj KmerCountsIterator) Ok() bool {
  return obj.i < obj.Len()
}

func (obj KmerCountsIterator) GetKmer() KmerClass {
  return obj.Kmers[obj.i]
}

func (obj KmerCountsIterator) GetName() string {
  return obj.Kmers[obj.i].Name
}

func (obj KmerCountsIterator) GetCount() int {
  return obj.At(obj.i)
}

func (obj *KmerCountsIterator) Next() {
  obj.i++
}
