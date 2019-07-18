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

type Kmer struct {
  // K-mer ID
  K    int
  I    int
  // K-mer string representation
  Name string
}

/* -------------------------------------------------------------------------- */

type KmerList []Kmer

func (obj KmerList) Len() int {
  return len(obj)
}

func (obj KmerList) Less(i, j int) bool {
  if obj[i].K != obj[j].K {
    return obj[i].K < obj[j].K
  } else {
    return obj[i].I < obj[j].I
  }
}

func (obj KmerList) Swap(i, j int) {
  obj[i], obj[j] = obj[j], obj[i]
}

func (obj KmerList) Sort() {
  sort.Sort(obj)
}

func (obj KmerList) Union(b ...KmerList) KmerList {
  m := make(KmerSet)
  for _, id := range obj {
    m[id] = struct{}{}
  }
  for _, bi := range b {
    for _, id := range bi {
      m[id] = struct{}{}
    }
  }
  return m.AsList()
}

/* -------------------------------------------------------------------------- */

type KmerSet map[Kmer]struct{}

func (obj KmerSet) AsList() KmerList {
  r := make(KmerList, len(obj))
  i := 0
  for id, _ := range obj {
    r[i] = id; i++
  }
  r.Sort()
  return r
}

/* -------------------------------------------------------------------------- */

type KmerCounts struct {
  // this is a sorted list of k-mers and might contain more entries than
  // the counts map
  Kmers  KmerList
  Counts map[Kmer]int
}

func (obj KmerCounts) Len() int {
  return len(obj.Kmers)
}

func (obj KmerCounts) At(i int) int {
  if c, ok := obj.Counts[obj.Kmers[i]]; ok {
    return c
  } else {
    return 0
  }
}

func (obj KmerCounts) Iterate() KmerCountsIterator {
  return KmerCountsIterator{obj, 0}
}

/* -------------------------------------------------------------------------- */

type KmerCountsList struct {
  Kmers    KmerList
  Counts []map[Kmer]int
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

/* -------------------------------------------------------------------------- */

type KmerCountsIterator struct {
  KmerCounts
  i int
}

func (obj KmerCountsIterator) Ok() bool {
  return obj.i < obj.Len()
}

func (obj KmerCountsIterator) GetKmer() Kmer {
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
