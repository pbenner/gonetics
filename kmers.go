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

type KmerId struct {
  K    int
  I    int
  Name string
}

/* -------------------------------------------------------------------------- */

type KmerIdList []KmerId

func (obj KmerIdList) Len() int {
  return len(obj)
}

func (obj KmerIdList) Less(i, j int) bool {
  if obj[i].K != obj[j].K {
    return obj[i].K < obj[j].K
  } else {
    return obj[i].I < obj[j].I
  }
}

func (obj KmerIdList) Swap(i, j int) {
  obj[i].K, obj[j].K = obj[j].K, obj[i].K
  obj[i].I, obj[j].I = obj[j].I, obj[i].I
}

func (obj KmerIdList) Sort() {
  sort.Sort(obj)
}

func (obj KmerIdList) Union(b ...KmerIdList) KmerIdList {
  m := make(KmerIdSet)
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

type KmerIdSet map[KmerId]struct{}

func (obj KmerIdSet) AsList() KmerIdList {
  r := make(KmerIdList, len(obj))
  i := 0
  for id, _ := range obj {
    r[i] = id; i++
  }
  r.Sort()
  return r
}

/* -------------------------------------------------------------------------- */

type KmerCounts struct {
  Ids    KmerIdList
  Counts map[KmerId]int
}

func (obj KmerCounts) Len() int {
  return len(obj.Ids)
}

func (obj KmerCounts) At(i int) int {
  if c, ok := obj.Counts[obj.Ids[i]]; ok {
    return c
  } else {
    return 0
  }
}

/* -------------------------------------------------------------------------- */

type KmerCountsList struct {
  Ids      KmerIdList
  Counts []map[KmerId]int
}

func (obj *KmerCountsList) Append(r ...KmerCounts) KmerCountsList {
  idLists := make([]KmerIdList, len(r))
  for i, ri := range r {
    idLists[i] = ri.Ids
  }
  ids    := obj.Ids.Union(idLists...)
  counts := obj.Counts
  for _, ri := range r {
    counts = append(counts, ri.Counts)
  }
  return KmerCountsList{Ids: ids, Counts: counts}
}

func (obj *KmerCountsList) Len() int {
  return len(obj.Counts)
}

func (obj *KmerCountsList) At(i int) KmerCounts {
  return KmerCounts{Ids: obj.Ids, Counts: obj.Counts[i]}
}

/* -------------------------------------------------------------------------- */

type KmerCountsIterator struct {
  KmerCounts
  i int
}

func (obj KmerCountsIterator) Ok() bool {
  return obj.i < obj.Len()
}

func (obj KmerCountsIterator) GetId() KmerId {
  return obj.Ids[obj.i]
}

func (obj KmerCountsIterator) GetName() string {
  return obj.Ids[obj.i].Name
}

func (obj KmerCountsIterator) GetCount() int {
  return obj.At(obj.i)
}
