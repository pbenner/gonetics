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

/* K-mer equivalence class
 * -------------------------------------------------------------------------- */

type KmerClass struct {
  // K-mer ID
  K    int
  I    int
  // K-mer string representation
  Name string
}

func (obj KmerClass) Equals(b KmerClass) bool {
  if obj.K != b.K {
    return false
  }
  if obj.I != b.I {
    return false
  }
  return true
}

/* -------------------------------------------------------------------------- */

type KmerList []KmerClass

func (obj KmerList) Clone() KmerList {
  r := make(KmerList, len(obj))
  copy(r, obj)
  return r
}

func (obj KmerList) Equals(b KmerList) bool {
  if len(obj) != len(b) {
    return false
  }
  for i := 0; i < len(obj); i++ {
    if !obj[i].Equals(b[i]) {
      return false
    }
  }
  return true
}

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

type KmerSet map[KmerClass]struct{}

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
  Counts map[KmerClass]int
}

func (obj KmerCounts) Len() int {
  return len(obj.Kmers)
}

func (obj KmerCounts) N() int {
  return len(obj.Counts)
}

func (obj KmerCounts) At(i int) int {
  if c, ok := obj.Counts[obj.Kmers[i]]; ok {
    return c
  } else {
    return 0
  }
}

func (obj KmerCounts) GetCount(kmer KmerClass) int {
  if c, ok := obj.Counts[kmer]; ok {
    return c
  } else {
    return 0
  }
}

func (obj KmerCounts) GetKmer(i int) KmerClass {
  return obj.Kmers[i]
}

func (obj KmerCounts) Iterate() KmerCountsIterator {
  return KmerCountsIterator{obj, 0}
}
