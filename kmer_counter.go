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
import "strings"

/* -------------------------------------------------------------------------- */

type KmerCounter struct {
  KmerCatalogue
  // map of observed k-mers (i.e. k-mers without any ambiguous characters)
  // to matching k-mers
  kmap      []map[int][]int
  frozen      bool
}

/* -------------------------------------------------------------------------- */

func NewKmerCounter(n, m int, comp, rev, rc bool, maxAmbiguous []int, al ComplementableAlphabet) (*KmerCounter, error) {
  r := KmerCounter{}
  if s, err := NewKmerCatalogue(n, m, comp, rev, rc, maxAmbiguous, al); err != nil {
    return nil, err
  } else {
    r.KmerCatalogue = *s
  }
  r.kmap = make([]map[int][]int, m-n+1)
  for k := n; k <= m; k++ {
    r.kmap[k-n] = make(map[int][]int)
  }
  return &r, nil
}

/* -------------------------------------------------------------------------- */

func (obj *KmerCounter) Clone() *KmerCounter {
  r := KmerCounter{}
  r.KmerCatalogue = *obj.KmerCatalogue.Clone()
  r.kmap = make([]map[int][]int, len(obj.kmap))
  for i := 0; i < len(obj.kmap); i++ {
    for k, v := range obj.kmap[i] {
      r.kmap[i][k] = v
    }
  }
  return &r
}

/* -------------------------------------------------------------------------- */

func (obj *KmerCounter) generateMatchingKmersRec(dest, src []byte, m map[int]struct{}, i int) {
  if i == len(src) {
    r := obj.KmerCatalogue.GetKmerClass(string(dest))
    m[r.I] = struct{}{}
  } else {
    x, _ := obj.al.Matching(src[i])
    for _, k := range x {
      if ok, _ := obj.al.IsWildcard(k); ok && (i == 0 || i == len(src)-1) {
        continue
      }
      dest[i] = k
      obj.generateMatchingKmersRec(dest, src, m, i+1)
    }
  }
}

func (obj *KmerCounter) generateMatchingKmers(dest, src []byte, m map[int]struct{}) {
  obj.generateMatchingKmersRec(dest, src, m, 0)
}

/* -------------------------------------------------------------------------- */

func (obj *KmerCounter) addKmer(c []byte, id int) []int {
  k := len(c)
  d := make([]byte, k)
  m := make(map[int]struct{})
  i := []int{}
  for _, kmer := range obj.KmerCatalogue.GetKmerClass(string(c)).Elements {
    obj.generateMatchingKmers(d, []byte(kmer), m)
  }
  for id, _ := range m {
    i = append(i, id)
  }
  obj.kmap[k-obj.n][id] = i
  return i
}

func (obj *KmerCounter) matchingKmers(c []byte) []int {
  r := obj.KmerCatalogue.GetKmerClass(string(c))
  if i, ok := obj.kmap[r.K-obj.n][r.I]; ok {
    return i
  } else {
    if obj.frozen {
      return []int(nil)
    } else {
      return obj.addKmer(c, r.I)
    }
  }
}

/* -------------------------------------------------------------------------- */

func (obj *KmerCounter) countKmers(sequence []byte, k int) KmerCounts {
  c := []byte(strings.ToLower(string(sequence)))
  r := make(map[int]int)
  // loop over sequence
  for i := 0; i < len(c); i++ {
    if i+k-1 < len(c) {
      for _, j := range obj.matchingKmers(c[i:i+k]) {
        r[j] += 1
      }
    }
  }
  // construct a list of k-mers
  kmers  := make(KmerClassList, len(r))
  counts := make(map[KmerClassId]int)
  i      := 0
  for id, c := range r {
    kmers[i] = obj.KmerCatalogue.GetKmerClassFromId(k, id)
    counts[kmers[i].KmerClassId] = c
    i++
  }
  kmers.Sort()
  return KmerCounts{Kmers: kmers, Counts: counts}
}

func (obj *KmerCounter) CountKmers(sequence []byte) KmerCounts {
  kmers  := KmerClassList{}
  counts := make(map[KmerClassId]int)
  for k := obj.n; k <= obj.m; k++ {
    kmerCounts := obj.countKmers(sequence, k)
    kmers = append(kmers, kmerCounts.Kmers...)
    for kmerId, c := range kmerCounts.Counts {
      counts[kmerId] = c
    }
  }
  return KmerCounts{Kmers: kmers, Counts: counts}
}

func (obj *KmerCounter) identifyKmers(sequence []byte, k int) KmerCounts {
  c := []byte(strings.ToLower(string(sequence)))
  r := make(map[int]struct{})
  // loop over sequence
  for i := 0; i < len(c); i++ {
    // loop over all k-mers
    if i+k-1 < len(c) {
      for _, j := range obj.matchingKmers(c[i:i+k]) {
        r[j] = struct{}{}
      }
    }
  }
  // construct a list of k-mers
  kmers  := make(KmerClassList, len(r))
  counts := make(map[KmerClassId]int)
  i      := 0
  for id, _ := range r {
    kmers[i] = obj.KmerCatalogue.GetKmerClassFromId(k, id)
    counts[kmers[i].KmerClassId] = 1
    i++
  }
  kmers.Sort()
  return KmerCounts{Kmers: kmers, Counts: counts}
}

func (obj *KmerCounter) IdentifyKmers(sequence []byte) KmerCounts {
  kmers  := KmerClassList{}
  counts := make(map[KmerClassId]int)
  for k := obj.n; k <= obj.m; k++ {
    kmerCounts := obj.identifyKmers(sequence, k)
    kmers = append(kmers, kmerCounts.Kmers...)
    for kmerId, c := range kmerCounts.Counts {
      counts[kmerId] = c
    }
  }
  return KmerCounts{Kmers: kmers, Counts: counts}
}

/* -------------------------------------------------------------------------- */

func (obj *KmerCounter) Freeze() {
  obj.frozen = true
}
