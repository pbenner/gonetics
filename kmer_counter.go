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

func NewKmerCounter(n, m int, comp, rev, rc bool, maxAmbiguous []int, al ComplementableAlphabet, kmers ...KmerClass) (*KmerCounter, error) {
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
  if len(kmers) > 0 {
  }
  return &r, nil
}

/* -------------------------------------------------------------------------- */

func (obj *KmerCounter) Clone() *KmerCounter {
  r := KmerCounter{}
  r.KmerCatalogue = *obj.KmerCatalogue.Clone()
  r.kmap = make([]map[int][]int, len(obj.kmap))
  for i := 0; i < len(obj.kmap); i++ {
    r.kmap[i] = make(map[int][]int)
    for k, v := range obj.kmap[i] {
      r.kmap[i][k] = v
    }
  }
  return &r
}

/* -------------------------------------------------------------------------- */

func (obj *KmerCounter) generalizeKmerRec(dest, src []byte, m map[int]struct{}, i int) {
  if i == len(src) {
    r := obj.KmerCatalogue.GetKmerClass(string(dest))
    m[r.I] = struct{}{}
  } else {
    x, _ := obj.Alphabet.Matching(src[i])
    for _, k := range x {
      if ok, _ := obj.Alphabet.IsWildcard(k); ok && (i == 0 || i == len(src)-1) {
        continue
      }
      dest[i] = k
      obj.generalizeKmerRec(dest, src, m, i+1)
    }
  }
}

func (obj *KmerCounter) generalizeKmer(dest, src []byte, m map[int]struct{}) {
  obj.generalizeKmerRec(dest, src, m, 0)
}

/* -------------------------------------------------------------------------- */

func (obj *KmerCounter) instantiateKmerRec(dest, src []byte, m map[int]struct{}, i int) {
  if i == len(src) {
    // compute equivalence class but do not store it in
    // the kmer catalogue!
    r := obj.KmerCatalogue.EquivalenceClass(string(dest))
    m[r.I] = struct{}{}
  } else {
    x, _ := obj.Alphabet.Bases(src[i])
    for _, k := range x {
      dest[i] = k
      obj.instantiateKmerRec(dest, src, m, i+1)
    }
  }
}

func (obj *KmerCounter) instantiateKmer(dest, src []byte, m map[int]struct{}) {
  obj.instantiateKmerRec(dest, src, m, 0)
}

/* -------------------------------------------------------------------------- */

func (obj *KmerCounter) addObservedKmer(kmer KmerClass) []int {
  d := make([]byte, kmer.K)
  m := make(map[int]struct{})
  i := []int{}
  for _, kmer := range kmer.Elements {
    obj.generalizeKmer(d, []byte(kmer), m)
  }
  for id, _ := range m {
    i = append(i, id)
  }
  obj.kmap[kmer.K-obj.N][kmer.I] = i
  return i
}

func (obj *KmerCounter) addKmer(kmer KmerClass) []int {
  d := make([]byte, kmer.K)
  m := make(map[int]struct{})
  i := []int{}
  for _, kmer := range kmer.Elements {
    obj.instantiateKmer(d, []byte(kmer), m)
  }
  for id, _ := range m {
    obj.kmap[kmer.K-obj.N][kmer.I] = append(obj.kmap[kmer.K-obj.N][kmer.I], id)
  }
  return i
}

func (obj *KmerCounter) matchingKmers(c []byte) []int {
  r := obj.KmerCatalogue.GetKmerClass(string(c))
  if i, ok := obj.kmap[r.K-obj.N][r.I]; ok {
    return i
  } else {
    if obj.frozen {
      return []int(nil)
    } else {
      return obj.addObservedKmer(r)
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
  for k := obj.N; k <= obj.M; k++ {
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
  for k := obj.N; k <= obj.M; k++ {
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
