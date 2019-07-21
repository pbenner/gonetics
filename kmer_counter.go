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
  kmap      []map[int][]int // for each k, map k-mer [ID] (with no
                            // ambiguous characters) to matching k-mers [ID]
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

func (obj *KmerCounter) generateMatchingKmersRec(dest, src []byte, m map[int]struct{}, i int) {
  if i == len(src) {
    id := obj.KmerCatalogue.GetId(string(dest))
    m[id] = struct{}{}
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
  for _, kmer := range obj.KmerCatalogue.GetNames(string(c)) {
    obj.generateMatchingKmers(d, []byte(kmer), m)
  }
  for id, _ := range m {
    i = append(i, id)
  }
  obj.kmap[k-obj.n][id] = i
  return i
}

func (obj *KmerCounter) matchingKmers(c []byte) []int {
  k := len(c)
  i := obj.KmerCatalogue.GetId(string(c))
  if r, ok := obj.kmap[k-obj.n][i]; ok {
    return r
  } else {
    if obj.frozen {
      return []int(nil)
    } else {
      return obj.addKmer(c, i)
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
  kmers  := make(KmerList, len(r))
  counts := make(map[KmerClass]int)
  i      := 0
  for id, c := range r {
    kmers[i] = KmerClass{K: k, I: id, Name: obj.KmerCatalogue.IdToName(k, id) }
    counts[kmers[i]] = c
    i++
  }
  kmers.Sort()
  return KmerCounts{Kmers: kmers, Counts: counts}
}

func (obj *KmerCounter) CountKmers(sequence []byte) KmerCounts {
  kmers  := KmerList{}
  counts := make(map[KmerClass]int)
  for k := obj.n; k <= obj.m; k++ {
    kmerCounts := obj.countKmers(sequence, k)
    kmers = append(kmers, kmerCounts.Kmers...)
    for kmer, c := range kmerCounts.Counts {
      counts[kmer] = c
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
  kmers  := make(KmerList, len(r))
  counts := make(map[KmerClass]int)
  i      := 0
  for id, _ := range r {
    kmers[i] = KmerClass{K: k, I: id, Name: obj.KmerCatalogue.IdToName(k, id) }
    counts[kmers[i]] = 1
    i++
  }
  kmers.Sort()
  return KmerCounts{Kmers: kmers, Counts: counts}
}

func (obj *KmerCounter) IdentifyKmers(sequence []byte) KmerCounts {
  kmers  := KmerList{}
  counts := make(map[KmerClass]int)
  for k := obj.n; k <= obj.m; k++ {
    kmerCounts := obj.identifyKmers(sequence, k)
    kmers = append(kmers, kmerCounts.Kmers...)
    for kmer, c := range kmerCounts.Counts {
      counts[kmer] = c
    }
  }
  return KmerCounts{Kmers: kmers, Counts: counts}
}

/* -------------------------------------------------------------------------- */

func (obj *KmerCounter) MinKmerSize() int {
  return obj.n
}

func (obj *KmerCounter) MaxKmerSize() int {
  return obj.m
}

func (obj *KmerCounter) MaxAmbiguous() []int {
  return obj.ma
}

func (obj *KmerCounter) Complement() bool {
  return obj.complement
}

func (obj *KmerCounter) Reverse() bool {
  return obj.reverse
}

func (obj *KmerCounter) Revcomp() bool {
  return obj.revcomp
}

func (obj *KmerCounter) Alphabet() ComplementableAlphabet {
  return obj.al
}

func (obj *KmerCounter) Freeze() {
  obj.frozen = true
}
