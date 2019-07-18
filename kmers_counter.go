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

type KmersCounter struct {
  KmersSet
  kmap      []map[int][]int // for each k, map k-mer [ID] (with no
                            // ambiguous characters) to matching k-mers [ID]
}

/* -------------------------------------------------------------------------- */

func NewKmersCounter(n, m int, comp, rev, rc bool, maxAmbiguous []int, al ComplementableAlphabet) (*KmersCounter, error) {
  r := KmersCounter{}
  if s, err := NewKmersSet(n, m, comp, rev, rc, maxAmbiguous, al); err != nil {
    return nil, err
  } else {
    r.KmersSet = *s
  }
  r.kmap = make([]map[int][]int, m-n+1)
  for k := n; k <= m; k++ {
    r.kmap[k-n] = make(map[int][]int)
  }
  return &r, nil
}

/* -------------------------------------------------------------------------- */

func (obj *KmersCounter) generateMatchingKmersRec(dest, src []byte, m map[int]struct{}, i int) {
  if i == len(src) {
    id := obj.KmersSet.GetId(string(dest))
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

func (obj *KmersCounter) generateMatchingKmers(dest, src []byte, m map[int]struct{}) {
  obj.generateMatchingKmersRec(dest, src, m, 0)
}

/* -------------------------------------------------------------------------- */

func (obj *KmersCounter) addKmer(c []byte, id int) []int {
  k := len(c)
  d := make([]byte, k)
  m := make(map[int]struct{})
  i := []int{}
  for _, kmer := range obj.KmersSet.GetNames(string(c)) {
    obj.generateMatchingKmers(d, []byte(kmer), m)
  }
  for id, _ := range m {
    i = append(i, id)
  }
  obj.kmap[k-obj.n][id] = i
  return i
}

func (obj *KmersCounter) matchingKmers(c []byte) []int {
  k := len(c)
  i := obj.KmersSet.GetId(string(c))
  if r, ok := obj.kmap[k-obj.n][i]; ok {
    return r
  } else {
    return obj.addKmer(c, i)
  }
}

func (obj *KmersCounter) MatchingKmers(c []byte) KmerList {
  k := len(c)
  i := obj.matchingKmers(c)
  r := make(KmerList, len(i))
  for j, _ := range i {
    r[j].K    = k
    r[j].I    = i[j]
    r[j].Name = obj.KmersSet.IdToName(k, i[j])
  }
  return r
}

/* -------------------------------------------------------------------------- */

func (obj *KmersCounter) countKmers(sequence []byte, k int) KmerCounts {
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
  counts := make(map[Kmer]int)
  i      := 0
  for id, c := range r {
    kmers[i] = Kmer{K: k, I: id, Name: obj.KmersSet.IdToName(k, id) }
    counts[kmers[i]] = c
    i++
  }
  kmers.Sort()
  return KmerCounts{Kmers: kmers, Counts: counts}
}

func (obj *KmersCounter) CountKmers(sequence []byte) KmerCounts {
  kmers  := KmerList{}
  counts := make(map[Kmer]int)
  for k := obj.n; k <= obj.m; k++ {
    kmerCounts := obj.countKmers(sequence, k)
    kmers = append(kmers, kmerCounts.Kmers...)
    for kmer, c := range kmerCounts.Counts {
      counts[kmer] = c
    }
  }
  return KmerCounts{Kmers: kmers, Counts: counts}
}

func (obj *KmersCounter) identifyKmers(sequence []byte, k int) KmerCounts {
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
  counts := make(map[Kmer]int)
  i      := 0
  for id, _ := range r {
    kmers[i] = Kmer{K: k, I: id, Name: obj.KmersSet.IdToName(k, id) }
    counts[kmers[i]] = 1
    i++
  }
  kmers.Sort()
  return KmerCounts{Kmers: kmers, Counts: counts}
}

func (obj *KmersCounter) IdentifyKmers(sequence []byte) KmerCounts {
  kmers  := KmerList{}
  counts := make(map[Kmer]int)
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

func (obj *KmersCounter) MinKmerSize() int {
  return obj.n
}

func (obj *KmersCounter) MaxKmerSize() int {
  return obj.m
}

func (obj *KmersCounter) MaxAmbiguous() []int {
  return obj.ma
}

func (obj *KmersCounter) Complement() bool {
  return obj.complement
}

func (obj *KmersCounter) Reverse() bool {
  return obj.reverse
}

func (obj *KmersCounter) Revcomp() bool {
  return obj.revcomp
}

func (obj *KmersCounter) Alphabet() ComplementableAlphabet {
  return obj.al
}
