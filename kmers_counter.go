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

/* -------------------------------------------------------------------------- */

func (obj *KmersCounter) countKmers(sequence []byte, k int) ([]int, []int) {
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
  ids    := []int{}
  counts := []int{}
  for id, c := range r {
    ids    = append(ids   , id)
    counts = append(counts, c)
  }
  sortIntPairs{ids, counts}.Sort()
  return ids, counts
}

func (obj *KmersCounter) CountKmers(sequence []byte) ([]string, []int) {
  names  := []string{}
  counts := []int{}
  for k := obj.n; k <= obj.m; k++ {
    ids, c := obj.countKmers(sequence, k)
    for i, id := range ids {
      names  = append(names , obj.KmersSet.IdToName(k, id))
      counts = append(counts, c[i])
    }
  }
  return names, counts
}

func (obj *KmersCounter) identifyKmers(sequence []byte, k int) []int {
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
  ids := []int{}
  for id, _ := range r {
    ids = append(ids, id)
  }
  sort.Ints(ids)
  return ids
}

func (obj *KmersCounter) IdentifyKmers(sequence []byte) []string {
  names := []string{}
  for k := obj.n; k <= obj.m; k++ {
    ids := obj.identifyKmers(sequence, k)
    for _, id := range ids {
      names = append(names , obj.KmersSet.IdToName(k, id))
    }
  }
  return names
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
