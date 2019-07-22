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

import "fmt"
import "sort"

/* K-mer equivalence class
 * -------------------------------------------------------------------------- */

type KmerClassId struct {
  K int
  I int
}

/* -------------------------------------------------------------------------- */

type KmerClass struct {
  KmerClassId
  // list of equivalent K-mers
  Elements []string
}

/* -------------------------------------------------------------------------- */

func NewKmerClass(k, i int, elements []string) KmerClass {
  return KmerClass{KmerClassId{K: k, I: i}, elements}
}

/* -------------------------------------------------------------------------- */

func (obj KmerClass) String() string {
  if len(obj.Elements) == 0 {
    return ""
  }
  str := obj.Elements[0]
  for i := 1; i < len(obj.Elements); i++ {
    str = fmt.Sprintf("%s|%s", str, obj.Elements[i])
  }
  return str
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

func (obj KmerClass) Matches(b KmerClass, alphabet ComplementableAlphabet) bool {
  for _, elem1 := range obj.Elements {
    for _, elem2 := range b.Elements {
      if Kmer(elem1).Matches(Kmer(elem2), alphabet) {
        return true
      }
    }
  }
  return false
}

func (obj KmerClass) CountAmbiguous(alphabet ComplementableAlphabet) int {
  m := 0
  s := []byte(obj.Elements[0])
  for i := 0; i < len(s); i++ {
    if ok, err := alphabet.IsAmbiguous(s[i]); err != nil {
      panic(err)
    } else {
      if ok {
        m++
      }
    }
  }
  return m
}

func (obj KmerClass) CountWildcard(alphabet ComplementableAlphabet) int {
  m := 0
  s := []byte(obj.Elements[0])
  for i := 0; i < len(s); i++ {
    if ok, err := alphabet.IsWildcard(s[i]); err != nil {
      panic(err)
    } else {
      if ok {
        m++
      }
    }
  }
  return m
}

/* -------------------------------------------------------------------------- */

type KmerClassList []KmerClass

func (obj KmerClassList) Clone() KmerClassList {
  r := make(KmerClassList, len(obj))
  copy(r, obj)
  return r
}

func (obj KmerClassList) Equals(b KmerClassList) bool {
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

func (obj KmerClassList) Len() int {
  return len(obj)
}

func (obj KmerClassList) Less(i, j int) bool {
  if obj[i].K != obj[j].K {
    return obj[i].K < obj[j].K
  } else {
    return obj[i].I < obj[j].I
  }
}

func (obj KmerClassList) Swap(i, j int) {
  obj[i], obj[j] = obj[j], obj[i]
}

func (obj KmerClassList) Sort() {
  sort.Sort(obj)
}

func (obj KmerClassList) Union(b ...KmerClassList) KmerClassList {
  m := make(KmerClassSet)
  for _, class := range obj {
    m[class.KmerClassId] = class.Elements
  }
  for _, bi := range b {
    for _, class := range bi {
      m[class.KmerClassId] = class.Elements
    }
  }
  return m.AsList()
}

/* -------------------------------------------------------------------------- */

type KmerClassSet map[KmerClassId][]string

func (obj KmerClassSet) AsList() KmerClassList {
  r := make(KmerClassList, len(obj))
  i := 0
  for id, elements := range obj {
    r[i].K = id.K
    r[i].I = id.I
    r[i].Elements = elements
    i++
  }
  r.Sort()
  return r
}
