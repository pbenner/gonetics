/* Copyright (C) 2016 Philipp Benner
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

/* endPointList
 * -------------------------------------------------------------------------- */

type findNearestHits struct {
    queryHits []int
  subjectHits []int
  distances   []int
}

func (obj findNearestHits) Len() int {
  return len(obj.queryHits)
}

func (obj findNearestHits) Less(i, j int) bool {
  if obj.queryHits[i] == obj.queryHits[j] {
    if obj.distances[i] == obj.distances[j] {
      return obj.subjectHits[i] < obj.subjectHits[j]
    } else {
      return obj.distances[i] < obj.distances[j]
    }
  } else {
    return obj.queryHits[i] < obj.queryHits[j]
  }
}

func (obj findNearestHits) Swap(i, j int) {
  obj.  queryHits[i], obj.  queryHits[j] = obj.  queryHits[j], obj.  queryHits[i]
  obj.subjectHits[i], obj.subjectHits[j] = obj.subjectHits[j], obj.subjectHits[i]
  obj.  distances[i], obj.  distances[j] = obj.  distances[j], obj.  distances[i]
}

/* -------------------------------------------------------------------------- */

func (r1 *endPoint) distance(r2 *endPoint) (int, int) {
  sign := -1
  if r1.getStart() > r2.getStart() {
    r1, r2 = r2, r1
    sign   = 1
  }
  if r1.getEnd() >= r2.getStart() {
    return 0, sign
  }
  d1 := r2.getStart() - r1.getEnd()
  d2 := r2.getEnd()   - r1.getEnd()
  if d1 < d2 {
    return d1, sign
  } else {
    return d2, sign
  }
}

/* FindNearest
 * -------------------------------------------------------------------------- */

// For every query region find the k nearest subject regions including
// all overlapping regions.
func FindNearest(query, subject GRanges, k int) ([]int, []int, []int) {

  n :=   query.Length()
  m := subject.Length()

    queryHits := []int{}
  subjectHits := []int{}
    distances := []int{}

  rmap := make(map[string]endPointList)
  // fill map
  for i := 0; i < n; i++ {
    start := endPoint{query.Ranges[i].From,  nil, nil, i, true}
    end   := endPoint{query.Ranges[i].To, &start, nil, i, true}
    start.end = &end

    entry := rmap[query.Seqnames[i]]
    entry  = append(entry, start)
    entry  = append(entry, end)
    rmap[query.Seqnames[i]] = entry
  }
  for i := 0; i < m; i++ {
    start := endPoint{subject.Ranges[i].From,  nil, nil, i, false}
    end   := endPoint{subject.Ranges[i].To, &start, nil, i, false}
    start.end = &end

    entry := rmap[subject.Seqnames[i]]
    entry  = append(entry, start)
    entry  = append(entry, end)
    rmap[subject.Seqnames[i]] = entry
  }
  // sort map entries
  for _, entry := range rmap {
    sort.Sort(entry)
  }
  // find k closest subjects
  for _, entry := range rmap {
    q, s := findOverlapsEntry(nil, nil, entry)
    for j := 0; j < len(q); j++ {
        queryHits = append(  queryHits, q[j])
      subjectHits = append(subjectHits, s[j])
        distances = append(  distances, 0)
    }
    n := len(entry)
    for i, r := range entry {
      // loop over query start regions
      if r.isQuery && r.end != nil {
        i1 := i-1 // looping to the left
        i2 := i+1 // looping to the right
        // find k nearest neighbors
        for j := 0; j < k && (i1 >= 0 || i2 < len(entry)); j++ {
          ir := -1
          dr := -1
          // find next subject end to the left
          for ; i1 >= 0; i1-- {
            if !entry[i1].isQuery && entry[i1].isEnd() {
              break
            }
          }
          // find next subject start to the right (and drop overlaps)
          for ; i2 < n; i2++ {
            if !entry[i2].isQuery && entry[i2].isStart() && entry[i2].position > r.getEnd() {
              break
            }
          }
          // are there two elements to compare?
          if i1 >= 0 && i2 < n {
            d1, s1 := r.distance(&entry[i1])
            d2, s2 := r.distance(&entry[i2])
            if d1 <= d2 {
              dr = d1*s1; ir = i1; i1--
            } else {
              dr = d2*s2; ir = i2; i2++
            }
          } else {
            if i1 >= 0 {
              d1, s1 := r.distance(&entry[i1])
              dr = d1*s1; ir = i1; i1--
            }
            if i2 <  n {
              d2, s2 := r.distance(&entry[i2])
              dr = d2*s2; ir = i2; i2++
            }
          }
          if ir != -1 {
              queryHits = append(  queryHits, entry[i ].srcIdx)
            subjectHits = append(subjectHits, entry[ir].srcIdx)
              distances = append(  distances, dr)
          }
        }
      }
    }
  }
  sort.Sort(findNearestHits{queryHits, subjectHits, distances})

  return queryHits, subjectHits, distances
}
