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

import "fmt"
import "sort"

/* -------------------------------------------------------------------------- */

type endPoint struct {
  position  int
  start    *endPoint
  end      *endPoint
  srcIdx    int
  isQuery   bool
}

func (obj *endPoint) isStart() bool {
  return obj.start == nil
}

func (obj *endPoint) isEnd() bool {
  return obj.end == nil
}

func (obj *endPoint) getStart() int {
  if obj.start == nil {
    return obj.position
  } else {
    return obj.start.position
  }
}

func (obj *endPoint) getEnd() int {
  if obj.end == nil {
    return obj.position
  } else {
    return obj.end.position
  }
}

func (obj *endPoint) String() string {
  return fmt.Sprintf("<%d,%d>", obj.getStart(), obj.getEnd())
}

/* endPointList
 * -------------------------------------------------------------------------- */

type endPointList []endPoint

func NewEndPointList() endPointList {
  r := []endPoint{}
  return r
}

func (r endPointList) Len() int {
  return len(r)
}

func (r endPointList) Less(i, j int) bool {
  if r[i].position != r[j].position {
    return r[i].position < r[j].position
  }
  // if the positions equal, give priority to the beginning of a range
  if r[i].start == nil && r[j].start != nil {
    return true
  } else {
    return false
  }
}

func (r endPointList) Swap(i, j int) {
  r[i], r[j] = r[j], r[i]
}

func (s *endPointList) Append(r endPoint) {
  (*s) = append(*s, r)
}

func (s *endPointList) Remove(r endPoint) {
  for i := 0; i < len(*s); i++ {
    if (*s)[i] == r {
      (*s) = append((*s)[0:i], (*s)[i+1:len(*s)]...)
    }
  }
}

/* FindOverlaps
 * -------------------------------------------------------------------------- */

func findOverlapsEntry(queryHits, subjectHits []int, entry endPointList) ([]int, []int) {
    queryList := NewEndPointList()
  subjectList := NewEndPointList()
  for _, r := range entry {
    if r.isQuery {
      if r.start == nil {
        queryList.Append(r)
        // record overlaps (all elements in subjectList overlap with this
        // position)
        for i := 0; i < len(subjectList); i++ {
            queryHits = append(  queryHits, r.srcIdx)
          subjectHits = append(subjectHits, subjectList[i].srcIdx)
        }
      } else {
        // remove start position from stack
        queryList.Remove(*r.start)
      }
    } else {
      if r.start == nil {
        subjectList.Append(r)
        // record overlaps (all elements in queryList overlap with this
        // position)
        for i := 0; i < len(queryList); i++ {
            queryHits = append(  queryHits, queryList[i].srcIdx)
          subjectHits = append(subjectHits, r.srcIdx)
        }
      } else {
        // remove start position from stack
        subjectList.Remove(*r.start)
      }
    }
  }
  return queryHits, subjectHits
}

func FindOverlaps(query, subject GRanges) ([]int, []int) {

  n :=   query.Length()
  m := subject.Length()

    queryHits := []int{}
  subjectHits := []int{}

  rmap := make(map[string]endPointList)
  // fill map
  for i := 0; i < n; i++ {
    start := endPoint{query.Ranges[i].From,  nil, nil, i, true}
    end   := endPoint{query.Ranges[i].To, &start, nil, i, true}
    entry := rmap[query.Seqnames[i]]
    entry  = append(entry, start)
    entry  = append(entry, end)
    rmap[query.Seqnames[i]] = entry
  }
  for i := 0; i < m; i++ {
    start := endPoint{subject.Ranges[i].From,  nil, nil, i, false}
    end   := endPoint{subject.Ranges[i].To, &start, nil, i, false}
    entry := rmap[subject.Seqnames[i]]
    entry  = append(entry, start)
    entry  = append(entry, end)
    rmap[subject.Seqnames[i]] = entry
  }
  // sort map entries
  for _, entry := range rmap {
    sort.Sort(entry)
  }
  // find overlaps
  for _, entry := range rmap {
    queryHits, subjectHits = findOverlapsEntry(queryHits, subjectHits, entry)
  }
  return queryHits, subjectHits
}
