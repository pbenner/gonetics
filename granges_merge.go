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

func (obj GRanges) merge(seqname string, entry endPointList) GRanges {
  seqnames := []string{}
  from     := []int{}
  to       := []int{}
  for i := 0; i < len(entry); i++ {
    // next item must be a start position
    if entry[i].start != nil {
      panic("internal error")
    }
    r_from := entry[i].position
    r_to   := entry[i].position+1
    // k: number open intervals
    for k := 1; k > 0; {
      // go to next item
      i += 1
      if entry[i].start == nil {
        // start of an interval
        k   += 1
        r_to = entry[i].position
      } else {
        // end of an interval
        k   -= 1
        r_to = entry[i].position+1
      }
    }
    seqnames = append(seqnames, seqname)
    from     = append(from,     r_from)
    to       = append(to,       r_to)
  }
  return obj.Append(NewGRanges(seqnames, from, to, nil))
}

func (obj GRanges) Merge(granges ...GRanges) GRanges {
  r    := GRanges{}
  rmap := make(map[string]endPointList)
  // fill map
  for _, g := range append(granges, obj) {
    for i := 0; i < g.Length(); i++ {
      start := endPoint{g.Ranges[i].From,  nil, nil, i, true}
      end   := endPoint{g.Ranges[i].To-1, &start, nil, i, true}
      entry := rmap[g.Seqnames[i]]
      entry  = append(entry, start)
      entry  = append(entry, end)
      rmap[g.Seqnames[i]] = entry
    }
  }
  seqnames := []string{}
  for key, _ := range rmap {
    seqnames = append(seqnames, key)
  }
  sort.Strings(seqnames)
  // sort map entries
  for _, entry := range rmap {
    sort.Sort(entry)
  }
  // find overlaps
  for _, seqname := range seqnames {
    r = r.merge(seqname, rmap[seqname])
  }
  return r
}
