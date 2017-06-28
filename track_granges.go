/* Copyright (C) 2016, 2017 Philipp Benner
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

/* -------------------------------------------------------------------------- */

func (track GenericTrack) GRanges(name string) (GRanges, error) {
  binSize  := track.GetBinSize()
  seqnames := []string{}
  from     := []int{}
  to       := []int{}
  values   := []float64{}
  for _, name := range track.GetSeqNames() {
    sequence, err := track.GetSequence(name); if err != nil {
      return GRanges{}, err
    }
    if sequence.NBins() == 0 {
      continue
    }
    // current values
    c_from := 0
    c_to   := binSize
    c_val  := sequence.AtBin(0)
    for i := 1; i < sequence.NBins(); i++ {
      if v := sequence.AtBin(i); v != c_val {
        seqnames = append(seqnames, name)
        from     = append(from,   c_from)
        to       = append(to,     c_to)
        values   = append(values, c_val)
        c_from   = c_to
        c_to     = c_from + binSize
        c_val    = v
      } else {
        c_to    += binSize
      }
    }
    // append last result
    seqnames = append(seqnames, name)
    from     = append(from,   c_from)
    to       = append(to,     c_to)
    values   = append(values, c_val)
  }
  r := NewGRanges(seqnames, from, to, nil)
  r.AddMeta(name, values)
  return r, nil
}
