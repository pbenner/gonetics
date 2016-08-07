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
import "bufio"
import "bytes"

/* -------------------------------------------------------------------------- */

type GRangesRow struct {
  Seqname   string
  Range     Range
  Strand    byte
  MetaRow
}

/* -------------------------------------------------------------------------- */

func NewGRangesRow(granges GRanges, i int) GRangesRow {
  return GRangesRow{
    granges.Seqnames[i],
    granges.Ranges  [i],
    granges.Strand  [i],
    granges.Meta.Row(i) }
}

/* -------------------------------------------------------------------------- */

func (r GRangesRow) String() string {
  var buffer bytes.Buffer
  writer := bufio.NewWriter(&buffer)

  if r.Strand == '*' {
    fmt.Fprintf(writer, "%s:%d-%d", r.Seqname, r.Range.From, r.Range.To)
  } else {
    fmt.Fprintf(writer, "%s:%d-%d (%c)", r.Seqname, r.Range.From, r.Range.To, r.Strand)
  }
  writer.Flush()

  return buffer.String()
}
