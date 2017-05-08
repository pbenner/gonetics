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

/* -------------------------------------------------------------------------- */

type Track interface {
  GetName     ()               string
  GetBinsize  ()               int
  GetSequence (seqname string) ([]float64, error)
  GetGenome   ()               Genome
  GetSeqNames ()               []string
  GetSlice    (r GRangesRow)   ([]float64, error)
  CloneTrack  ()               Track
  At (seqname string, position int) (float64, error)
  Set(seqname string, position int, value float64) error
  ReadBigWig  (filename, name string, f BinSummaryStatistics, binsize int) error
  WriteBigWig (filename, description string, genome Genome, args... interface{}) error
}
