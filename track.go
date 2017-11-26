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

type TrackSequence struct {
  sequence []float64
  binSize    int
}

func (obj TrackSequence) At(i int) float64 {
  return obj.sequence[i/obj.binSize]
}

func (obj TrackSequence) AtBin(i int) float64 {
  return obj.sequence[i]
}

func (obj TrackSequence) NBins() int {
  return len(obj.sequence)
}

func (obj TrackSequence) GetBinSize() int {
  return obj.binSize
}

/* -------------------------------------------------------------------------- */

type TrackMutableSequence struct {
  TrackSequence
}

func (obj TrackSequence) Set(i int, v float64) {
  obj.sequence[i/obj.binSize] = v
}

func (obj TrackSequence) SetBin(i int, v float64) {
  obj.sequence[i] = v
}

/* -------------------------------------------------------------------------- */

type Track interface {
  GetName     ()               string
  GetBinSize  ()               int
  GetSequence (seqname string) (TrackSequence, error)
  GetGenome   ()               Genome
  GetSeqNames ()               []string
  GetSlice    (r GRangesRow)   ([]float64, error)
  CloneTrack  ()               Track
}

type MutableTrack interface {
  Track
  GetMutableSequence(seqname string) (TrackMutableSequence, error)
  CloneMutableTrack ()               MutableTrack
}

/* -------------------------------------------------------------------------- */

type GenericTrack struct {
  Track
}

type GenericMutableTrack struct {
  MutableTrack
}
