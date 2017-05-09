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

/* -------------------------------------------------------------------------- */

type TMapType map[string][]float64

// A track is a container for experimental data mapped to genomic
// locations. The data is binned in order to reduce memory usage.
// The first position in a sequence is numbered 0.
type SimpleTrack struct {
  Name    string
  Genome  Genome
  Data    TMapType
  Binsize int
}

/* constructor
 * -------------------------------------------------------------------------- */

func NewSimpleTrack(name string, seqnames []string, sequences [][]float64, binsize int) (SimpleTrack, error) {
  if len(seqnames) != len(sequences) {
    return SimpleTrack{}, fmt.Errorf("invalid arguments")
  }
  data    := make(TMapType)
  lengths := make([]int, len(sequences))
  for i, sequence := range sequences {
    lengths[i] = binsize*len(sequence)
    data[seqnames[i]] = sequence
  }
  genome := NewGenome(seqnames, lengths)

  return SimpleTrack{name, genome, data, binsize}, nil
}

func AllocSimpleTrack(name string, genome Genome, binsize int) SimpleTrack {
  data := make(TMapType)

  for i := 0; i < genome.Length(); i++ {
    // by convention drop the last positions if they do not fully
    // cover the last bin (i.e. round down), this is required by
    // wig related tools
    data[genome.Seqnames[i]] = make([]float64,
      divIntDown(genome.Lengths[i], binsize))
  }
  return SimpleTrack{name, genome, data, binsize}
}

func EmptySimpleTrack(name string) SimpleTrack {
  data := make(TMapType)
  return SimpleTrack{name, Genome{}, data, 0}
}

/* -------------------------------------------------------------------------- */

func (track SimpleTrack) Clone() SimpleTrack {
  name    := track.Name
  binsize := track.Binsize
  data    := make(TMapType)
  genome  := track.Genome.Clone()

  for name, sequence := range track.Data {
    t := make([]float64, len(sequence))
    copy(t, sequence)
    data[name] = t
  }
  return SimpleTrack{name, genome, data, binsize}
}

func (track SimpleTrack) CloneTrack() Track {
  return track.Clone()
}

func (track SimpleTrack) CloneMutableTrack() MutableTrack {
  return track.Clone()
}

/* access methods
 * -------------------------------------------------------------------------- */

func (track SimpleTrack) GetBinsize() int {
  return track.Binsize
}

func (track SimpleTrack) Index(position int) int {
  if position < 0 {
    panic("negative position")
  }
  return position/track.Binsize
}

func (track SimpleTrack) GetName() string {
  return track.Name
}

func (track SimpleTrack) GetSeqNames() []string {
  return track.Genome.Seqnames
}

func (track SimpleTrack) GetGenome() Genome {
  return track.Genome
}

func (track SimpleTrack) GetSequence(query string) (TrackSequence, error) {
  if seq, ok := track.Data[query]; !ok {
    return TrackSequence{}, fmt.Errorf("sequence `%s' not found", query)
  } else {
    return TrackSequence{seq, track.Binsize}, nil
  }
}

func (track SimpleTrack) GetMutableSequence(query string) (TrackMutableSequence, error) {
  for name, seq := range track.Data {
    if name == query {
      return TrackMutableSequence{TrackSequence{seq, track.Binsize}}, nil
    }
  }
  return TrackMutableSequence{}, fmt.Errorf("sequence `%s' not found", query)
}

func (track SimpleTrack) GetSlice(r GRangesRow) ([]float64, error) {
  seq, ok := track.Data[r.Seqname]
  if !ok {
    return nil, fmt.Errorf("GetSlice(): invalid seqname `%s'", r.Seqname)
  }
  from := r.Range.From/track.Binsize
  to   := r.Range.To  /track.Binsize
  if from >= len(seq) {
    return nil, nil
  }
  if to < 0 {
    return nil, nil
  }
  if from < 0 {
    from = 0
  }
  if to > len(seq) {
    to = len(seq)
  }
  return seq[from:to], nil
}
