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

type LazyTrack struct {
  Name       string
  BinSize    int
  BinOverlap int
  Bwr        BigWigReader
  Filename   string
  BinSumStat BinSummaryStatistics
  Init       float64
}

/* constructor
 * -------------------------------------------------------------------------- */

func NewLazyTrack(filename, name string, f BinSummaryStatistics, binSize, binOverlap int, init float64) (LazyTrack, error) {
  bwr, err := NewBigWigReader(filename); if err != nil {
    return LazyTrack{}, err
  }
  return LazyTrack{name, binSize, binOverlap, *bwr, filename, f, init}, nil
}

/* -------------------------------------------------------------------------- */

func (track LazyTrack) Clone() LazyTrack {
  track, err := NewLazyTrack(track.Filename, track.Name, track.BinSumStat, track.BinSize, track.BinOverlap, track.Init); if err != nil {
    panic(err)
  }
  return track
}

func (track LazyTrack) CloneTrack() Track {
  return track.Clone()
}

/* access methods
 * -------------------------------------------------------------------------- */

func (track LazyTrack) GetBinsize() int {
  return track.BinSize
}

func (track LazyTrack) GetName() string {
  return track.Name
}

func (track LazyTrack) GetSeqNames() []string {
  return track.Bwr.Genome.Seqnames
}

func (track LazyTrack) GetGenome() Genome {
  return track.Bwr.Genome
}

func (track LazyTrack) GetSequence(query string) (TrackSequence, error) {
  if seq, _, err := track.Bwr.QuerySequence(query, track.BinSumStat, track.BinSize, track.BinOverlap, track.Init); err != nil {
    return TrackSequence{}, err
  } else {
    return TrackSequence{seq, track.BinSize}, nil
  }
}

func (track LazyTrack) GetSlice(r GRangesRow) ([]float64, error) {
  binFrom := r.Range.From/track.BinSize
  binTo   := r.Range.To  /track.BinSize
  seq     := make([]float64, binTo-binFrom)
  for r := range track.Bwr.Query(r.Seqname, r.Range.From, r.Range.To, track.BinSize) {
    if r.Error != nil {
      return nil, r.Error
    }
    for i := r.From; i < r.To; i += track.BinSize {
      seq[i/track.BinSize] = r.Sum/r.Valid
    }
  }
  return seq, nil
}

/* -------------------------------------------------------------------------- */

func (track *LazyTrack) ReadBigWig(filename, name string, f BinSummaryStatistics, binSize, binOverlap int, init float64) error {
  if tmp, err := NewLazyTrack(filename, name, f, binSize, binOverlap, init); err != nil {
    return err
  } else {
    *track = tmp
  }
  return nil
}
