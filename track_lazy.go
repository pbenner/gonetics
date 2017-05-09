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
  Binsize    int
  Bwr        BigWigReader
  Filename   string
  BinSumStat BinSummaryStatistics
}

/* constructor
 * -------------------------------------------------------------------------- */

func NewLazyTrack(filename, name string, binsize int) (LazyTrack, error) {
  bwr, err := NewBigWigReader(filename); if err != nil {
    return LazyTrack{}, err
  }
  return LazyTrack{name, binsize, *bwr, filename, BinMean}, nil
}

/* -------------------------------------------------------------------------- */

func (track LazyTrack) Clone() LazyTrack {
  track, err := NewLazyTrack(track.Filename, track.Name, track.Binsize); if err != nil {
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
  return track.Binsize
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
  if seq, err := track.Bwr.QuerySequence(query, track.BinSumStat, track.Binsize); err != nil {
    return TrackSequence{}, err
  } else {
    return TrackSequence{seq, track.Binsize}, nil
  }
}

func (track LazyTrack) GetSlice(r GRangesRow) ([]float64, error) {
  binFrom := r.Range.From/track.Binsize
  binTo   := r.Range.To  /track.Binsize
  seq     := make([]float64, binTo-binFrom)
  for r := range track.Bwr.Query(r.Seqname, r.Range.From, r.Range.To, track.Binsize) {
    if r.Error != nil {
      return nil, r.Error
    }
    for i := r.From; i < r.To; i += track.Binsize {
      seq[i/track.Binsize] = r.Sum/r.Valid
    }
  }
  return seq, nil
}
