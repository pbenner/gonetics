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
import "math"
import "io"
import "os"

/* -------------------------------------------------------------------------- */

type LazyTrack struct {
  Name       string
  BinSize    int
  BinOverlap int
  Bwr        BigWigReader
  BinSumStat BinSummaryStatistics
  Init       float64
}

/* constructor
 * -------------------------------------------------------------------------- */

func NewLazyTrack(reader io.ReadSeeker, name string, f BinSummaryStatistics, binSize, binOverlap int, init float64) (LazyTrack, error) {
  bwr, err := NewBigWigReader(reader); if err != nil {
    return LazyTrack{}, err
  }
  if binSize == 0 {
    for record := range bwr.Query(".*", 0, math.MaxInt64, binSize) {
      if record.Error != nil {
        return LazyTrack{}, record.Error
      }
      binSize = record.To - record.From
      record.Quit()
    }
  }
  return LazyTrack{name, binSize, binOverlap, *bwr, f, init}, nil
}

/* -------------------------------------------------------------------------- */

func (track LazyTrack) Clone() LazyTrack {
  r := track
  return r
}

func (track LazyTrack) CloneTrack() Track {
  return track.Clone()
}

/* access methods
 * -------------------------------------------------------------------------- */

func (track LazyTrack) GetBinSize() int {
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
  if seq, binSize, err := track.Bwr.QuerySequence(query, track.BinSumStat, track.BinSize, track.BinOverlap, track.Init); err != nil {
    return TrackSequence{}, err
  } else {
    return TrackSequence{seq, binSize}, nil
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

func (track *LazyTrack) FilterGenome(f func(name string, length int) bool) {
  track.Bwr.Genome = track.Bwr.Genome.Filter(f)
}

/* -------------------------------------------------------------------------- */

func (track *LazyTrack) ReadBigWig(reader io.ReadSeeker, name string, f BinSummaryStatistics, binSize, binOverlap int, init float64) error {
  if tmp, err := NewLazyTrack(reader, name, f, binSize, binOverlap, init); err != nil {
    return err
  } else {
    *track = tmp
  }
  return nil
}

/* -------------------------------------------------------------------------- */

type LazyTrackFile struct {
  LazyTrack
  f *os.File
}

func (obj *LazyTrackFile) ImportBigWig(filename, name string, s BinSummaryStatistics, binSize, binOverlap int, init float64) error {
  f, err := os.Open(filename)
  if err != nil {
    return err
  }

  if tmp, err := NewLazyTrack(f, name, s, binSize, binOverlap, init); err != nil {
    return err
  } else {
    obj.LazyTrack = tmp
  }
  return nil
}

func (obj *LazyTrackFile) Close() error {
  return obj.f.Close()
}
