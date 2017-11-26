/* Copyright (C) 2016-2017 Philipp Benner
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


package main

/* -------------------------------------------------------------------------- */

import   "fmt"
import   "log"
import   "math"
import   "strconv"
import   "strings"
import   "os"

import   "github.com/pborman/getopt"

import . "github.com/pbenner/gonetics"

/* -------------------------------------------------------------------------- */

type Config struct {
  Meta    []string
  Exclude    string
  Regions    string
  KNearest   int
  BinSize    int
  BinOver    int
  BinStat    BinSummaryStatistics
  TrackInit  float64
  Verbose    int
}

/* -------------------------------------------------------------------------- */

func PrintStderr(config Config, level int, format string, args ...interface{}) {
  if config.Verbose >= level {
    fmt.Fprintf(os.Stderr, format, args...)
  }
}

/* -------------------------------------------------------------------------- */

func importBed3(config Config, filename string) GRanges {
  granges := GRanges{}
  PrintStderr(config, 1, "Reading bed file `%s'... ", filename)
  if err := granges.ImportBed3(filename); err != nil {
    PrintStderr(config, 1, "failed\n")
    log.Fatal(err)
  } else {
    PrintStderr(config, 1, "done\n")
  }
  return granges
}

func importBed6(config Config, filename string) GRanges {
  granges := GRanges{}
  PrintStderr(config, 1, "Reading bed file `%s'... ", filename)
  if err := granges.ImportBed6(filename); err != nil {
    PrintStderr(config, 1, "failed\n")
    log.Fatal(err)
  } else {
    PrintStderr(config, 1, "done\n")
  }
  return granges
}

func importTable(config Config, filename string, names, types []string) GRanges {
  granges := GRanges{}
  PrintStderr(config, 1, "Reading table `%s'... ", filename)
  if err := granges.ImportTable(filename, names, types); err != nil {
    PrintStderr(config, 1, "failed\n")
    log.Fatal(err)
  } else {
    PrintStderr(config, 1, "done\n")
  }
  return granges
}

func exportTable(config Config, granges GRanges, filename string, header, strand, compress bool, args ...interface{}) {
  PrintStderr(config, 1, "Writing table `%s'... ", filename)
  if err := granges.ExportTable(filename, header, strand, compress, args...); err != nil {
    PrintStderr(config, 1, "failed\n")
    log.Fatal(err)
  } else {
    PrintStderr(config, 1, "done\n")
  }
}

func importLazyTrack(config Config, trackFilename string) (LazyTrackFile, error) {
  track := LazyTrackFile{}
  PrintStderr(config, 1, "Lazy importing track `%s'... ", trackFilename)
  if err := track.ImportBigWig(trackFilename, "", config.BinStat, config.BinSize, config.BinOver, config.TrackInit); err != nil {
    PrintStderr(config, 1, "failed\n")
    return track, err
  }
  PrintStderr(config, 1, "done\n")
  return track, nil
}

/* -------------------------------------------------------------------------- */

func getBinSummaryStatistics(str string) BinSummaryStatistics {
  switch str {
  case "mean":
    return BinMean
  case "discrete mean":
    return BinDiscreteMean
  case "min":
    return BinMin
  case "max":
    return BinMax
  }
  log.Fatal("invalid bin summary statistics: %s", str)
  return nil
}

/* -------------------------------------------------------------------------- */

func findNearestRegions(config Config, r GRanges) GRanges {
  if config.Regions == "" {
    return r
  }
  regions := GRanges{}

  if strings.HasSuffix(config.Regions, ".bed") {
    regions = importBed6(config, config.Regions)
  } else {
    meta_tmp := append(config.Meta, "name")
    types := []string{}
    for i := 0; i < len(meta_tmp); i++ {
      types = append(types, "[]string")
    }
    regions = importTable(config, config.Regions, meta_tmp, types)
  }

  subjectNames := regions.GetMetaStr("name")
  if len(subjectNames) == 0 {
    log.Fatal("regions file `%s' has no name column")
  }
  queryHits, subjectHits, distances := FindNearest(r, regions, config.KNearest)

  regionsNames     := make([][]string, r.Length())
  regionsDistances := make([][]int, r.Length())

  // add regions names and distances
  for i := 0; i < len(queryHits); i++ {
    qi :=   queryHits[i]
    si := subjectHits[i]
    regionsNames    [qi] = append(regionsNames    [qi], subjectNames[si])
    regionsDistances[qi] = append(regionsDistances[qi], distances[i])
  }
  r.AddMeta("names",     regionsNames)
  r.AddMeta("distances", regionsDistances)

  // add meta columns
  for j := 0; j < len(config.Meta); j++ {
    col := make([][]string, r.Length())
    met := regions.GetMetaStr(config.Meta[j])
    if len(met) == 0 {
      log.Fatalf("regions file `%s' has no column named `%s'", config.Regions, config.Meta[j])
    }
    for i := 0; i < len(queryHits); i++ {
      qi :=   queryHits[i]
      si := subjectHits[i]
      col[qi] = append(col[qi], met[si])
    }
    r.AddMeta(config.Meta[j], col)
  }

  return r
}

/* -------------------------------------------------------------------------- */

func allPositive(sequences []TrackSequence, thresholds []float64, i int) bool {
  for j := 0; j < len(sequences); j++ {
    if math.IsNaN(sequences[j].AtBin(i)) || sequences[j].AtBin(i) <= thresholds[j] {
      return false
    }
  }
  return true
}

func allSum(sequences []TrackSequence, i int) float64 {
  sum := 0.0
  for j := 0; j < len(sequences); j++ {
    sum += sequences[j].AtBin(i)
  }
  return sum
}

func getJointPeaks(tracks []Track, thresholds []float64) (GRanges, error) {
  if len(tracks) != len(thresholds) {
    return GRanges{}, fmt.Errorf("GetJointPeaks(): invalid arguments")
  }
  if len(tracks) == 0 {
    return GRanges{}, nil
  }
  seqnames := []string{}
  from     := []int{}
  to       := []int{}
  strand   := []byte{}
  test     := [][]float64{}

  for _, name := range tracks[0].GetSeqNames() {
    s, err := tracks[0].GetSequence(name); if err != nil {
      return GRanges{}, err
    }
    sequences := []TrackSequence{s}
    binsize   := s.GetBinSize()
    seqlen    := s.NBins()
    // check remaining track for consistency
    for j := 1; j < len(tracks); j++ {
      if sequence, err := tracks[j].GetSequence(name); err != nil {
        return GRanges{}, fmt.Errorf("reading sequence from track `%d' failed: %v", j+1, err)
      } else {
        if sequence.GetBinSize() != binsize {
          return GRanges{}, fmt.Errorf("tracks `1' and `%d' have different bin sizes (`%d' and `%d')", j+1, binsize, sequence.GetBinSize())
        }
        if sequence.NBins() != seqlen {
          return GRanges{}, fmt.Errorf("sequence `%s' on track `1' and `%d' have different lengths (`%d' and `%d')", name, j+1, seqlen, sequence.NBins())
        }
        sequences = append(sequences, sequence)
      }
    }
    for i := 0; i < seqlen; i++ {
      if allPositive(sequences, thresholds, i) {
        // peak begins here
        i_from := i
        // maximum value
        v_max  := allSum(sequences, i)
        // position of the maximum value
        i_max  := i
        // increment until either the sequence ended or
        // the value drops below the threshold
        for i < seqlen && allPositive(sequences, thresholds, i) {
          if sum := allSum(sequences, i); sum > v_max {
            // update maximum position and value
            i_max = i
            v_max = sum
          }
          i += 1
        }
        tmp := make([]float64, len(sequences))
        for j := 0; j < len(sequences); j++ {
          tmp[j] = sequences[j].AtBin(i_max)
        }
        // save peak
        seqnames = append(seqnames, name)
        from     = append(from, i_from*binsize)
        to       = append(to,   i     *binsize)
        test     = append(test, tmp)
      }
    }
  }
  peaks := NewGRanges(seqnames, from, to, strand)
  peaks.AddMeta("test", test)
  // sum up test results for sorting rows
  peaks.ReduceFloat("test","test.sum", func(x []float64) float64 {
    sum := 0.0
    for i := 0; i < len(x); i++ {
      sum += x[i]
    }
    return sum
  })
  peaks, _ = peaks.Sort("test.sum", true)
  peaks.DeleteMeta("test.sum")

  return peaks, nil
}

/* -------------------------------------------------------------------------- */

func positive(config Config, filenameOut string, filenameIn []string, thresholds []float64) {

  tracks  := []Track{}
  exclude := GRanges{}
  if config.Exclude != "" {
    exclude = importBed3(config, config.Exclude)
  }

  for i := 0; i < len(filenameIn); i++ {
    track, err := importLazyTrack(config, filenameIn[i]); if err != nil {
      log.Fatal(err)
    }
    tracks = append(tracks, track)
  }

  if granges, err := getJointPeaks(tracks, thresholds); err != nil {
    log.Fatal(err)
  } else {
    granges.RemoveOverlapsWith(exclude)
    granges = findNearestRegions(config, granges)
    exportTable(config, granges, filenameOut, true, false, false, OptionPrintScientific{true})
  }
}

/* -------------------------------------------------------------------------- */

func main() {

  config    := Config{}
  threshold := 0.0

  options := getopt.New()

  optBinSize   := options.    IntLong("bin-size",       0 ,      0, "bin size")
  optBinStat   := options. StringLong("bin-summary",    0 , "mean", "bin summary statistic [mean (default), max, min, discrete mean]")
  optBinOver   := options.    IntLong("bin-overlap",    0 ,      0, "number of overlapping bins when computing the summary")
  optExclude   := options. StringLong("exclude",        0 ,     "", "bed file with excluded regions")
  optRegions   := options. StringLong("regions",        0 ,     "", "augment positive list with names of k nearest regions [format: bed6 with name column or table]")
  optKNearest  := options.    IntLong("k-nearest",      0 ,      1, "number of nearest regions")
  optMeta      := options. StringLong("regions-meta",   0 ,     "", "comma separated list of column names added to the resulting table")
  optThreshold := options. StringLong("threshold",      0 ,  "1.0", "default threshold value for all tracks")
  optVerbose   := options.CounterLong("verbose",       'v',         "verbose level [-v or -vv]")
  optHelp      := options.   BoolLong("help",          'h',         "print help")

  options.SetParameters("<RESULT.table> <INPUT.bw[:THRESHOLD]>...")
  options.Parse(os.Args)

  if *optHelp {
    options.PrintUsage(os.Stdout)
    os.Exit(0)
  }
  if len(options.Args()) < 2 {
    options.PrintUsage(os.Stderr)
    os.Exit(1)
  }
  if *optMeta != "" {
    config.Meta = strings.Split(*optMeta, ",")
  }
  if t, err := strconv.ParseFloat(*optThreshold, 64); err != nil {
    log.Fatal(err)
  } else {
    threshold = t
  }
  config.Verbose  = *optVerbose
  config.Exclude  = *optExclude
  config.Regions  = *optRegions
  config.KNearest = *optKNearest
  config.BinSize  = *optBinSize
  config.BinOver  = *optBinOver
  config.BinStat  = getBinSummaryStatistics(*optBinStat)

  filenameOut := options.Args()[0]
  filenameIn  := [] string{}
  thresholds  := []float64{}

  // parse input files
  for i := 1; i < len(options.Args()); i++ {
    tmp := strings.Split(options.Args()[i], ":")
    if len(tmp) == 1 {
      filenameIn = append(filenameIn, tmp[0])
      thresholds = append(thresholds, threshold)
    } else
    if len(tmp) == 2 {
      t, err := strconv.ParseFloat(tmp[1], 64)
      if err != nil {
        log.Fatal(err)
      }
      filenameIn = append(filenameIn, tmp[0])
      thresholds = append(thresholds, t)
    } else {
      log.Fatalf("invalid argument `%s'", options.Args()[i])
    }
  }
  positive(config, filenameOut, filenameIn, thresholds)
}
