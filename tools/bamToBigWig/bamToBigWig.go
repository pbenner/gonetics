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
import   "path/filepath"
import   "strconv"
import   "strings"
import   "os"

import   "github.com/pborman/getopt"
import . "github.com/pbenner/gonetics"

import   "gonum.org/v1/plot"
import   "gonum.org/v1/plot/plotter"
import   "gonum.org/v1/plot/plotutil"
import   "gonum.org/v1/plot/vg"

/* -------------------------------------------------------------------------- */

type Config struct {
  Verbose                int
  BinSummaryStatistics   string
  BWZoomLevels         []int
  BinningMethod          string
  BinSize                int
  BinOverlap             int
  TrackInit              float64
  PairedEnd              bool
  NormalizeTrack         string
  ShiftReads          [2]int
  LogScale               bool
  Pseudocounts        [2]float64
  Fraglen                int
  FraglenRange        [2]int
  FilterMapQ             int
  FilterReadLengths   [2]int
  FilterDuplicates       bool
  FilterStrand           byte
  SmoothenControl        bool
  SmoothenSizes        []int
  SmoothenMin            float64
  SaveFraglen            bool
  SaveCrossCorr          bool
  SaveCrossCorrPlot      bool
}

/* i/o
 * -------------------------------------------------------------------------- */

func PrintStderr(config Config, level int, format string, args ...interface{}) {
  if config.Verbose >= level {
    fmt.Fprintf(os.Stderr, format, args...)
  }
}

/* utility
 * -------------------------------------------------------------------------- */

func parseFilename(filename string) (string, int) {
  if tmp := strings.Split(filename, ":"); len(tmp) == 2 {
    t, err := strconv.ParseInt(tmp[1], 10, 64)
    if err != nil {
      log.Fatal(err)
    }
    return tmp[0], int(t)
  } else
  if len(tmp) >= 2 {
    log.Fatal("invalid input file description `%s'", filename)
  }
  return filename, 0
}

/* read filters
 * -------------------------------------------------------------------------- */

func filterDuplicates(config Config, reads GRanges) GRanges {
  if config.FilterDuplicates == false {
    return reads
  }
  PrintStderr(config, 1, "Filtering reads (duplicates)... ")
  idx := []int{}
  if config.PairedEnd {
    flag1 := reads.GetMetaInt("flag1")
    flag2 := reads.GetMetaInt("flag2")
    for i := 0; i < reads.Length(); i++ {
      if BamFlag(flag1[i]).Duplicate() || BamFlag(flag2[i]).Duplicate() {
        idx = append(idx, i)
      }
    }
  } else {
    flag := reads.GetMetaInt("flag")
    for i := 0; i < reads.Length(); i++ {
      if BamFlag(flag[i]).Duplicate() {
        idx = append(idx, i)
      }
    }
  }
  PrintStderr(config, 1, "done\n")
  PrintStderr(config, 1, "Filtered out %d reads (%.2f%%)\n", len(idx), 100.0*float64(len(idx))/float64(reads.Length()))
  return reads.Remove(idx)
}

func filterStrand(config Config, reads GRanges) GRanges {
  if config.FilterStrand == '*' {
    return reads
  }
  PrintStderr(config, 1, "Filtering reads (strand: %c)... ", config.FilterStrand)
  idx := []int{}
  for i := 0; i < reads.Length(); i++ {
    if reads.Strand[i] != config.FilterStrand {
      idx = append(idx, i)
    }
  }
  PrintStderr(config, 1, "done\n")
  PrintStderr(config, 1, "Filtered out %d reads (%.2f%%)\n", len(idx), 100.0*float64(len(idx))/float64(reads.Length()))
  return reads.Remove(idx)
}

func filterMapQ(config Config, reads GRanges) GRanges {
  if config.FilterMapQ <= 0 {
    return reads
  }
  PrintStderr(config, 1, "Filtering reads (minimum mapping quality: %d)... ", config.FilterMapQ)
  idx := []int{}
  if config.PairedEnd {
    mapq1 := reads.GetMetaInt("mapq1")
    mapq2 := reads.GetMetaInt("mapq2")
    if len(mapq1) == 0 || len(mapq2) == 0 {
      PrintStderr(config, 1, "failed\n")
      log.Fatalf("no mapping quality available")
    }
    for i := 0; i < reads.Length(); i++ {
      if mapq1[i] < config.FilterMapQ || mapq2[i] < config.FilterMapQ {
        idx = append(idx, i)
      }
    }
  } else {
    mapq := reads.GetMetaInt("score")
    if len(mapq) == 0 {
      mapq = reads.GetMetaInt("mapq")
      if len(mapq) == 0 {
        PrintStderr(config, 1, "failed\n")
        log.Fatalf("no mapping quality available")
      }
    }
    for i := 0; i < reads.Length(); i++ {
      if mapq[i] < config.FilterMapQ {
        idx = append(idx, i)
      }
    }
  }
  PrintStderr(config, 1, "done\n")
  PrintStderr(config, 1, "Filtered out %d reads (%.2f%%)\n", len(idx), 100.0*float64(len(idx))/float64(reads.Length()))
  return reads.Remove(idx)
}

func filterReadLength(config Config, reads GRanges) GRanges {
  if config.FilterReadLengths[0] == 0 && config.FilterReadLengths[1] == 0 {
    return reads
  }
  PrintStderr(config, 1, "Filtering reads (admissible read length: %v) ... ", config.FilterReadLengths)
  idx := []int{}
  for i := 0; i < reads.Length(); i++ {
    len := reads.Ranges[i].To - reads.Ranges[i].From
    if len < config.FilterReadLengths[0] ||
      (len > config.FilterReadLengths[1] && config.FilterReadLengths[1] != 0) {
      idx = append(idx, i)
    }
  }
  PrintStderr(config, 1, "done\n")
  PrintStderr(config, 1, "Filtered out %d reads (%.2f%%) with non-admissible length\n", len(idx), 100.0*float64(len(idx))/float64(reads.Length()))
  return reads.Remove(idx)
}

func shiftReads(config Config, reads GRanges) GRanges {
  if config.ShiftReads[0] == 0 && config.ShiftReads[1] == 0 {
    return reads
  }
  PrintStderr(config, 1, "Shifting reads reads (forward strand: %d, reverse strand: %d) ... ",
    config.ShiftReads[0], config.ShiftReads[1])
  reads = reads.Clone()
  for i := 0; i < reads.Length(); i++ {
    if reads.Strand[i] == '+' {
      reads.Ranges[i].From += config.ShiftReads[0]
      reads.Ranges[i].To   += config.ShiftReads[0]
    } else
    if reads.Strand[i] == '-' {
      reads.Ranges[i].From += config.ShiftReads[1]
      reads.Ranges[i].To   += config.ShiftReads[1]
    }
    if reads.Ranges[i].From < 0 {
      reads.Ranges[i].To   -= reads.Ranges[i].From
      reads.Ranges[i].From  = 0
    }
  }
  PrintStderr(config, 1, "done\n")
  return reads
}

/* fragment length estimation
 * -------------------------------------------------------------------------- */

func saveFraglen(config Config, filename string, fraglen int) {
  basename := strings.TrimRight(filename, filepath.Ext(filename))
  filename  = fmt.Sprintf("%s.fraglen.txt", basename)

  f, err := os.Create(filename)
  if err != nil {
    log.Fatalf("opening `%s' failed: %v", filename, err)
  }
  defer f.Close()

  fmt.Fprintf(f, "%d\n", fraglen)

  PrintStderr(config, 1, "Wrote fragment length estimate to `%s'", filename)
}

func saveCrossCorr(config Config, filename string, x []int, y []float64) {
  basename := strings.TrimRight(filename, filepath.Ext(filename))
  filename  = fmt.Sprintf("%s.fraglen.table", basename)

  f, err := os.Create(filename)
  if err != nil {
    log.Fatalf("opening `%s' failed: %v", filename, err)
  }
  defer f.Close()

  for i := 0; i < len(x); i++ {
    fmt.Fprintf(f, "%d %f", x[i], y[i])
  }
  PrintStderr(config, 1, "Wrote crosscorrelation to `%s'", filename)
}

func saveCrossCorrPlot(config Config, filename string, x []int, y []float64) {
  basename := strings.TrimRight(filename, filepath.Ext(filename))
  filename  = fmt.Sprintf("%s.fraglen.pdf", basename)

  xy := make(plotter.XYs, len(x))
  for i := 0; i < len(x); i++ {
    xy[i].X = float64(x[i])+1
    xy[i].Y = y[i]
  }
  p, err := plot.New()
  if err != nil {
    panic(err)
  }
  p.Title.Text = ""
  p.X.Label.Text = "shift"
  p.Y.Label.Text = "cross-correlation"
  p.X.Scale = plot.LogScale{}
  p.Y.Scale = plot.LogScale{}
  p.X.Tick.Marker = plot.LogTicks{}
  p.Y.Tick.Marker = plot.LogTicks{}

  err = plotutil.AddLines(p, xy)
  if err != nil {
    panic(err)
  }
  if err := p.Save(8*vg.Inch, 4*vg.Inch, filename); err != nil {
    panic(err)
  }
  PrintStderr(config, 1, "Wrote crosscorrelation plot to `%s'", filename)
}

func estimateFraglen(config Config, filename string, genome Genome) int {
  PrintStderr(config, 1, "Reading tags from `%s'... ", filename)
  reads := GRanges{}
  if err := reads.ImportBamSingleEnd(filename, BamReaderOptions{}); err != nil {
    PrintStderr(config, 1, "failed\n")
    log.Fatal(err)
  }
  PrintStderr(config, 1, "done\n")

  // first round of filtering
  reads = filterReadLength(config, reads)
  reads = filterDuplicates(config, reads)
  reads = filterMapQ(config, reads)

  // estimate fragment length
  PrintStderr(config, 1, "Estimating mean fragment length... ")
  if fraglen, x, y, err := EstimateFragmentLength(reads, genome, 2000, config.BinSize, config.FraglenRange); err != nil {
    PrintStderr(config, 1, "failed\n")
    log.Fatalf("estimating read length failed: %v", err)
    return 0
  } else {
    PrintStderr(config, 1, "done\n")
    PrintStderr(config, 1, "Estimated mean fragment length: %d\n", fraglen)

    if config.SaveFraglen {
      saveFraglen(config, filename, fraglen)
    }
    if config.SaveCrossCorr {
      saveCrossCorr(config, filename, x, y)
    }
    if config.SaveCrossCorrPlot {
      saveCrossCorrPlot(config, filename, x, y)
    }
    return fraglen
  }
}

/* -------------------------------------------------------------------------- */

func bamToBigWig(config Config, filenameTrack string, filenamesTreatment, filenamesControl []string, fraglenTreatment, fraglenControl []int, genome Genome) {

  // treatment data
  track1 := AllocSimpleTrack("treatment", genome, config.BinSize)

  // number of reads
  n_treatment := 0
  n_control   := 0

  for i, filename := range filenamesTreatment {
    fraglen := fraglenTreatment[i]

    PrintStderr(config, 1, "Reading treatment tags from `%s'... ", filename)
    treatment := GRanges{}
    if config.PairedEnd {
      if err := treatment.ImportBamPairedEnd(filename, BamReaderOptions{}); err != nil {
        PrintStderr(config, 1, "failed\n")
        log.Fatal(err)
      }
    } else {
      if err := treatment.ImportBamSingleEnd(filename, BamReaderOptions{}); err != nil {
        PrintStderr(config, 1, "failed\n")
        log.Fatal(err)
      }
    }
    PrintStderr(config, 1, "done\n")

    // first round of filtering
    treatment = filterReadLength(config, treatment)
    treatment = filterDuplicates(config, treatment)
    treatment = filterMapQ(config, treatment)
    // second round of filtering
    treatment = filterStrand(config, treatment)
    treatment = shiftReads(config, treatment)

    n_treatment += treatment.Length()

    switch config.BinningMethod {
    case "simple":
      GenericMutableTrack{track1}.AddReads(treatment, fraglen, false)
    case "overlap":
      GenericMutableTrack{track1}.AddReads(treatment, fraglen, true)
    default:
      log.Fatal("invalid binning method `%s'", config.BinningMethod)
    }
  }
  if config.NormalizeTrack == "rpm" {
    PrintStderr(config, 1, "Normalizing treatment track (rpm)... ")
    c := float64(1000000)/float64(n_treatment)
    GenericMutableTrack{track1}.Map(track1, func(name string, i int, x float64) float64 {
      return c*x
    })
    // adapt pseudocounts!
    config.Pseudocounts[0] *= c
    PrintStderr(config, 1, "done\n")
  }

  if len(filenamesControl) > 0 {
    // control data
    track2 := AllocSimpleTrack("control", genome, config.BinSize)

    for i, filename := range filenamesControl {
      fraglen  := fraglenControl[i]

      PrintStderr(config, 1, "Reading control tags from `%s'... ", filename)
      control := GRanges{}
      if config.PairedEnd {
        if err := control.ImportBamPairedEnd(filename, BamReaderOptions{}); err != nil {
          PrintStderr(config, 1, "failed\n")
          log.Fatal(err)
        }
      } else {
        if err := control.ImportBamSingleEnd(filename, BamReaderOptions{}); err != nil {
          PrintStderr(config, 1, "failed\n")
          log.Fatal(err)
        }
      }
      PrintStderr(config, 1, "done\n")

      // first round of filtering
      control = filterReadLength(config, control)
      control = filterDuplicates(config, control)
      control = filterMapQ(config, control)
      // second round of filtering
      control = filterStrand(config, control)
      control = shiftReads(config, control)

      n_control = control.Length()

      switch config.BinningMethod {
      case "simple":
        GenericMutableTrack{track2}.AddReads(control, fraglen, false)
      case "overlap":
        GenericMutableTrack{track2}.AddReads(control, fraglen, true)
      default:
        log.Fatal("invalid binning method `%s'", config.BinningMethod)
      }
    }
    if config.NormalizeTrack == "rpm" {
      PrintStderr(config, 1, "Normalizing control track (rpm)... ")
      c := float64(1000000)/float64(n_control)
      GenericMutableTrack{track2}.Map(track2, func(name string, i int, x float64) float64 {
        return c*x
      })
      // adapt pseudocounts!
      config.Pseudocounts[1] *= c
      PrintStderr(config, 1, "done\n")
    }
    if config.SmoothenControl {
      GenericMutableTrack{track2}.Smoothen(config.SmoothenMin, config.SmoothenSizes)
    }
    PrintStderr(config, 1, "Combining treatment and control tracks... ")
    if err := (GenericMutableTrack{track1}).Normalize(track1, track2, config.Pseudocounts[0], config.Pseudocounts[1], config.LogScale); err != nil {
      PrintStderr(config, 1, "failed\n")
      log.Fatalf("normalizing track failed: %v", err)
    }
    PrintStderr(config, 1, "done\n")
  } else {
    // no control data
    if config.Pseudocounts[0] != 0.0 {
      PrintStderr(config, 1, "Adding pseudocount `%f'... ", config.Pseudocounts[0])
      GenericMutableTrack{track1}.Map(track1, func(name string, i int, x float64) float64 { return x+config.Pseudocounts[0] })
      PrintStderr(config, 1, "done\n")
    }
    if config.LogScale {
      PrintStderr(config, 1, "Log-transforming data... ")
      GenericMutableTrack{track1}.Map(track1, func(name string, i int, x float64) float64 { return math.Log(x) })
      PrintStderr(config, 1, "done\n")
    }
  }
  PrintStderr(config, 1, "Writing track `%s'... ", filenameTrack)
  parameters := DefaultBigWigParameters()
  parameters.ReductionLevels = config.BWZoomLevels
  if err := (GenericTrack{track1}).ExportBigWig(filenameTrack, genome, parameters); err != nil {
    PrintStderr(config, 1, "failed\n")
  } else {
    PrintStderr(config, 1, "done\n")
  }
}

/* -------------------------------------------------------------------------- */

func main() {

  var config Config

  options := getopt.New()

  // bigWig options
  optBWZoomLevels      := options. StringLong("bigwig-zoom-levels",         0 , "", "comma separated list of BigWig zoom levels")
  // read options
  optShiftReads        := options. StringLong("shift-reads",                0 , "", "shift reads on the positive strand by `x' bps and those on the negative strand by `y' bps [format: x,y]")
  optPairedEnd         := options.   BoolLong("paired-end",                 0 ,     "reads are paired-end")
  // options for filterering reads
  optFilterStrand      := options. StringLong("filter-strand",              0 , "", "use reads on either the forward `+' or reverse `-' strand")
  optReadLength        := options. StringLong("filter-read-lengths",        0 , "", "feasible range of read-lengths [format: min:max]")
  optFilterMapQ        := options.    IntLong("filter-mapq",                0 ,  0, "filter reads for minimum mapping quality (default: 0)")
  optFilterDuplicates  := options.   BoolLong("filter-duplicates",          0 ,     "remove reads marked as duplicates")
  // track options
  optBinningMethod     := options. StringLong("binning-method",             0 , "", "binning method (i.e. simple or overlap [default])")
  optBinSize           := options.    IntLong("bin-size",                   0 , -1, "track bin size")
  optNormalizeTrack    := options. StringLong("normalize-track",            0 , "", "normalize track with the specified method (i.e. rpm)")
  optPseudocounts      := options. StringLong("pseudocounts",               0 , "", "pseudocounts added to treatment and control signal (default: `0,0')")
  optSmoothenControl   := options.   BoolLong("smoothen-control",           0 ,     "smoothen control with an adaptive window method")
  optSmoothenSizes     := options. StringLong("smoothen-window-sizes",      0 , "", "feasible window sizes for the smoothening method [format: s1,s2,...]")
  optSmoothenMin       := options. StringLong("smoothen-min-counts",        0 , "", "minimum number of counts for the smoothening method")
  optLogScale          := options.   BoolLong("log-scale",                  0 ,     "log-transform data")
  // options for estimating and setting fragment lengths
  optFraglen           := options.    IntLong("fragment-length",            0 , -1, "fragment length for all input files (reads are extended to the given length)")
  optFraglenRange      := options. StringLong("fragment-length-range",      0 , "", "feasible range of fragment lengths (format from:to)")
  optEstimateFraglen   := options.   BoolLong("estimate-fragment-length",   0 ,     "use crosscorrelation to estimate the fragment length")
  optSaveFraglen       := options.   BoolLong("save-fraglen",               0 ,     "save estimated fragment length in a file named <BAM_BASENAME>.fraglen.txt")
  optSaveCrossCorr     := options.   BoolLong("save-crosscorrelation",      0 ,     "save crosscorrelation between forward and reverse strands in a file named <BAM_BASENAME>.fraglen.table")
  optSaveCrossCorrPlot := options.   BoolLong("save-crosscorrelation-plot", 0 ,     "save crosscorrelation plot in a file names <BAM_BASENAME>.fraglen.pdf")
  // generic options
  optVerbose           := options.CounterLong("verbose",                   'v',     "verbose level [-v or -vv]")
  optHelp              := options.   BoolLong("help",                      'h',     "print help")

  options.SetParameters("<TREATMENT1.bam[:FRAGLEN],TREATMENT2.bam[:FRAGLEN],...> <CONTROL1.bam[:FRAGLEN],CONTROL2.bam[:FRAGLEN],...> <RESULT.bw>")
  options.Parse(os.Args)

  // parse options
  //////////////////////////////////////////////////////////////////////////////
  if *optHelp {
    options.PrintUsage(os.Stdout)
    os.Exit(0)
  }
  if *optVerbose != 0 {
    config.Verbose = *optVerbose
  }
  if len(options.Args()) != 2 && len(options.Args()) != 3 {
    options.PrintUsage(os.Stderr)
    os.Exit(1)
  }
  if *optBinSize != -1 {
    if *optBinSize < 1 {
      options.PrintUsage(os.Stderr)
      os.Exit(1)
    }
    config.BinSize = *optBinSize
  }
  if *optBinningMethod != "" {
    config.BinningMethod = *optBinningMethod
  }
  if *optPseudocounts != "" {
    tmp := strings.Split(*optPseudocounts, ",")
    if len(tmp) != 2 {
      options.PrintUsage(os.Stderr)
      os.Exit(1)
    }
    t1, err := strconv.ParseFloat(tmp[0], 64)
    if err != nil {
      log.Fatal(err)
    }
    t2, err := strconv.ParseFloat(tmp[1], 64)
    if err != nil {
      log.Fatal(err)
    }
    config.Pseudocounts[0] = t1
    config.Pseudocounts[1] = t2
  }
  if *optLogScale {
    config.LogScale = true
  }
  if *optFraglen != -1 {
    if *optFraglen < 0 {
      options.PrintUsage(os.Stderr)
      os.Exit(1)
    }
    config.Fraglen = *optFraglen
  }
  if *optPairedEnd && *optEstimateFraglen {
    log.Fatal("cannot estimate fragment length for paired-end reads")
  }
  if *optReadLength != "" {
    tmp := strings.Split(*optReadLength, ":")
    if len(tmp) != 2 {
      options.PrintUsage(os.Stderr)
      os.Exit(1)
    }
    t1, err := strconv.ParseInt(tmp[0], 10, 64)
    if err != nil {
      log.Fatal(err)
    }
    t2, err := strconv.ParseInt(tmp[1], 10, 64)
    if err != nil {
      log.Fatal(err)
    }
    if t1 > t2 {
      options.PrintUsage(os.Stderr)
      os.Exit(1)
    }
    config.FilterReadLengths[0] = int(t1)
    config.FilterReadLengths[1] = int(t2)
  }
  if *optFilterMapQ < 0 {
      options.PrintUsage(os.Stderr)
      os.Exit(1)
  } else {
    config.FilterMapQ = *optFilterMapQ
  }
  if *optFilterStrand != "" {
    switch *optFilterStrand {
    case "+": config.FilterStrand = '+'
    case "-": config.FilterStrand = '-'
    default:
      options.PrintUsage(os.Stderr)
      os.Exit(1)
    }
  }
  if *optShiftReads != "" {
    tmp := strings.Split(*optShiftReads, ",")
    if len(tmp) != 2 {
      options.PrintUsage(os.Stderr)
      os.Exit(1)
    }
    t1, err := strconv.ParseInt(tmp[0], 10, 64)
    if err != nil {
      log.Fatal(err)
    }
    t2, err := strconv.ParseInt(tmp[1], 10, 64)
    if err != nil {
      log.Fatal(err)
    }
    config.ShiftReads[0] = int(t1)
    config.ShiftReads[1] = int(t2)
  }
  if *optSmoothenControl {
    config.SmoothenControl = true
  }
  if *optSmoothenSizes != "" {
    config.SmoothenSizes = []int{}
    tmp := strings.Split(*optSmoothenSizes, ",")
    for i := 0; i < len(tmp); i++ {
      t, err := strconv.ParseInt(tmp[i], 10, 64)
      if err != nil {
        log.Fatal(err)
      }
      config.SmoothenSizes = append(config.SmoothenSizes, int(t))
    }
  }
  if *optSmoothenMin != "" {
    t, err := strconv.ParseFloat(*optSmoothenMin, 64)
    if err != nil {
      log.Fatal(err)
    }
    config.SmoothenMin = t
  }
  if *optBWZoomLevels != "" {
    tmp := strings.Split(*optBWZoomLevels, ",")
    config.BWZoomLevels = []int{}
    for i := 0; i < len(tmp); i++ {
      if t, err := strconv.ParseInt(tmp[i], 10, 64); err != nil {
        log.Fatal(err)
      } else {
        config.BWZoomLevels = append(config.BWZoomLevels, int(t))
      }
    }
  }
  if *optNormalizeTrack != "" {
    switch *optNormalizeTrack {
    case "rpm":
    default:
      log.Fatal("invalid normalization method `%s'", *optNormalizeTrack)
    }
    config.NormalizeTrack = *optNormalizeTrack
  }
  if *optFraglenRange != "" {
    tmp := strings.Split(*optFraglenRange, ":")
    if len(tmp) != 2 {
      options.PrintUsage(os.Stderr)
      os.Exit(1)
    }
    t1, err := strconv.ParseInt  (tmp[0], 10, 64)
    if err != nil {
      log.Fatalf("parsing fragment length range failed: %v", err)
    }
    t2, err := strconv.ParseInt  (tmp[1], 10, 64)
    if err != nil {
      log.Fatalf("parsing fragment length range failed: %v", err)
    }
    config.FraglenRange[0] = int(t1)
    config.FraglenRange[1] = int(t2)
  }
  config.FilterDuplicates  = *optFilterDuplicates
  config.PairedEnd         = *optPairedEnd
  config.SaveFraglen       = *optSaveFraglen
  config.SaveCrossCorr     = *optSaveCrossCorr
  config.SaveCrossCorrPlot = *optSaveCrossCorrPlot

  // parse arguments
  //////////////////////////////////////////////////////////////////////////////
  filenamesTreatment := strings.Split(options.Args()[0], ",")
  filenamesControl   := []string{}
  filenameTrack      := ""
  if len(options.Args()) == 3 {
    filenamesControl = strings.Split(options.Args()[1], ",")
    filenameTrack = options.Args()[2]
  } else {
    filenameTrack = options.Args()[1]
  }

  fraglenTreatment := make([]int, len(filenamesTreatment))
  fraglenControl   := make([]int, len(filenamesControl))
  // split filename:fraglen
  for i, filename := range filenamesTreatment {
    filenamesTreatment[i], fraglenTreatment[i] = parseFilename(filename)
  }
  for i, filename := range filenamesControl {
    filenamesControl[i], fraglenControl[i] = parseFilename(filename)
  }

  // read genome
  //////////////////////////////////////////////////////////////////////////////
  var genome Genome

  for _, filename := range append(filenamesTreatment, filenamesControl...) {
    g, err := BamImportGenome(filename); if err != nil {
      log.Fatal(err)
    }
    if genome.Length() == 0 {
      genome = g
    } else {
      if !genome.Equals(g) {
        log.Fatal("bam genomes are not equal")
      }
    }
  }

  // fragment length estimation
  //////////////////////////////////////////////////////////////////////////////
  if *optEstimateFraglen {
    for i, filename := range filenamesTreatment {
      if fraglenTreatment[i] != 0 {
        continue
      }
      fraglenTreatment[i] = estimateFraglen(config, filename, genome)
    }
    for i, filename := range filenamesControl {
      if fraglenControl[i] != 0 {
        continue
      }
      fraglenControl[i] = estimateFraglen(config, filename, genome)
    }
  }

  // bam -> bigWig
  //////////////////////////////////////////////////////////////////////////////
  bamToBigWig(config, filenameTrack, filenamesTreatment, filenamesControl, fraglenTreatment, fraglenControl, genome)
}
