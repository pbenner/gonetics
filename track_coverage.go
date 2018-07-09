/* Copyright (C) 2016-2018 Philipp Benner
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

import   "fmt"
import   "io"
import   "io/ioutil"
import   "math"

/* -------------------------------------------------------------------------- */

type OptionLogger struct {
  Value io.Writer
}

type OptionBinningMethod struct {
  Value string
}

type OptionBinSize struct {
  Value int
}

type OptionBinOverlap struct {
  Value int
}

type OptionNormalizeTrack struct {
  Value string
}

type OptionShiftReads struct {
  Value [2]int
}

type OptionPairedAsSingleEnd struct {
  Value bool
}

type OptionPairedEndStrandSpecific struct {
  Value bool
}

type OptionLogScale struct {
  Value bool
}

type OptionPseudocounts struct {
  Value [2]float64
}

type OptionEstimateFraglen struct {
  Value bool
}

type OptionFraglenRange struct {
  Value [2]int
}

type OptionFraglenBinSize struct {
  Value int
}

type OptionFilterChroms struct {
  Value []string
}

type OptionFilterMapQ struct {
  Value int
}

type OptionFilterReadLengths struct {
  Value [2]int
}

type OptionFilterDuplicates struct {
  Value bool
}

type OptionFilterStrand struct {
  Value byte
}

type OptionFilterPairedEnd struct {
  Value bool
}

type OptionFilterSingleEnd struct {
  Value bool
}

type OptionSmoothenControl struct {
  Value bool
}

type OptionSmoothenSizes struct {
  Value []int
}

type OptionSmoothenMin struct {
  Value float64
}

/* -------------------------------------------------------------------------- */

type BamCoverageConfig struct {
  Logger               io.Writer
  BinningMethod           string
  BinSize                 int
  BinOverlap              int
  NormalizeTrack          string
  ShiftReads           [2]int
  PairedAsSingleEnd       bool
  PairedEndStrandSpecific bool
  LogScale                bool
  Pseudocounts         [2]float64
  EstimateFraglen         bool
  FraglenRange         [2]int
  FraglenBinSize          int
  FilterChroms          []string
  FilterMapQ              int
  FilterReadLengths    [2]int
  FilterDuplicates        bool
  FilterStrand            byte
  FilterPairedEnd         bool
  FilterSingleEnd         bool
  SmoothenControl         bool
  SmoothenSizes         []int
  SmoothenMin             float64
}

func BamCoverageDefaultConfig() BamCoverageConfig {
  config := BamCoverageConfig{}
  // set default values
  config.Logger                  = ioutil.Discard
  config.BinningMethod           = "simple"
  config.BinSize                 = 10
  config.BinOverlap              = 0
  config.PairedAsSingleEnd       = false
  config.PairedEndStrandSpecific = false
  config.EstimateFraglen         = false
  config.FraglenRange            = [2]int{-1, -1}
  config.FraglenBinSize          = 10
  config.FilterReadLengths       = [2]int{0,0}
  config.FilterMapQ              = 0
  config.FilterDuplicates        = false
  config.FilterStrand            = '*'
  config.FilterPairedEnd         = false
  config.FilterSingleEnd         = false
  config.LogScale                = false
  config.Pseudocounts            = [2]float64{0.0, 0.0}
  config.SmoothenControl         = false
  config.SmoothenSizes           = []int{}
  config.SmoothenMin             = 20.0
  return config
}

/* -------------------------------------------------------------------------- */

type fraglenEstimate struct {
  Fraglen   int
  X       []int
  Y       []float64
}

/* read filters
 * -------------------------------------------------------------------------- */

// treat all paired end reads as single end reads, this allows
// to extend/crop paired end reads when adding them to the track
// with AddReads()
func filterPairedAsSingleEnd(config BamCoverageConfig, chanIn ReadChannel) ReadChannel {
  if config.PairedAsSingleEnd == false {
    return chanIn
  }
  chanOut := make(chan Read)
  go func() {
    for r := range chanIn {
      r.PairedEnd = false; chanOut <- r
    }
    close(chanOut)
  }()
  return chanOut
}

func filterPairedEnd(config BamCoverageConfig, chanIn ReadChannel) ReadChannel {
  if config.FilterPairedEnd == false {
    return chanIn
  }
  chanOut := make(chan Read)
  go func() {
    n := 0
    m := 0
    for r := range chanIn {
      if r.PairedEnd {
        chanOut <- r; m++
      }
      n++
    }
    if n != 0 {
      fmt.Fprintf(config.Logger, "Filtered out %d unpaired reads (%.2f%%)\n", n-m, 100.0*float64(n-m)/float64(n))
    }
    close(chanOut)
  }()
  return chanOut
}

func filterSingleEnd(config BamCoverageConfig, veto bool, chanIn ReadChannel) ReadChannel {
  if config.FilterSingleEnd == false && !veto {
    return chanIn
  }
  chanOut := make(chan Read)
  go func() {
    n := 0
    m := 0
    for r := range chanIn {
      if !r.PairedEnd {
        chanOut <- r; m++
      }
      n++
    }
    if n != 0 {
      fmt.Fprintf(config.Logger, "Filtered out %d paired reads (%.2f%%)\n", n-m, 100.0*float64(n-m)/float64(n))
    }
    close(chanOut)
  }()
  return chanOut
}

func filterDuplicates(config BamCoverageConfig, chanIn ReadChannel) ReadChannel {
  if config.FilterDuplicates == false {
    return chanIn
  }
  chanOut := make(chan Read)
  go func() {
    n := 0
    m := 0
    for r := range chanIn {
      if !r.Duplicate {
        chanOut <- r; m++
      }
      n++
    }
    if n != 0 {
      fmt.Fprintf(config.Logger, "Filtered out %d duplicates (%.2f%%)\n", n-m, 100.0*float64(n-m)/float64(n))
    }
    close(chanOut)
  }()
  return chanOut
}

func filterStrand(config BamCoverageConfig, chanIn ReadChannel) ReadChannel {
  if config.FilterStrand == '*' {
    return chanIn
  }
  chanOut := make(chan Read)
  go func() {
    n := 0
    m := 0
    for r := range chanIn {
      if r.Strand == config.FilterStrand {
        chanOut <- r; m++
      }
      n++
    }
    if n != 0 {
      fmt.Fprintf(config.Logger, "Filtered out %d reads not on strand %c (%.2f%%)\n", n-m, config.FilterStrand, 100.0*float64(n-m)/float64(n))
    }
    close(chanOut)
  }()
  return chanOut
}

func filterMapQ(config BamCoverageConfig, chanIn ReadChannel) ReadChannel {
  if config.FilterMapQ <= 0 {
    return chanIn
  }
  chanOut := make(chan Read)
  go func() {
    n := 0
    m := 0
    for r := range chanIn {
      if r.MapQ >= config.FilterMapQ {
        chanOut <- r; m++
      }
      n++
    }
    if n != 0 {
      fmt.Fprintf(config.Logger, "Filtered out %d reads with mapping quality lower than %d (%.2f%%)\n", n-m, config.FilterMapQ, 100.0*float64(n-m)/float64(n))
    }
    close(chanOut)
  }()
  return chanOut
}

func filterReadLength(config BamCoverageConfig, chanIn ReadChannel) ReadChannel {
  if config.FilterReadLengths[0] == 0 && config.FilterReadLengths[1] == 0 {
    return chanIn
  }
  chanOut := make(chan Read)
  go func() {
    n := 0
    m := 0
    for r := range chanIn {
      len := r.Range.To - r.Range.From
      if len >= config.FilterReadLengths[0] &&
        (len <= config.FilterReadLengths[1] || config.FilterReadLengths[1] == 0) {
        chanOut <- r; m++
      }
      n++
    }
    if n != 0 {
      fmt.Fprintf(config.Logger, "Filtered out %d reads with non-admissible length (%.2f%%)\n", n-m, 100.0*float64(n-m)/float64(n))
    }
    close(chanOut)
  }()
  return chanOut
}

func shiftReads(config BamCoverageConfig, chanIn ReadChannel) ReadChannel {
  if config.ShiftReads[0] == 0 && config.ShiftReads[1] == 0 {
    return chanIn
  }
  chanOut := make(chan Read)
  go func() {
    for r := range chanIn {
      if r.Strand == '+' {
        r.Range.From += config.ShiftReads[0]
        r.Range.To   += config.ShiftReads[0]
      } else
      if r.Strand == '-' {
        r.Range.From += config.ShiftReads[1]
        r.Range.To   += config.ShiftReads[1]
      }
      if r.Range.From < 0 {
        r.Range.To   -= r.Range.From
        r.Range.From  = 0
      }
      chanOut <- r
    }
    fmt.Fprintf(config.Logger, "Shifted reads (forward strand: %d, reverse strand: %d)\n",
      config.ShiftReads[0], config.ShiftReads[1])
    close(chanOut)
  }()
  return chanOut
}

/* fragment length estimation
 * -------------------------------------------------------------------------- */

func estimateFraglen(config BamCoverageConfig, filename string, genome Genome) (fraglenEstimate, error) {
  var reads ReadChannel

  fmt.Fprintf(config.Logger, "Reading tags from `%s'...\n", filename)
  if bam, err := OpenBamFile(filename, BamReaderOptions{}); err != nil {
    fmt.Fprintf(config.Logger, "failed\n")
    return fraglenEstimate{}, err
  } else {
    defer bam.Close()
    reads = bam.ReadSimple(false, false)
  }

  // first round of filtering
  reads = filterSingleEnd(config, true, reads)
  reads = filterReadLength(config, reads)
  reads = filterDuplicates(config, reads)
  reads = filterMapQ(config, reads)

  // estimate fragment length
  fmt.Fprintf(config.Logger, "Estimating mean fragment length...\n")
  if fraglen, x, y, err := EstimateFragmentLength(reads, genome, 2000, config.FraglenBinSize, config.FraglenRange); err != nil {
    fmt.Fprintf(config.Logger, "failed\n")
    return fraglenEstimate{0, x, y}, err
  } else {
    fmt.Fprintf(config.Logger, "Estimated mean fragment length: %d\n", fraglen)
    return fraglenEstimate{fraglen, x, y}, nil
  }
}

/* -------------------------------------------------------------------------- */

func bamCoverage(config BamCoverageConfig, filenameTrack string, filenamesTreatment, filenamesControl []string, fraglenTreatment, fraglenControl []int, genome Genome) (SimpleTrack, error) {

  // treatment data
  track1 := AllocSimpleTrack("treatment", genome, config.BinSize)

  // number of reads
  n_treatment := 0
  n_control   := 0

  for i, filename := range filenamesTreatment {
    fraglen := fraglenTreatment[i]

    var treatment ReadChannel
    fmt.Fprintf(config.Logger, "Reading treatment tags from `%s'...\n", filename)
    if bam, err := OpenBamFile(filename, BamReaderOptions{}); err != nil {
      fmt.Fprintf(config.Logger, "failed\n")
      return SimpleTrack{}, err
    } else {
      defer bam.Close()
      treatment = bam.ReadSimple(!config.PairedAsSingleEnd, config.PairedEndStrandSpecific)
    }

    // first round of filtering
    treatment = filterPairedEnd(config, treatment)
    treatment = filterSingleEnd(config, false, treatment)
    treatment = filterPairedAsSingleEnd(config, treatment)
    treatment = filterReadLength(config, treatment)
    treatment = filterDuplicates(config, treatment)
    treatment = filterMapQ(config, treatment)
    // second round of filtering
    treatment = filterStrand(config, treatment)
    treatment = shiftReads(config, treatment)

    n_treatment += GenericMutableTrack{track1}.AddReads(treatment, fraglen, config.BinningMethod)
  }
  if config.NormalizeTrack == "rpkm" {
    fmt.Fprintf(config.Logger, "Normalizing treatment track (rpkm)... ")
    c := float64(1000000)/(float64(n_treatment)*float64(config.BinSize))
    GenericMutableTrack{track1}.Map(track1, func(name string, i int, x float64) float64 {
      return c*x
    })
    // adapt pseudocounts!
    config.Pseudocounts[0] *= c
    fmt.Fprintf(config.Logger, "done\n")
  }
  if config.NormalizeTrack == "cpm" {
    fmt.Fprintf(config.Logger, "Normalizing treatment track (cpm)... ")
    c := float64(1000000)/float64(n_treatment)
    GenericMutableTrack{track1}.Map(track1, func(name string, i int, x float64) float64 {
      return c*x
    })
    // adapt pseudocounts!
    config.Pseudocounts[0] *= c
    fmt.Fprintf(config.Logger, "done\n")
  }

  if len(filenamesControl) > 0 {
    // control data
    track2 := AllocSimpleTrack("control", genome, config.BinSize)

    for i, filename := range filenamesControl {
      fraglen  := fraglenControl[i]

      var control ReadChannel
      fmt.Fprintf(config.Logger, "Reading treatment tags from `%s'...\n", filename)
      if bam, err := OpenBamFile(filename, BamReaderOptions{}); err != nil {
        fmt.Fprintf(config.Logger, "failed\n")
        return SimpleTrack{}, err
      } else {
        defer bam.Close()
        control = bam.ReadSimple(!config.PairedAsSingleEnd, config.PairedEndStrandSpecific)
      }

      // first round of filtering
      control = filterPairedEnd(config, control)
      control = filterSingleEnd(config, false, control)
      control = filterPairedAsSingleEnd(config, control)
      control = filterReadLength(config, control)
      control = filterDuplicates(config, control)
      control = filterMapQ(config, control)
      // second round of filtering
      control = filterStrand(config, control)
      control = shiftReads(config, control)

      n_control += GenericMutableTrack{track2}.AddReads(control, fraglen, config.BinningMethod)
    }
    if config.NormalizeTrack == "rpkm" {
      fmt.Fprintf(config.Logger, "Normalizing control track (rpkm)... ")
      c := float64(1000000)/(float64(n_control)*float64(config.BinSize))
      GenericMutableTrack{track2}.Map(track2, func(name string, i int, x float64) float64 {
        return c*x
      })
      // adapt pseudocounts!
      config.Pseudocounts[1] *= c
      fmt.Fprintf(config.Logger, "done\n")
    }
    if config.NormalizeTrack == "cpm" {
      fmt.Fprintf(config.Logger, "Normalizing control track (cpm)... ")
      c := float64(1000000)/float64(n_control)
      GenericMutableTrack{track2}.Map(track2, func(name string, i int, x float64) float64 {
        return c*x
      })
      // adapt pseudocounts!
      config.Pseudocounts[1] *= c
      fmt.Fprintf(config.Logger, "done\n")
    }
    if config.SmoothenControl {
      GenericMutableTrack{track2}.Smoothen(config.SmoothenMin, config.SmoothenSizes)
    }
    fmt.Fprintf(config.Logger, "Combining treatment and control tracks... ")
    if err := (GenericMutableTrack{track1}).Normalize(track1, track2, config.Pseudocounts[0], config.Pseudocounts[1], config.LogScale); err != nil {
      fmt.Fprintf(config.Logger, "failed\n")
      return SimpleTrack{}, err
    }
    fmt.Fprintf(config.Logger, "done\n")
  } else {
    // no control data
    if config.Pseudocounts[0] != 0.0 {
      fmt.Fprintf(config.Logger, "Adding pseudocount `%f'... ", config.Pseudocounts[0])
      GenericMutableTrack{track1}.Map(track1, func(name string, i int, x float64) float64 { return x+config.Pseudocounts[0] })
      fmt.Fprintf(config.Logger, "done\n")
    }
    if config.LogScale {
      fmt.Fprintf(config.Logger, "Log-transforming data... ")
      GenericMutableTrack{track1}.Map(track1, func(name string, i int, x float64) float64 { return math.Log(x) })
      fmt.Fprintf(config.Logger, "done\n")
    }
  }
  if len(config.FilterChroms) != 0 {
    fmt.Fprintf(config.Logger, "Removing all reads from `%v'... ", config.FilterChroms)
    for _, chr := range config.FilterChroms {
      if s, err := track1.GetMutableSequence(chr); err == nil {
        for i := 0; i < s.NBins(); i++ {
          s.SetBin(i, 0.0)
        }
      }
    }
    fmt.Fprintf(config.Logger, "done\n")
  }
  return track1, nil
}

/* -------------------------------------------------------------------------- */

func BamCoverage(filenameTrack string, filenamesTreatment, filenamesControl []string, fraglenTreatment, fraglenControl []int, options ...interface{}) (SimpleTrack, []fraglenEstimate, []fraglenEstimate, error) {

  config := BamCoverageDefaultConfig()

  // parse options
  //////////////////////////////////////////////////////////////////////////////

  for _, option := range options {
    switch opt := option.(type) {
    case OptionLogger:
      config.Logger = opt.Value
    case OptionBinningMethod:
      config.BinningMethod = opt.Value
    case OptionBinSize:
      config.BinSize = opt.Value
    case OptionBinOverlap:
      config.BinOverlap = opt.Value
    case OptionNormalizeTrack:
      config.NormalizeTrack = opt.Value
    case OptionShiftReads:
      config.ShiftReads = opt.Value
    case OptionPairedAsSingleEnd:
      config.PairedAsSingleEnd = opt.Value
    case OptionPairedEndStrandSpecific:
      config.PairedEndStrandSpecific = opt.Value
    case OptionLogScale:
      config.LogScale = opt.Value
    case OptionPseudocounts:
      config.Pseudocounts = opt.Value
    case OptionEstimateFraglen:
      config.EstimateFraglen = opt.Value
    case OptionFraglenRange:
      config.FraglenRange = opt.Value
    case OptionFraglenBinSize:
      config.FraglenBinSize = opt.Value
    case OptionFilterChroms:
      config.FilterChroms = opt.Value
    case OptionFilterMapQ:
      config.FilterMapQ = opt.Value
    case OptionFilterReadLengths:
      config.FilterReadLengths = opt.Value
    case OptionFilterDuplicates:
      config.FilterDuplicates = opt.Value
    case OptionFilterStrand:
      config.FilterStrand = opt.Value
    case OptionFilterPairedEnd:
      config.FilterPairedEnd = opt.Value
    case OptionFilterSingleEnd:
      config.FilterSingleEnd = opt.Value
    case OptionSmoothenControl:
      config.SmoothenControl = opt.Value
    case OptionSmoothenSizes:
      config.SmoothenSizes = opt.Value
    case OptionSmoothenMin:
      config.SmoothenMin = opt.Value
    default:
      return SimpleTrack{}, nil, nil, fmt.Errorf("BamCoverage(): invalid option: %v", opt)
    }
  }

  // read genome
  //////////////////////////////////////////////////////////////////////////////
  var genome Genome

  for _, filename := range append(filenamesTreatment, filenamesControl...) {
    g, err := BamImportGenome(filename); if err != nil {
      return SimpleTrack{}, nil, nil, err
    }
    if genome.Length() == 0 {
      genome = g
    } else {
      if !genome.Equals(g) {
        return SimpleTrack{}, nil, nil, fmt.Errorf("bam genomes are not equal")
      }
    }
  }

  treatmentFraglenEstimates := []fraglenEstimate{}
    controlFraglenEstimates := []fraglenEstimate{}

  // check fraglen arguments
  //////////////////////////////////////////////////////////////////////////////
  if len(fraglenTreatment) == 0 {
    fraglenTreatment = make([]int, len(filenamesTreatment))
  }
  if len(fraglenControl) == 0 {
    fraglenControl = make([]int, len(filenamesControl))
  }

  // fragment length estimation
  //////////////////////////////////////////////////////////////////////////////
  if config.EstimateFraglen {
    for i, filename := range filenamesTreatment {
      if fraglenTreatment[i] != 0 {
        continue
      }
      if estimate, err := estimateFraglen(config, filename, genome); err != nil {
        return SimpleTrack{}, treatmentFraglenEstimates, controlFraglenEstimates, fmt.Errorf("Estimating fragment length for `%s' failed: %v", filename, err)
      } else {
        treatmentFraglenEstimates = append(treatmentFraglenEstimates, estimate)
        fraglenTreatment[i] = estimate.Fraglen
      }
    }
    for i, filename := range filenamesControl {
      if fraglenControl[i] != 0 {
        continue
      }
      if estimate, err := estimateFraglen(config, filename, genome); err != nil {
        return SimpleTrack{}, controlFraglenEstimates, controlFraglenEstimates, fmt.Errorf("Estimating fragment length for `%s' failed: %v", filename, err)
      } else {
        controlFraglenEstimates = append(controlFraglenEstimates, estimate)
        fraglenControl[i] = estimate.Fraglen
      }
    }
  }
  //////////////////////////////////////////////////////////////////////////////
  if result, err := bamCoverage(config, filenameTrack, filenamesTreatment, filenamesControl, fraglenTreatment, fraglenControl, genome); err != nil {
    return SimpleTrack{}, treatmentFraglenEstimates, controlFraglenEstimates, err
  } else {
    return result, treatmentFraglenEstimates, controlFraglenEstimates, err
  }
}
