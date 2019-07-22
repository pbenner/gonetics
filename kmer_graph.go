/* Copyright (C) 2019 Philipp Benner
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

type KmerGraph struct {
  catalogue KmerCatalogue
  nodes     map[KmerClassId]*KmerGraphNode
}

/* -------------------------------------------------------------------------- */

type KmerGraphNode struct {
  Kmer KmerClass
  // shorter kmers
  Intra []*KmerGraphNode
  // longer kmers
  Extra []*KmerGraphNode
  // matching kmers (those with more uncertain entries)
  Infra []*KmerGraphNode
  // instances of kmers (those with less uncertain entries)
  Supra []*KmerGraphNode
}

/* -------------------------------------------------------------------------- */

func (obj *KmerGraph) RelatedKmers(kmer string) KmerClassList {
  cl, ok := obj.catalogue.GetKmerClassIfPresent(kmer)
  if !ok {
    return nil
  }
  if node, ok := obj.nodes[cl.KmerClassId]; !ok {
    return nil
  } else {
    r := KmerClassList(nil)
    for _, n := range node.Intra {
      r = append(r, n.Kmer)
    }
    for _, n := range node.Extra {
      r = append(r, n.Kmer)
    }
    for _, n := range node.Infra {
      r = append(r, n.Kmer)
    }
    for _, n := range node.Supra {
      r = append(r, n.Kmer)
    }
    return r
  }
}

/* -------------------------------------------------------------------------- */

func (obj *KmerGraph) newNode(kmer KmerClass) *KmerGraphNode {
  if node, ok := obj.nodes[kmer.KmerClassId]; ok {
    return node
  }
  r := KmerGraphNode{}
  m := kmer.CountAmbiguous(obj.catalogue.al)
  if m == 0 {
    // this k-mer has no ambiguous positions, fill
    // intra, extra, and infra k-mers
    r.Intra = obj.computeIntra(kmer)
    r.Extra = obj.computeExtra(kmer)
    //r.Infra = obj.computeInfra(kmer)
  } else {
    // this k-mer has ambiguous entries, fill
    // infra and supra k-mers
    //r.Infra = obj.computeInfra(kmer)
    //r.Supra = obj.computeSupra(kmer)
  }
  obj.nodes[kmer.KmerClassId] = &r
  return &r
}

/* -------------------------------------------------------------------------- */

func (obj *KmerGraph) computeIntra(kmer KmerClass) []*KmerGraphNode {
  if obj.catalogue.n > kmer.K-1 {
    return nil
  }
  // two possible intra-k-mers
  s1 := kmer.Elements[0][0:kmer.K-1]
  s2 := kmer.Elements[0][1:kmer.K]
  r  := []*KmerGraphNode(nil)
  if cl, ok := obj.catalogue.GetKmerClassIfPresent(s1); ok {
    r = append(r, obj.newNode(cl))
  }
  if cl, ok := obj.catalogue.GetKmerClassIfPresent(s2); ok {
    r = append(r, obj.newNode(cl))
  }
  return r
}

func (obj *KmerGraph) computeExtra(kmer KmerClass) []*KmerGraphNode {
  if obj.catalogue.m < kmer.K+1 {
    return nil
  }
  r   := []*KmerGraphNode(nil)
  it1 := NewKmerCylinderIterator(kmer.K+1, 0, obj.catalogue.al, 0, kmer.Elements[0])
  it2 := NewKmerCylinderIterator(kmer.K+1, 0, obj.catalogue.al, 1, kmer.Elements[0])
  for ; it1.Ok(); it1.Next() {
    if cl, ok := obj.catalogue.GetKmerClassIfPresent(it1.Get()); ok {
      r = append(r, obj.newNode(cl))
    }
  }
  for ; it2.Ok(); it2.Next() {
    if cl, ok := obj.catalogue.GetKmerClassIfPresent(it2.Get()); ok {
      r = append(r, obj.newNode(cl))
    }
  }
  return r
}
