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

func (obj *KmerGraph) GetNode(kmer string) *KmerGraphNode {
  cl, ok := obj.catalogue.GetKmerClassIfPresent(kmer)
  if !ok {
    return nil
  }
  if node, ok := obj.nodes[cl.KmerClassId]; !ok {
    return nil
  } else {
    return node
  }
}

func (obj *KmerGraph) RelatedKmers(kmer string) KmerClassList {
  if node := obj.GetNode(kmer); node == nil {
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

func (obj *KmerGraph) constructGraph() {
  n := obj.catalogue.n
  m := obj.catalogue.m
  // loop over k-mer sizes, smaller k-mers must be added first
  for k := n; k <= m; k++ {
    // first add observed k-mers (i.e. those without any
    // ambiguous characters)
    for i, elements := range obj.catalogue.elements[k-n] {
      kmer := NewKmerClass(k, i, elements)
      if kmer.CountAmbiguous(obj.catalogue.al) == 0 {
        obj.newNode(kmer)
      }
    }
    // add k-mers with ambiguous characters to the graph
    for i, elements := range obj.catalogue.elements[k-n] {
      kmer := NewKmerClass(k, i, elements)
      if kmer.CountAmbiguous(obj.catalogue.al) != 0 {
        // TODO
      }
    }
  }
}

func (obj *KmerGraph) newNode(kmer KmerClass) *KmerGraphNode {
  if node, ok := obj.nodes[kmer.KmerClassId]; ok {
    return node
  }
  r := KmerGraphNode{}
  m := kmer.CountAmbiguous(obj.catalogue.al)
  if m == 0 {
    // this k-mer has no ambiguous positions, fill
    // intra, extra, and infra k-mers
    obj.computeIntraAndExtra(&r, kmer)
    //r.Infra = obj.computeInfra(kmer)
  } else {
    // this k-mer has ambiguous entries, fill
    // infra and supra k-mers
    //obj.computeInfra(kmer)
    //obj.computeSupra(kmer)
  }
  obj.nodes[kmer.KmerClassId] = &r
  return &r
}

/* -------------------------------------------------------------------------- */

func (obj *KmerGraph) computeIntraAndExtra(node1 *KmerGraphNode, kmer KmerClass) {
  if obj.catalogue.n > kmer.K-1 {
    return
  }
  // two possible intra-k-mers
  s1 := kmer.Elements[0][0:kmer.K-1]
  s2 := kmer.Elements[0][1:kmer.K]
  if node2 := obj.GetNode(s1); node2 != nil {
    node2.Extra = append(node2.Extra, node1)
    node1.Intra = append(node1.Intra, node2)
  }
  if node2 := obj.GetNode(s2); node2 != nil {
    node2.Extra = append(node2.Extra, node1)
    node1.Intra = append(node1.Intra, node2)
  }
}
