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

func NewKmerGraph(kmers KmerClassList, rel KmerEquivalenceRelation) KmerGraph {
  catalogue := newKmerCatalogue(rel)
  // insert k-mers into the catalogue
  for _, kmer := range kmers {
    catalogue.AddKmerClass(kmer)
  }
  return NewKmerGraphFromCatalogue(catalogue)
}

func NewKmerGraphFromCatalogue(catalogue *KmerCatalogue) KmerGraph {
  r := KmerGraph{}
  r.catalogue = *catalogue
  r.nodes     = make(map[KmerClassId]*KmerGraphNode)
  r.constructGraph()
  return r
}

/* -------------------------------------------------------------------------- */

func (obj KmerGraph) GetNode(kmer string) *KmerGraphNode {
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

func (obj KmerGraph) RelatedKmers(kmer string) KmerClassList {
  if node := obj.GetNode(kmer); node == nil {
    return nil
  } else {
    r := KmerClassList(nil)
    for _, n := range node.Intra {
      r = append(r, n.Kmer)
    }
    {
      r = append(r, obj.relatedKmersSupra(node)...)
      r = append(r, obj.relatedKmersInfra(node)...)
    }
    for _, n := range node.Extra {
      r = append(r, n.Kmer)
    }
    r.Sort()
    return r
  }
}

func (obj KmerGraph) relatedKmersInfra(node *KmerGraphNode) KmerClassList {
  r := KmerClassList(nil)
  for _, n := range node.Infra {
    r = append(r, n.Kmer)
    r = append(r, obj.relatedKmersInfra(n)...)
  }
  return r
}

func (obj KmerGraph) relatedKmersSupra(node *KmerGraphNode) KmerClassList {
  r := KmerClassList(nil)
  for _, n := range node.Supra {
    r = append(r, n.Kmer)
    r = append(r, obj.relatedKmersSupra(n)...)
  }
  return r
}

/* -------------------------------------------------------------------------- */

func (obj *KmerGraph) constructGraphLoop(k, n, m int) {
  kmers := make(  [] KmerClassList, k)
  nodes := make([][]*KmerGraphNode, k)
  // sort k-mers by number of ambiguous characters
  for i, elements := range obj.catalogue.elements[k-n] {
    kmer := NewKmerClass(k, i, elements)
    j    := kmer.CountAmbiguous(obj.catalogue.Alphabet)
    kmers[j] = append(kmers[j], kmer)
  }
  // first add observed k-mers (i.e. those without any
  // ambiguous characters)
  for _, kmer := range kmers[0] {
    nodes[0] = append(nodes[0], obj.newNode1(kmer))
  }
  // add k-mers with ambiguous characters to the graph
  for j := 1; j < k; j++ {
    // loop over k-mers with j ambiguous characters
    for _, kmer := range kmers[j] {
      nodes[j] = append(nodes[j], obj.newNode2(kmer, nodes[j-1]))
    }
  }
}

func (obj *KmerGraph) constructGraph() {
  n := obj.catalogue.N
  m := obj.catalogue.M
  // loop over k-mer sizes, smaller k-mers must be added first
  for k := n; k <= m; k++ {
    obj.constructGraphLoop(k, n, m)
  }
}

func (obj *KmerGraph) newNode1(kmer KmerClass) *KmerGraphNode {
  r := KmerGraphNode{Kmer: kmer}
  obj.computeIntraAndExtra(&r)
  obj.nodes[kmer.KmerClassId] = &r
  return &r
}

func (obj *KmerGraph) newNode2(kmer KmerClass, nodes []*KmerGraphNode) *KmerGraphNode {
  r := KmerGraphNode{Kmer: kmer}
  obj.computeInfraAndSupra(&r, nodes)
  obj.nodes[kmer.KmerClassId] = &r
  return &r
}

/* -------------------------------------------------------------------------- */

func (obj *KmerGraph) computeIntraAndExtra(node1 *KmerGraphNode) {
  if obj.catalogue.N > node1.Kmer.K-1 {
    return
  }
  // two possible intra-k-mers
  s1 := node1.Kmer.Elements[0][0:node1.Kmer.K-1]
  s2 := node1.Kmer.Elements[0][1:node1.Kmer.K]
  if node2 := obj.GetNode(s1); node2 != nil {
    node2.Extra = append(node2.Extra, node1)
    node1.Intra = append(node1.Intra, node2)
  }
  if node2 := obj.GetNode(s2); node2 != nil {
    node2.Extra = append(node2.Extra, node1)
    node1.Intra = append(node1.Intra, node2)
  }
}

func (obj *KmerGraph) computeInfraAndSupra(node1 *KmerGraphNode, nodes []*KmerGraphNode) {
  for _, node2 := range nodes {
    if node1.Kmer.Matches(node2.Kmer, obj.catalogue.Alphabet) {
      node1.Supra = append(node1.Supra, node2)
      node2.Infra = append(node2.Infra, node1)
    }
  }
}
