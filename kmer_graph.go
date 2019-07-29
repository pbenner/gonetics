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
  clusters  map[int][]KmerGraphNode
}

/* -------------------------------------------------------------------------- */

type KmerGraphNode struct {
  Kmer      KmerClass
  ClusterId uint64
  Rank       int
  // less specific k-mers 
  Infra []*KmerGraphNode
  // more specific k-mers 
  Supra []*KmerGraphNode
  // k-mers of same length with wildcards at equivalent positions
  Intra []*KmerGraphNode
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
    for _, n := range node.Infra {
      r = append(r, n.Kmer)
    }
    for _, n := range node.Supra {
      r = append(r, n.Kmer)
    }
    r.Sort()
    return r
  }
}

/* -------------------------------------------------------------------------- */

// number of non-ambiguous characters
func (obj *KmerGraph) rank(kmer KmerClass) int {
  j := 0
  s := []byte(kmer.Elements[0])
  for i := 0; i < len(s); i++ {
    if ok, err := obj.catalogue.Alphabet.IsWildcard(s[i]); err != nil {
      panic(err)
    } else {
      if !ok {
        j += 1
      }
    }
  }
  return j-1
}

func (obj *KmerGraph) clusterId(kmer KmerClass) uint64 {
  j := uint64(0)
  s := []byte(kmer.Elements[0])
  for i := uint(0); i < uint(len(s)); i++ {
    if ok, err := obj.catalogue.Alphabet.IsWildcard(s[i]); err != nil {
      panic(err)
    } else {
      if !ok {
        j += (1<<i)
      }
    }
  }
  return j
}

func (obj *KmerGraph) constructGraphLoop(n, m int) {
  nodes := make([]map[uint64][]*KmerGraphNode, m+1)
  for i := 0; i < len(nodes); i++ {
    nodes[i] = make(map[uint64][]*KmerGraphNode)
  }
  // sort k-mers by rank
  for k := n; k <= m; k++ {
    for i, elements := range obj.catalogue.elements[k-n] {
      kmer := NewKmerClass(k, i, elements)
      node := obj.newNode(kmer)
      nodes[node.Rank][node.ClusterId] = append(nodes[node.Rank][node.ClusterId], node)
    }
  }
  // remove k-mer lists without any entries
  for j := 0; j < len(nodes); j++ {
    if len(nodes[j]) == 0 {
      nodes = append(nodes[0:j], nodes[j+1:]...)
    }
  }
  // connect nodes
  for j := 1; j < len(nodes); j++ {
    // loop over k-mers with rank j
    for clusterId1, cluster1 := range nodes[j] {
      for clusterId2, cluster2 := range nodes[j-1] {
        if obj.relatedClusters(clusterId1, clusterId2) {
          obj.connectClusters(cluster1, cluster2)
        }
      }
    }
  }
}

func (obj *KmerGraph) constructGraph() {
  n := obj.catalogue.N
  m := obj.catalogue.M
  obj.constructGraphLoop(n, m)
}

func (obj *KmerGraph) newNode(kmer KmerClass) *KmerGraphNode {
  r := KmerGraphNode{Kmer: kmer}
  r.ClusterId = obj.clusterId(kmer)
  r.Rank      = obj.rank     (kmer)
  obj.nodes[kmer.KmerClassId] = &r
  return &r
}

func (obj *KmerGraph) relatedClusters(clusterId1, clusterId2 uint64) bool {
  for k := uint(0); k < uint(obj.catalogue.M+1); k++ {
    a := clusterId1 & (1<<k) != 0
    b := clusterId2 & (1<<k) != 0
    if !(a && b || !b) {
      return false
    }
  }
  return true
}

func (obj *KmerGraph) connectNodes(node1, node2 *KmerGraphNode) {
  if len(node2.Kmer.Elements[0]) <= len(node1.Kmer.Elements[0]) {
    if node2.Kmer.Matches(node1.Kmer, obj.catalogue.Alphabet) {
      node1.Supra = append(node1.Supra, node2)
      node2.Infra = append(node2.Infra, node1)
    }
  }
}

func (obj *KmerGraph) connectClusters(cluster1, cluster2 []*KmerGraphNode) {
  for i := 0; i < len(cluster1); i++ {
    for j := 0; j < len(cluster2); j++ {
      obj.connectNodes(cluster1[i], cluster2[j])
    }
  }
}
