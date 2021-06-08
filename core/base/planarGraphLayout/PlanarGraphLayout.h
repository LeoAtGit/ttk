/// \ingroup base
/// \typename ttk::PlanarGraphLayout
/// \author Jonas Lukasczyk (jl@jluk.de)
/// \date 01.12.2018
///
/// \brief TTK %planarGraphLayout processing package.
///
/// %PlanarGraphLayout is a TTK processing package that computes a planar graph
/// layout of a \b vtkUnstructuredGrid. To improve the quality of the layout it
/// is possible to pass additional field data to the algorithm:\n \b 1) \b
/// Sequences: Points are positioned along the x-axis based on a sequence (e.g.,
/// time indicies or scalar values). \b 1) \b Sizes: Points cover space on the
/// y-axis based on their size. \b 1) \b Branches: Points with the same branch
/// label are positioned on straight lines. \b 1) \b Levels: The layout of
/// points with the same level label are computed individually and afterwards
/// nested based on the level hierarchy. This makes it possible to draw nested
/// graphs where each level is a layer of the resulting graph.
///
/// \b Related \b publication: \n
/// 'Nested Tracking Graphs'
/// Jonas Lukasczyk, Gunther Weber, Ross Maciejewski, Christoph Garth, and Heike
/// Leitte. Computer Graphics Forum (Special Issue, Proceedings Eurographics /
/// IEEE Symposium on Visualization). Vol. 36. No. 3. 2017.
///

#pragma once

#include <map>
#include <limits>
#include <set>
#include <iterator>
// base code includes
#include <Debug.h>

namespace ttk {

  template <typename DT, typename IT>
  struct Branch {
    DT leaf, root;
    std::map<ttk::LongSimplexId,IT> branchPoints;
    std::set<ttk::LongSimplexId> vertices;
    Branch() {
      leaf = -1;
      root = -1; 
    }

    DT GetMaxValues(const DT &sequenceList) {
      return sequenceList[leaf]; 
    }
    
  };
  
  class PlanarGraphLayout : virtual public Debug {

  public:
    PlanarGraphLayout();
    ~PlanarGraphLayout();

    enum class ALGORITHM {
      DOT = 0,
      MERGE_TREE_OPTIMIZATION = 1
    };

    template <typename DT, typename IT>
    int computeGraphLayout(
      // Output
      float *layout,

      // Input
      const LongSimplexId *connectivityList,
      const size_t &nPoints,
      const size_t &nEdges,
      const DT *pointSequences=nullptr,
      const float *sizes=nullptr,
      const IT *branches=nullptr,
      const IT *levels=nullptr
    ) const;

    template <typename DT, typename IT>
    int GenerateLayout(IT b,
		       float* layout,
		       const std::vector<Branch<DT,IT>> branchList,
		       const DT* pointSequences) const;
    
    template <typename DT, typename IT>
    int computeMergeTreeLayout(
      // Output
      float *layout,

      // Input
      const LongSimplexId *connectivityList,
      const size_t &nPoints,
      const size_t &nEdges,
      const DT *pointSequences,
      const IT *branches
    ) const;

    template <typename IT>
    int extractLevel(
      // Output
      std::vector<size_t> &nodeIndicies,
      std::vector<size_t> &edgeIndicies,

      // Input
      const LongSimplexId *connectivityList,
      const size_t &nPoints,
      const size_t &nEdges,
      const IT &level,
      const IT *levels) const;

    template <typename IT, typename DT>
    int computeDotString(
      // Output
      std::string &dotString,

      // Input
      const LongSimplexId *connectivityList,
      const DT *pointSequences,
      const float *sizes,
      const IT *branches,
      const std::vector<size_t> &nodeIndicies,
      const std::vector<size_t> &edgeIndicies,
      const std::map<DT, size_t> &sequenceValueToIndexMap) const;

    template <typename IT>
    int computeSlots(
      // Output
      float *layout,

      // Input
      const LongSimplexId *connectivityList,
      const size_t &nPoints,
      const size_t &nEdges,
      const float *sizes,
      const IT *levels,
      const IT &nLevels) const;

    // Compute Dot Layout
    int computeDotLayout(
      // Output
      float *layout,

      // Input
      const std::vector<size_t> &nodeIndicies,
      const std::string &dotString) const;


    
  };
} // namespace ttk

// =============================================================================
// Extract Level
// =============================================================================
template <typename IT>
int ttk::PlanarGraphLayout::extractLevel(
  // Output
  std::vector<size_t> &nodeIndicies,
  std::vector<size_t> &edgeIndicies,

  // Input
  const LongSimplexId *connectivityList,
  const size_t &nPoints,
  const size_t &nEdges,
  const IT &level,
  const IT *levels) const {

  // If levels==nullptr then return all points and edges
  if(levels == nullptr) {
    nodeIndicies.resize(nPoints);
    for(size_t i = 0; i < nPoints; i++)
      nodeIndicies[i] = i;

    edgeIndicies.resize(nEdges);
    for(size_t i = 0; i < nEdges; i++)
      edgeIndicies[i] = i;

    return 1;
  }

  // Get nodes at level
  for(size_t i = 0; i < nPoints; i++)
    if(levels[i] == level)
      nodeIndicies.push_back(i);

  // Get edges at level
  size_t nEdges3 = nEdges * 3;
  for(size_t i = 0; i < nEdges3; i += 3) {
    auto n0l = levels[connectivityList[i + 1]];
    auto n1l = levels[connectivityList[i + 2]];
    if(n0l == level && n0l == n1l)
      edgeIndicies.push_back(i / 3);
  }

  return 1;
}

// =============================================================================
// Compute Dot String
// =============================================================================
template <typename IT, typename DT>
int ttk::PlanarGraphLayout::computeDotString(
  // Output
  std::string &dotString,

  // Input
  const LongSimplexId *connectivityList,
  const DT *pointSequences,
  const float *sizes,
  const IT *branches,
  const std::vector<size_t> &nodeIndicies,
  const std::vector<size_t> &edgeIndicies,
  const std::map<DT, size_t> &sequenceValueToIndexMap) const {

  Timer t;

  this->printMsg("Generating DOT String", 0, debug::LineMode::REPLACE);

  bool useSequences = pointSequences != nullptr;
  bool useSizes = sizes != nullptr;
  bool useBranches = branches != nullptr;

  std::string headString = "digraph g {rankdir=LR;";
  std::string nodeString = "";
  std::string edgeString = "";
  std::string rankString = "";

  // lambda functions that generate string representations of nodes
  auto sl = [](size_t s) { return "\"s" + std::to_string(s) + "\""; };
  auto nl = [](size_t id) { return std::to_string(id); };

  // ---------------------------------------------------------------------------
  // Nodes
  // ---------------------------------------------------------------------------
  {
    // Set default node style
    nodeString += "node[label=\"\",shape=box,width=1,height=1];";

    // If useSizes then map size to node height
    if(useSizes)
      for(auto &i : nodeIndicies)
        nodeString += nl(i) + "[height=" + std::to_string(sizes[i]) + "];";
  }

  // ---------------------------------------------------------------------------
  // Ranks
  // ---------------------------------------------------------------------------
  if(useSequences) {
    size_t nSequenceValues = sequenceValueToIndexMap.size();

    // Sequence Chain
    {
      edgeString += sl(0);
      for(size_t s = 1; s < nSequenceValues; s++)
        edgeString += "->" + sl(s);
      edgeString += "[weight=1];";
    }

    // Collect nodes with the same sequence index
    std::vector<std::vector<size_t>> sequenceIndexToPointIndexMap(
      nSequenceValues);
    for(auto &i : nodeIndicies)
      sequenceIndexToPointIndexMap
        [sequenceValueToIndexMap.find(pointSequences[i])->second]
          .push_back(i);

    // Compute individual ranks
    for(size_t s = 0; s < nSequenceValues; s++) {
      rankString += "{rank=same " + sl(s);

      auto &nodes = sequenceIndexToPointIndexMap[s];
      for(auto &i : nodes)
        rankString += " " + nl(i);

      rankString += "}";
    }
  }

  // ---------------------------------------------------------------------------
  // Edges
  // ---------------------------------------------------------------------------
  {
    for(auto &edgeIndex : edgeIndicies) {
      size_t temp = edgeIndex * 3;
      auto &i0 = connectivityList[temp + 1];
      auto &i1 = connectivityList[temp + 2];
      edgeString += nl(i0) + "->" + nl(i1);

      if(useBranches) {
        auto b0 = branches[i0];
        auto b1 = branches[i1];
        edgeString += b0 == b1 ? "[weight=1]" : "[weight=0]";
      }

      edgeString += ";";
    }
  }

  // ---------------------------------------------------------------------------
  // Finalize
  // ---------------------------------------------------------------------------

  // Build Dot String
  { dotString = headString + nodeString + edgeString + rankString + "}"; }

  // Print Status
  this->printMsg("Generating DOT String", 1, t.getElapsedTime());
  this->printMsg("\n" + dotString + "\n", debug::Priority::VERBOSE);

  return 1;
}

// =============================================================================
// Compute Slots
// =============================================================================
template <typename IT>
int ttk::PlanarGraphLayout::computeSlots(
  // Output
  float *layout,

  // Input
  const LongSimplexId *connectivityList,
  const size_t &nPoints,
  const size_t &nEdges,
  const float *sizes,
  const IT *levels,
  const IT &nLevels) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(sizes == nullptr || levels == nullptr) {
    return -1;
  }
#endif // TTK_ENABLE_KAMIKAZE

  Timer t;
  this->printMsg("Computing Slots", 0, debug::LineMode::REPLACE);

  // Comparator that sorts children based on layout.y
  struct ChildrenComparator {
    const float *layout_;

    ChildrenComparator(const float *layout) : layout_(layout){};

    inline bool operator()(const size_t &i, const size_t &j) {
      return layout_[i * 2 + 1] < layout_[j * 2 + 1];
    }
  };

  auto comparator = ChildrenComparator(layout);

  // ---------------------------------------------------------------------------
  // Compute Children
  // ---------------------------------------------------------------------------
  std::vector<std::vector<size_t>> nodeIndexChildrenIndexMap(nPoints);

  size_t nEdges3 = nEdges * 3;
  for(size_t i = 0; i < nEdges3; i += 3) {
    auto n0 = connectivityList[i + 1];
    auto n1 = connectivityList[i + 2];
    if((levels[n0] + 1) == levels[n1])
      nodeIndexChildrenIndexMap[n0].push_back(n1);
  }

  // ---------------------------------------------------------------------------
  // Adjust positions from bottom to top (skip last level)
  // ---------------------------------------------------------------------------
  for(IT l = 0; l < nLevels - 1; l++) {
    std::vector<size_t> nodeIndicies;
    std::vector<size_t> edgeIndicies;

    // get nodes at current level (parents)
    this->extractLevel<IT>(
      // Output
      nodeIndicies, edgeIndicies,

      // Input
      connectivityList, nPoints, nEdges, l, levels);

    // for each parent adjust position of children
    for(auto &parent : nodeIndicies) {
      auto &children = nodeIndexChildrenIndexMap[parent];
      size_t nChildren = children.size();
      if(nChildren < 1)
        continue;

      // sort children
      sort(children.begin(), children.end(), comparator);

      // size of parent
      float sizeParent = sizes[parent];

      // size of child
      float sizeChildren = 0;
      for(auto &child : children)
        sizeChildren += sizes[child];

      // gap space
      float gap = sizeParent - sizeChildren;
      float gapDelta = (gap / (nChildren + 1)) / 2;

      float y = layout[parent * 2 + 1] + sizeParent * 0.5 - gapDelta;
      for(auto &child : children) {
        float temp = gapDelta + sizes[child] / 2;
        layout[child * 2 + 1] = y - temp;
        y -= 2 * temp;
      }
    }
  }

  this->printMsg("Computing Slots", 1, t.getElapsedTime());

  return 1;
}

// =============================================================================
// Compute Layout
// =============================================================================
template <typename DT, typename IT>
int ttk::PlanarGraphLayout::computeGraphLayout(
  // Output
  float *layout,

  // Input
  const LongSimplexId *connectivityList,
  const size_t &nPoints,
  const size_t &nEdges,
  const DT *pointSequences,
  const float *sizes,
  const IT *branches,
  const IT *levels
) const {

  Timer t;

  // Init Input
  bool useSequences = pointSequences != nullptr;
  bool useSizes = sizes != nullptr;
  bool useBranches = branches != nullptr;
  bool useLevels = levels != nullptr;

  // Print Input
  {
    std::string modeS = "";
    if(useSequences)
      modeS += "Sequence + ";
    if(useSizes)
      modeS += "Size + ";
    if(useBranches)
      modeS += "Branches + ";
    if(useLevels)
      modeS += "Levels + ";

    this->printMsg(debug::Separator::L1);
    this->printMsg({{"#Nodes", std::to_string(nPoints)},
                    {"#Edges", std::to_string(nEdges)},
                    {"Mode", modeS.substr(0, modeS.length() - 3)}});
    this->printMsg(debug::Separator::L2);
  }


  if(useLevels && !useSizes) {
    this->printErr("'UseLevels' requires 'UseSizes'.");
    return 0;
  }

  // Global SequenceValue to SequenceIndex map
  std::map<DT, size_t> sequenceValueToIndexMap;
  if(useSequences) {
    for(size_t i = 0; i < nPoints; i++)
      sequenceValueToIndexMap[pointSequences[i]] = 0;
    size_t i = 0;
    for(auto &el : sequenceValueToIndexMap)
      el.second = i++;
  }

  // Get number of levels
  IT nLevels = 1;
  if(useLevels) {
    for(size_t i = 0; i < nPoints; i++)
      if(nLevels < levels[i])
        nLevels = levels[i];
    nLevels += 1;
  }

  // ---------------------------------------------------------------------------
  // Compute initial layout for each level
  // ---------------------------------------------------------------------------
  for(IT l = 0; l < nLevels; l++) {
    std::vector<size_t> nodeIndicies;
    std::vector<size_t> edgeIndicies;

    // Extract nodes and edges at certain level
    {
      int status = this->extractLevel<IT>(
        // Output
        nodeIndicies, edgeIndicies,

        // Input
        connectivityList, nPoints, nEdges, l, levels);
      if(status != 1)
        return 0;
    }

    // Compute Dot String
    std::string dotString;
    {
      int status = this->computeDotString<IT, DT>(
        // Output
        dotString,

        // Input
        connectivityList, pointSequences, sizes, branches, nodeIndicies,
        edgeIndicies, sequenceValueToIndexMap);
      if(status != 1)
        return 0;
    }

    // Compute Dot Layout
    {
      int status = this->computeDotLayout(layout, nodeIndicies, dotString);
      if(status != 1)
        return 0;
    }
  }

  // ---------------------------------------------------------------------------
  // If nLevels>1 then compute slots
  // ---------------------------------------------------------------------------
  if(nLevels > 1) {
    this->computeSlots<IT>(
      // Output
      layout,

      // Input
      connectivityList, nPoints, nEdges, sizes, levels, nLevels);
  }

  // ---------------------------------------------------------------------------
  // Print performance
  // ---------------------------------------------------------------------------
  this->printMsg(debug::Separator::L2);
  this->printMsg("Complete", 1, t.getElapsedTime());
  this->printMsg(debug::Separator::L1);

  return 1;
}

template <typename DT, typename IT>
int ttk::PlanarGraphLayout::GenerateLayout(IT b,
					   float* layout,
					   const std::vector<Branch<DT,IT>> branchList,
					   const DT* pointSequences) const {
  
  Branch<DT,IT> curr_branch = branchList[b];
  std::cout << "Branch " << b << std::endl;
  std::cout << "Leaf " << curr_branch.leaf << std::endl;
  std::cout << "Root " << curr_branch.root << std::endl;

  std::cout << "Number of vertices " << curr_branch.vertices.size() << std::endl; 
  
  for(ttk::LongSimplexId v : curr_branch.vertices) {
    layout[v * 2] = (float) b;
    layout[v * 2 + 1] = (float) pointSequences[v];
  }
  
  // go through the branches, generate a layout for them
  for(std::pair<ttk::LongSimplexId, IT> bp : curr_branch.branchPoints) {
    GenerateLayout(bp.second, layout, branchList, pointSequences);
  }
  
  return 1;
}

template <typename DT, typename IT>
int ttk::PlanarGraphLayout::computeMergeTreeLayout(
  // Output
  float *layout,

  // Input
  const ttk::LongSimplexId *connectivityList,
  const size_t &nPoints,
  const size_t &nEdges,
  const DT *pointSequences,
  const IT *branches
) const {

  Timer timer;

  const std::string msg = "Computing Merge Tree Layout";

  this->printMsg(msg, 0,0,this->threadNumber_,debug::LineMode::REPLACE);

  //get number of unique branches
  std::set<IT> branchNumbers;
  
  for(size_t i = 0; i < nPoints; i++) {
    branchNumbers.insert(branches[i]);
  }  
  
  std::vector<Branch<DT,IT>> branchList(branchNumbers.size());
  
  for(size_t i = 0; i < nEdges; i++) {
    //i*3 gives you number of vertices, always 2 in this case
    ttk::LongSimplexId v1 = connectivityList[i*3+1];
    ttk::LongSimplexId v2 = connectivityList[i*3+2];
    
    // scalar values of vertices
    DT s1 = pointSequences[v1];
    DT s2 = pointSequences[v2];

    // branchId of vertices
    IT b1 = branches[v1];
    IT b2 = branches[v2];
    
    if(b1 == b2) {
      Branch<DT,IT>* b = &branchList[b1];
      
      if(b->leaf != -1) {
	b->leaf = (s1 > pointSequences[(int) b->leaf]) ? v1 : b->leaf;
	b->leaf = (s2 > pointSequences[(int) b->leaf]) ? v2 : b->leaf;
      } else {
	b->leaf = (s1 > s2) ? v1 : v2;
      }
      
      if(b->root != -1) {
	b->root = (s1 < pointSequences[(int) b->root]) ? b->root : v1;
	b->root = (s2 < pointSequences[(int) b->root]) ? b->root : v2;
      } else {
	b->root = (s2 > s1) ? v2 : v1;
      }
      
      b->vertices.insert(v1);
      b->vertices.insert(v2);
      
    } else {
      //we have a branch point
      Branch<DT,IT>* branch1 = &branchList[b1];
      Branch<DT,IT>* branch2 = &branchList[b2];

      //branch from right to left
      if(s2 > s1) {
	branch1->branchPoints.insert(std::pair<ttk::LongSimplexId, IT>(v1,b2));
      } else {
	branch2->branchPoints.insert(std::pair<ttk::LongSimplexId, IT>(v2,b1));
      }

      branch1->vertices.insert(v1);
      branch2->vertices.insert(v2);
    }		 
  }

  GenerateLayout((IT) 0, layout, branchList, pointSequences);
  
  this->printMsg(msg, 1, timer.getElapsedTime());

  return 1;
}


