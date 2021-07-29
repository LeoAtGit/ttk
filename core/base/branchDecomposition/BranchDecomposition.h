/// \ingroup base
/// \class ttk::BranchDecomposition
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 29.07.2021
///
/// This module assigns to each vertex of a tracking graph a branch id based on
/// a given attribute. First, all birth nodes are assigned a unique branch id,
/// and then the algorithm iterates over every vertex in order of time and then
/// either inherits the branch id of its largest predecessor (but only if the
/// current vertex is also the largest successor of this predecessor), or the
/// vertex gets assinged a new unique branch id.
///
/// \b Related \b publications: \n
/// 'Nested Tracking Graphs'
/// Jonas Lukasczyk, Gunther Weber, Ross Maciejewski, Christoph Garth, and Heike
/// Leitte. Computer Graphics Forum (Special Issue, Proceedings Eurographics /
/// IEEE Symposium on Visualization). Vol. 36. No. 3. 2017.
///
/// 'Dynamic Nested Tracking Graphs'
/// Jonas Lukasczyk, Christoph Garth, Gunther H. Weber, Tim Biedert, Ross
/// Maciejewski, and Heike Leitte. IEEE Transactions on Visualization and
/// Computer Graphics, Vol. 26, No. 1, 2020.

#pragma once

// ttk common includes
#include <Debug.h>
#include <TrackingGraph.h>
#include <Triangulation.h>

namespace ttk {

  class BranchDecomposition : virtual public Debug {

  public:
    BranchDecomposition() {
      this->setDebugMsgPrefix("BranchDecomposition");
    };

    template <typename IT, typename DT>
    int computeBranchDecompositionByAttribute(int *branchId,
                                              ttk::TrackingGraph &trackingGraph,
                                              const IT *time,
                                              const DT *attribute,
                                              const int attributeAssociation) {
      ttk::Timer globalTimer;
      const int nNodes = trackingGraph.inEdges.size();

      std::vector<int> nodesSortedByTime(nNodes);
      // sort all nodes by time in ascending order
      {
        ttk::Timer timer;
        const std::string msg
          = "Sorting Nodes by Time (#" + std::to_string(nNodes) + ")";
        this->printMsg(msg, 0, 0, 1, ttk::debug::LineMode::REPLACE);

        for(int i = 0; i < nNodes; i++)
          nodesSortedByTime[i] = i;
        std::sort(nodesSortedByTime.begin(), nodesSortedByTime.end(),
                  [=](int a, int b) { return time[a] < time[b]; });

        this->printMsg(msg, 1, timer.getElapsedTime(), 1);
      }

      // Sorting Edges by attribute
      {
        ttk::Timer timer;
        const std::string msg = "Sorting Edges by Attribute";
        this->printMsg(
          msg, 0, 0, this->threadNumber_, ttk::debug::LineMode::REPLACE);

        const auto compareAttributeU = [=](const ttk::TrackingGraph::Edge &a,
                                           const ttk::TrackingGraph::Edge &b) {
          return attribute[a.u] > attribute[b.u];
        };
        const auto compareAttributeV = [=](const ttk::TrackingGraph::Edge &a,
                                           const ttk::TrackingGraph::Edge &b) {
          return attribute[a.v] > attribute[b.v];
        };
        const auto compareAttributeE = [=](const ttk::TrackingGraph::Edge &a,
                                           const ttk::TrackingGraph::Edge &b) {
          return attribute[a.e] > attribute[b.e];
        };

        if(attributeAssociation == 0) {
#pragma omp parallel for num_threads(this->threadNumber_)
          for(int i = 0; i < nNodes; i++) {
            auto &outEdges = trackingGraph.outEdges[i];
            auto &inEdges = trackingGraph.inEdges[i];
            if(outEdges.size() > 1)
              std::sort(outEdges.begin(), outEdges.end(), compareAttributeV);
            if(inEdges.size() > 1)
              std::sort(inEdges.begin(), inEdges.end(), compareAttributeU);
          }
        } else {
#pragma omp parallel for num_threads(this->threadNumber_)
          for(int i = 0; i < nNodes; i++) {
            auto &outEdges = trackingGraph.outEdges[i];
            auto &inEdges = trackingGraph.inEdges[i];
            if(outEdges.size() > 1)
              std::sort(outEdges.begin(), outEdges.end(), compareAttributeE);
            if(inEdges.size() > 1)
              std::sort(inEdges.begin(), inEdges.end(), compareAttributeE);
          }
        }

        this->printMsg(msg, 1, timer.getElapsedTime(), this->threadNumber_);
      }

      std::vector<int> maxPrevNode(nNodes);
      std::vector<int> maxNextNode(nNodes);

      // Initializing Branches
      {
        ttk::Timer timer;
        const std::string msg = "Initializing Branches";
        this->printMsg(
          msg, 0, 0, this->threadNumber_, ttk::debug::LineMode::REPLACE);

#pragma omp parallel for num_threads(this->threadNumber_)
        for(int i = 0; i < nNodes; i++) {
          maxPrevNode[i] = trackingGraph.inEdges[i].size() > 0
                             ? trackingGraph.inEdges[i][0].u
                             : -1;
          maxNextNode[i] = trackingGraph.outEdges[i].size() > 0
                             ? trackingGraph.outEdges[i][0].v
                             : -1;

          branchId[i] = trackingGraph.inEdges[i].size() < 1 ? i : -1;
        }

        this->printMsg(msg, 1, timer.getElapsedTime(), this->threadNumber_);
      }

      // Propagating Branchs
      {
        ttk::Timer timer;
        const std::string msg = "Propagating Branches";
        this->printMsg(msg, 0, 0, 1, ttk::debug::LineMode::REPLACE);

        // propagate branch id along graph
        for(int i = 0; i < nNodes; i++) {
          const auto &v = nodesSortedByTime[i];
          if(branchId[v] != -1)
            continue;

          int maxPrevNodeOfV = -1;
          for(const auto &e : trackingGraph.inEdges[v]) {
            if(maxNextNode[e.u] == v) {
              maxPrevNodeOfV = e.u;
              break;
            }
          }

          branchId[v] = maxPrevNodeOfV < 0 ? v : branchId[maxPrevNodeOfV];
        }

        this->printMsg(msg, 1, timer.getElapsedTime(), this->threadNumber_);
      }

      this->printMsg(ttk::debug::Separator::L2);
      this->printMsg("Complete", 1, globalTimer.getElapsedTime());
      this->printMsg(ttk::debug::Separator::L1);

      return 1;
    }
  }; // BranchDecomposition class

} // namespace ttk
