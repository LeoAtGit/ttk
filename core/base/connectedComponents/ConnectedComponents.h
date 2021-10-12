/// \ingroup base
/// \class ttk::ConnectedComponents
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 01.02.2019
///
/// \brief TTK %connectedComponents processing package.
///
/// This filter consumes a scalar field with a feature mask and computes for
/// each edge connected group of vertices with a non-background mask value a
/// so-called connected component via flood-filling, where the backgroud is
/// masked with values smaller-equal zero. The computed components store the
/// size, seed, and center of mass of each component. The flag
/// UseSeedIdAsComponentId controls if the resulting segmentation is either
/// labeled by the index of the component, or by its seed location (which can be
/// used as a deterministic component label).

#pragma once

// base code includes
#include <Debug.h>
#include <Triangulation.h>

typedef ttk::SimplexId TID;

namespace ttk {
  class ConnectedComponents : virtual public Debug {

  public:
    struct Component {
      int seed = -1;
      float center[3]{0, 0, 0};
      float size{0};
    };

    ConnectedComponents() {
      this->setDebugMsgPrefix("ConnectedComponents");
    };
    ~ConnectedComponents(){};

    int preconditionTriangulation(
      ttk::AbstractTriangulation *triangulation) const {
      return triangulation->preconditionVertexNeighbors();
    };

    template <class DT, class TT = ttk::AbstractTriangulation>
    int computeFloodFill(DT *labels,
                         std::vector<Component> &components,
                         std::vector<TID> &stack,

                         const TT *triangulation,
                         const TID seed) const {
      // get component id
      const TID componentId = components.size();

      constexpr DT featureLabel = -2;

      int stackSize = 1;
      stack[0] = seed;
      labels[seed] = componentId;

      TID cIndex;
      TID nIndex;
      float size = 0;
      float x, y, z;
      float center[3] = {0, 0, 0};

      while(stackSize > 0) {
        cIndex = stack[--stackSize];

        // update node data
        triangulation->getVertexPoint(cIndex, x, y, z);
        center[0] += x;
        center[1] += y;
        center[2] += z;
        size++;

        // add neihbors
        size_t nNeighbors = triangulation->getVertexNeighborNumber(cIndex);
        for(size_t i = 0; i < nNeighbors; i++) {
          triangulation->getVertexNeighbor(cIndex, i, nIndex);
          if(labels[nIndex] == featureLabel) {
            labels[nIndex] = componentId;
            stack[stackSize++] = nIndex;
          }
        }
      }
      center[0] /= size;
      center[1] /= size;
      center[2] /= size;

      // create component
      components.resize(componentId + 1);
      auto &c = components[componentId];
      std::copy(center, center + 3, c.center);
      c.size = size;
      c.seed = seed;

      return 1;
    }

    template <class DT, class TT = ttk::AbstractTriangulation>
    int computeConnectedComponents(std::vector<Component> &components,
                                   int *outputLabels,
                                   const DT *featureMask,
                                   const TT *triangulation,
                                   const bool useSeedAsComponentId
                                   = false) const {

      Timer timer;

      TID nVertices = triangulation->getNumberOfVertices();

      // init stack
      std::vector<TID> stack(nVertices);

      constexpr int backgroundLabel = -1;
      constexpr int featureLabel = -2;

      std::string msg = "Computing Connected Components";
      this->printMsg(msg, 0, 0, 1, ttk::debug::LineMode::REPLACE);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
      for(TID i = 0; i < nVertices; i++)
        outputLabels[i] = featureMask[i] <= 0 ? backgroundLabel : featureLabel;

      for(TID i = 0; i < nVertices; i++)
        if(outputLabels[i] == featureLabel)
          this->computeFloodFill<int, TT>(
            outputLabels, components, stack, triangulation, i);

      this->printMsg(msg, 1, timer.getElapsedTime(), 1);

      if(useSeedAsComponentId) {
        timer.reStart();
        msg = "Labeling Components by Seed Id";
        this->printMsg(
          msg, 0, 0, this->threadNumber_, ttk::debug::LineMode::REPLACE);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
        for(TID i = 0; i < nVertices; i++) {
          auto &cid = outputLabels[i];
          if(cid >= 0)
            cid = components[cid].seed;
        }

        this->printMsg(msg, 1, timer.getElapsedTime(), this->threadNumber_);
      }

      return 1;
    }
  };
} // namespace ttk