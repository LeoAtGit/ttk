/// \ingroup base
/// \class ttk::ConnectedComponents
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 01.02.2019
///
/// \brief TTK %connectedComponents processing package.
///
/// This filter consumes a scalar field containing feature labels and computes
/// for each edge connected group of vertices with the same label a so-called
/// component, where negative labels represent the background. The computed
/// components store the size and center of mass of each component. The module
/// also assigns a unqiue positive integer to each component and maps this id to
/// the segmentation.

#pragma once

// base code includes
#include <Debug.h>
#include <Triangulation.h>

typedef ttk::SimplexId TID;

namespace ttk {
  class ConnectedComponents : virtual public Debug {

  public:
    struct Component {
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

      return 1;
    }

    template <class DT, class TT = ttk::AbstractTriangulation>
    int computeConnectedComponents(std::vector<Component> &components,
                                   DT *outputLabels,
                                   const DT *inputLabels,
                                   const TT *triangulation) const {

      Timer timer;

      TID nVertices = triangulation->getNumberOfVertices();

      // init stack
      std::vector<TID> stack(nVertices);

      constexpr DT backgroundLabel = -1;
      constexpr DT featureLabel = -2;

      this->printMsg("Computing Connected Components", 0,
                     timer.getElapsedTime(), 1, ttk::debug::LineMode::REPLACE);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
      for(TID i = 0; i < nVertices; i++)
        outputLabels[i] = inputLabels[i] < 0 ? backgroundLabel : featureLabel;

      const TID progressF = nVertices / 10;
      TID progressC = 0;
      double progressT = 0;

      for(TID i = 0; i < nVertices; i++) {

        if(progressC++ > progressF) {
          progressC = 0;
          progressT += 0.1;
          this->printMsg("Computing Connected Components", progressT,
                         timer.getElapsedTime(), 1,
                         ttk::debug::LineMode::REPLACE);
        }

        if(outputLabels[i] == featureLabel) {
          this->computeFloodFill(outputLabels, components, stack,

                                 triangulation, i);
        }
      }

      this->printMsg(
        "Computing Connected Components", 1, timer.getElapsedTime(), 1);

      return 1;
    }
  };
} // namespace ttk