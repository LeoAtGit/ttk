#pragma once

#include <Triangulation.h>

namespace ttk {

  class MergeTreeRefinement : virtual public Debug {

    public:

      MergeTreeRefinement(){
          this->setDebugMsgPrefix("MergeTreeRefinement");
      };
      ~MergeTreeRefinement(){};
  };
}