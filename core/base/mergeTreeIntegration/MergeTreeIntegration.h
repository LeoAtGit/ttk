#pragma once

#include <Triangulation.h>

namespace ttk {

  class MergeTreeIntegration : virtual public Debug {

    public:

      MergeTreeIntegration(){
          this->setDebugMsgPrefix("MergeTreeIntegration");
      };
      ~MergeTreeIntegration(){};
  };
}