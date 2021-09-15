/// \ingroup base
/// \class ttk::PerlinVectorField
/// \author Emma Nilsson <emma.nilsson@liu.se>
/// \date 2021-08-31.
///
/// This module defines the %PerlinVectorField class that computes vector field
/// based on Perlin Noise.
///
/// \b Related \b publication: \n
/// 'PerlinVectorField'
/// Jonas Lukasczyk and Julien Tierny.
/// TTK Publications.
/// 2021.
///

#pragma once

// ttk common includes
#include <Debug.h>
#include <Triangulation.h>

namespace ttk {


  class PerlinVectorField : virtual public Debug {

  public:
    PerlinVectorField();

  }; // PerlinVectorField class

} // namespace ttk
