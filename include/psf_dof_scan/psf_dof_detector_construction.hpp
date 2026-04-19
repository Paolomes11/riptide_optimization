#ifndef RIPTIDE_PSF_DOF_DETECTOR_CONSTRUCTION_HPP
#define RIPTIDE_PSF_DOF_DETECTOR_CONSTRUCTION_HPP

#include "optimization/detector_construction.hpp"

#include <filesystem>
#include <string>

namespace riptide {

class PsfDofDetectorConstruction : public DetectorConstruction {
 public:
  using DetectorConstruction::DetectorConstruction;

  void ConstructSDandField() override {
  }
};

} // namespace riptide

#endif // RIPTIDE_PSF_DOF_DETECTOR_CONSTRUCTION_HPP
