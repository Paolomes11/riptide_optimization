#ifndef RIPTIDE_ACTION_INITIALIZATION_HPP
#define RIPTIDE_ACTION_INITIALIZATION_HPP

#include <G4VUserActionInitialization.hh>
#include <string>

namespace riptide {

class ActionInitialization : public G4VUserActionInitialization {
  std::string m_output_file;

 public:
  explicit ActionInitialization(const std::string& output_file);

  void BuildForMaster() const override;
  void Build() const override;
};

} // namespace riptide

#endif // RIPTIDE_ACTION_INITIALIZATION_HPP