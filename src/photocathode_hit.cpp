#include "photocathode_hit.hpp"

namespace riptide
{

    /// ThreadLocal è necessario per dichiarare l'allocatore come variabile unica
    G4ThreadLocal G4Allocator<PhotocathodeHit> *PhotocathodeHitAllocator = nullptr;

} // namespace riptide::optics
