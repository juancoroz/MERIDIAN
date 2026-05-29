//  mission_progress.cpp -- shared progress atomics, definition unit
#include "mission_progress.h"

namespace rocket6dof {
namespace progress {

std::atomic<double> current_t  {0.0};
std::atomic<double> total_t    {0.0};
std::atomic<bool>   is_running {false};
std::atomic<int>    last_rc    {0};

}  // namespace progress
}  // namespace rocket6dof
