///----------------------------------------
///     @file memo2_test_stub.cpp
///   @ingroup ASTAP++/tests
///    @brief Test-only no-op definition of astap::core::memo2_message.
///  @details The canonical implementation lives in core/platform.cpp, which
///           installs the process-wide progress sink and pulls in a large
///           slice of the library. Unit tests that exercise only numeric
///           routines link this shim instead, satisfying the linker without
///           that dependency. Progress output is irrelevant to those tests.
///----------------------------------------

#include <string_view>

namespace astap::core {

void memo2_message(std::string_view) {}

}  // namespace astap::core
