#include "test_utils.h"

#include <gtest/gtest.h>

// clang-format off
// #include "ExampleTests.h"
#include "FlippingTests.h"
#include "TracingTests.h"
#include "CommonRefinementTests.h"
// clang-format on

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
