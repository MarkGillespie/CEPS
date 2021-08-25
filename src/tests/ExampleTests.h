#include "test_utils.h"

class ExampleTest : public ::testing::Test {
public:
  static std::unique_ptr<double> globalDouble;
  double localDouble;

protected:
  static void SetUpTestSuite() {
    double *ten = new double();
    *ten = 10;
    globalDouble = std::unique_ptr<double>(ten);
  }

  void SetUp() override { localDouble = 10; }
};

std::unique_ptr<double> ExampleTest::globalDouble = nullptr;

TEST_F(ExampleTest, ExampleTestCase) {
  EXPECT_DOUBLE_EQ(*globalDouble, localDouble);
}
