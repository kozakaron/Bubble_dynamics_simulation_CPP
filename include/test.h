#ifndef TESTS_H
#define TESTS_H
#include <functional>
#include <vector>
#include <string>
#include <iostream>
#include <exception>
#include <cmath>

#include "common.h"
#include "parameters.h"

using std::cout;
using std::endl;
namespace testing{

// Test fixture class
class Tester{
public:

    using TestCase = std::function<std::string()>;
    Tester() {};
    ~Tester() {};
    virtual void set_up() = 0;
    virtual void tear_down() = 0;
    // Run all tests
    void run_tests();
    // Add a test case
    // Usage: add_test("test_name", [&](){ /* test code */ });
    void add_test(const std::string& test_name, TestCase test_case);
private:
    std::vector<std::string> test_case_names;
    std::vector<TestCase> test_cases;
};

#define ADD_TEST(tester, name, ...) \
    tester.add_test(name, [&]() -> std::string { \
        __VA_ARGS__ \
        return ""; \
    });

// Assertion macros
#define FAIL(message) \
    return message;

#define ASSERT_TRUE(condition) \
    if (!(condition)) return "Assert true failed: " #condition " != true";

#define ASSERT_FALSE(condition) \
    if ((condition)) return "Assert false failed: " #condition " != false";

#define ASSERT_EQUAL(a, b) \
    if ((a) != (b)) return "Assert equal failed: " #a " != " #b;

#define ASSERT_NOT_EQUAL(a, b) \
    if ((a) == (b)) return "Assert not equal failed: " #a " == " #b;

#define ASSERT_NEAR(a, b, tol) \
    if (!std::isfinite((a))) return "Assert equal failed: " #a " is not a number"; \
    if (!std::isfinite((b))) return "Assert equal failed: " #b " is not a number"; \
    if (std::abs((a) - (b)) > (tol)) \
        return "Assert near failed: " #a " = " + std::to_string(a) + "   !=   " #b " = " + std::to_string(b);

#define ASSERT_APPROX(a, b, tol) \
    if (!std::isfinite((a))) return "Assert approx failed: " #a " is not a number"; \
    if (!std::isfinite((b))) return "Assert approx failed: " #b " is not a number"; \
    if (std::abs((a) - (b)) / std::max(std::abs(a), std::abs(b)) > (tol)) \
        return "Assert approx failed: " #a " = " + std::to_string(a) + "   !=   " #b " = " + std::to_string(b);

#define ASSERT_APPROX_ARRAY(a, b, size, tol) \
    for (size_t i = 0; i < size; i++) \
    { \
        if (!std::isfinite((a[i]))) return "Assert approx failed: " #a "[" + std::to_string(i) + "] is not a number at "; \
        if (!std::isfinite((b[i]))) return "Assert approx failed: " #b "[" + std::to_string(i) + "] is not a number"; \
        if (std::abs((a[i]) - (b[i])) / std::max(std::abs(a[i]), std::abs(b[i])) > (tol)) \
            return "Assert approx failed: " #a "[" + std::to_string(i) + "] = " + std::to_string(a[i]) + "   !=   " #b "[" + std::to_string(i) + "] = " + std::to_string(b[i]); \
    }

} // namespace testing
#endif // TESTS_H