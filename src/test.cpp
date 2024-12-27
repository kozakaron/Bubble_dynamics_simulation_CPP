#ifdef TEST
#include "test.h"

size_t testing::Tester::total_passed = 0;
size_t testing::Tester::total_failed = 0;

testing::Tester::Tester(){
    this->test_group_name = "Unnamed test group";
    this->test_cases = {};
    this->test_case_names = {};
}

testing::Tester::Tester(std::string test_group_name){
    this->test_group_name = test_group_name;
    this->test_cases = {};
    this->test_case_names = {};
}

void testing::Tester::set_up(){
    if (ErrorHandler::get_error_count() != 0) ErrorHandler::print_errors();
    ErrorHandler::clear_errors();
    ErrorHandler::print_when_log = false;
}

void testing::Tester::tear_down(){
    if (ErrorHandler::get_error_count() != 0) ErrorHandler::print_errors();
    ErrorHandler::clear_errors();
    ErrorHandler::print_when_log = true;
}

void testing::Tester::run_tests(){
    size_t passed = 0;
    size_t failed = 0;
    std::cout << BOLD << this->test_group_name << RESET << std::endl;

    for (size_t i = 0; i < this->test_cases.size(); i++){
        set_up();

        TestCase &test_case = this->test_cases[i];
        std::string test_name = this->test_case_names[i];
        std::string ret = test_case();

        if (ret.empty()){
            std::cout << BOLD << GREEN << "    PASSED: " << RESET << test_name << std::endl;
            passed++;
        } else {
            std::cout << BOLD << RED <<   "    FAILED: " << RESET << test_name << std::endl;
            std::cout << "            " << ret << std::endl;
            failed++;
        }

        tear_down();
    }

    testing::Tester::total_passed += passed;
    testing::Tester::total_failed += failed;
}

void testing::Tester::add_test(const std::string& test_name, TestCase test_case){
    this->test_case_names.emplace_back(test_name);
    this->test_cases.emplace_back(test_case);
}

void testing::Tester::print_summary(){
    std::cout << "_____________________________________________" << std::endl;
    std::cout << BOLD << GREEN << "TOTAL PASSED: " << RESET << testing::Tester::total_passed << std::endl;
    std::cout << BOLD << RED << "TOTAL FAILED: " << RESET << testing::Tester::total_failed << std::endl;
}

#endif // TEST