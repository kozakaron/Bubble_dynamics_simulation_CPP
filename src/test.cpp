#ifdef TEST
#include "test.h"

void testing::Tester::run_tests(){
    size_t passed = 0;
    size_t failed = 0;

    for (size_t i = 0; i < this->test_cases.size(); i++){
        set_up();

        TestCase &test_case = this->test_cases[i];
        std::string test_name = this->test_case_names[i];
        std::string ret = test_case();

        if (ret.empty()){
            cout << BOLD << GREEN << "PASSED: " << RESET << test_name << endl;
            passed++;
        } else {
            cout << BOLD << RED <<   "FAILED: " << RESET << test_name << endl;
            cout << ret << endl;
            failed++;
        }

        tear_down();
    }

    cout << "_____________________________________________" << endl;
    cout << BOLD << GREEN << "PASSED: " << RESET << passed << " out of " << passed + failed << endl;
    cout << BOLD << RED << "FAILED: " << RESET << failed << " out of " << passed + failed << endl;
}

void testing::Tester::add_test(const std::string& test_name, TestCase test_case){
    this->test_case_names.emplace_back(test_name);
    this->test_cases.emplace_back(test_case);
}

#endif // TEST