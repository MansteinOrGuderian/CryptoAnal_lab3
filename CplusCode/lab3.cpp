#pragma once
#define NOMINMAX // Disable Windows min/max macros that conflict with GMP
#ifdef _WIN32
#include <windows.h>
#endif

#include <iostream>      // For console input/output (std::cout, cin)
#include <fstream>       // For file operations
#include <sstream>       // For string stream operations
#include <string>        // For string class
#include <unordered_map> // For hash table data structure
#include <vector>        // For dynamic arrays
#include <chrono>        // For time measurement
#include <random>        // For random number generation
#include <algorithm>     // For algorithms like std::swap

// GMP library with MSVC warning suppressions
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable: 4146) // unary minus operator applied to unsigned type
#pragma warning(disable: 4244) // conversion from 'mp_limb_t' to 'unsigned long'
#pragma warning(disable: 4267) // conversion from 'size_t' to 'int'
#pragma warning(disable: 4800) // forcing value to bool
#pragma warning(disable: 4018) // signed/unsigned mismatch
#endif

#include <gmpxx.h> // GNU Multiple Precision arithmetic library C++ wrapper

#ifdef _MSC_VER
#pragma warning(pop)
#endif

using namespace std;
using namespace chrono;

// Testing settings
const string MitM_test_path = "../var_data/test/MitM_RSA_256_var15.txt"; // ../var_data/MitM_RSA_2048_20_regular_hard_var15.txt
const uint32_t E_CONST = 65537; // Common public exponent for RSA

const string SE_test_path = "../var_data/test/SE_RSA_256_var15.txt"; // ../var_data/SE_RSA_1024_5_hard_var15.txt
const uint32_t SE_COUNT = 3; // 5

// Parsing input hex to BigInteger
unordered_map<string, mpz_class> read_txt(const string& path) {
    unordered_map<string, mpz_class> test_values;
    ifstream file(path);
    if (!file.is_open())
        throw runtime_error("Cannot open file: " + path);

    string line;
    while (getline(file, line)) {
        size_t pos = line.find(" = ");
        if (pos != string::npos) {
            string key = line.substr(0, pos);
            string val = line.substr(pos + 5); // Skip " = 0x"
            mpz_class value;
            value.set_str(val, 16);
            test_values[key] = value;
        }
    }

    return test_values;
}

int main() {

}
#ifdef _MSC_VER
#pragma warning(pop)
#endif
