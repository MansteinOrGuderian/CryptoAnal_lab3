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

const string SE_test_path = "../var_data/SE_RSA_1024_5_hard_var15.txt"; // ../var_data/test/SE_RSA_256_var15.txt
const uint32_t SE_COUNT = 5; // 3

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

// Extended Euclidean Algorithm for Bezout coefficients
std::pair<mpz_class, mpz_class> Extended_Euclidean_algorithm(const mpz_class& a, const mpz_class& b) {
    bool swap = false;
    if (a < b)
        swap = true;

    mpz_class m = a, n = b;
    if (swap)
        std::swap(m, n);

    mpz_class old_v = 1, v = 0;
    mpz_class old_u = 0, u = 1;

    while (n != 0) {
        mpz_class q = m / n;
        mpz_class temp = n;
        n = m - (q * n);
        m = temp;

        temp = v;
        v = old_v - (q * v);
        old_v = temp;

        temp = u;
        u = old_u - (q * u);
        old_u = temp;
    }

    if (swap)
        return make_pair(old_u, old_v);

    return make_pair(old_v, old_u);
}

// Chinese Remainder Theorem solver. System of x ≡ Cs[i] (mod Ns[i])
mpz_class Chinese_Remainder_Theorem_solver(const std::vector<mpz_class>& Ns, const std::vector<mpz_class>& Cs) {
    mpz_class PROD_N = 1;
    for (const auto& n : Ns)
        PROD_N *= n;

    std::vector<mpz_class> Ni, Mi;
    for (const auto& ni : Ns) {
        Ni.push_back(PROD_N / ni);
        Mi.push_back(Extended_Euclidean_algorithm(Ni.back(), ni).first);
    }

    mpz_class res = 0;
    for (size_t i = 0; i < Cs.size(); ++i)
        res = (res + (Cs[i] * Ni[i] * Mi[i])) % PROD_N;

    if (res < 0)
        res = res + PROD_N;

    return res;
}

// Small Exponent attack
duration<double, micro> Small_Exponent_attack_test() {
    std::cout << "Getting values from: '" << SE_test_path << "'" << std::endl;
    auto test_values = read_txt(SE_test_path);

    std::vector<mpz_class> Ns, Cs;

    std::cout << "Small Exponent attack started" << std::endl;
    for (uint32_t i = 1; i <= SE_COUNT; ++i) {
        string n_key = "N" + to_string(i);
        string c_key = "C" + to_string(i);

        Ns.push_back(test_values[n_key]);
        Cs.push_back(test_values[c_key]);

        std::cout << n_key << " = " << Ns.back() << std::endl;
        std::cout << c_key << " = " << Cs.back() << std::endl;
    }

    auto timer = steady_clock::now();

    mpz_class C = Chinese_Remainder_Theorem_solver(Ns, Cs);

    uint32_t e = SE_COUNT; // Small exponent is public info
    mpz_class M;
    mpz_root(M.get_mpz_t(), C.get_mpz_t(), e);

    std::cout << "SE message: " << M << '\n' << "SE message in hex : ";
    gmp_printf("%Zx\n", M); // output in hex
    return steady_clock::now() - timer;
}

int main() {
    try {
        std::cout << "\n--------------------------------------------\n" << std::endl;

        auto SE_time = Small_Exponent_attack_test();
        std::cout << "'Small exponent' execution time: " << SE_time.count() << " x 10^-6 s" << std::endl;

    }
    catch (const exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}

#ifdef _MSC_VER
#pragma warning(pop)
#endif
