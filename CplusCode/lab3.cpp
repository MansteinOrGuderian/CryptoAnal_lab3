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

#include <thread>        // For multi-threading
#include <numeric>       // For numeric operations
#include <mutex>         // For thread synchronization
#include <atomic>        // For atomic operations in threading

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

// ---------------------------------------------------------------------------
// ----------------------- Constants & Input Data ----------------------------
// ---------------------------------------------------------------------------

// Small exponent attack
const string SE_test_path = "../var_data/SE_RSA_1024_5_hard_var15.txt"; // ../var_data/test_SE_RSA_256_var15.txt 
const uint32_t SE_COUNT = 5; // 3 - for test data

// Man(Meet)-in-the-Middle attack
const string MitM_test_path = "../var_data/MitM_RSA_2048_20_regular_hard_var15.txt"; // ../var_data/test_MitM_RSA_256_var15.txt
const uint32_t E_CONST = 65537; // ../var_data/bonus_MitM_RSA_2048_56_var15.txt"
const uint32_t L_CONST = 56;

// For concurrency calculation in MitM
const uint32_t BLOCK_POWER = 24;
const size_t BLOCK_SIZE = 1ULL << BLOCK_POWER;
const size_t USER_BLOCK_SHIFT = 0;

// ---------------------------------------------------------------------------
// ------------------ Parsing input hex to Integers --------------------------
// ---------------------------------------------------------------------------

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


// ---------------------------------------------------------------------------
// ------------------------- Small Exponent attack ---------------------------
// ---------------------------------------------------------------------------

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

// ---------------------------------------------------------------------------
// -------------------------- Utility functions ------------------------------
// ---------------------------------------------------------------------------

// (a * x) mod n = 1
mpz_class modinv(const mpz_class& a, const mpz_class& n) {
    mpz_class result;
    mpz_invert(result.get_mpz_t(), a.get_mpz_t(), n.get_mpz_t());
    return result; // x
}

// base^exp mod mod
mpz_class modpow(const mpz_class& base, const mpz_class& exp, const mpz_class& mod) {
    mpz_class result;
    mpz_powm(result.get_mpz_t(), base.get_mpz_t(), exp.get_mpz_t(), mod.get_mpz_t());
    return result;
}

// Custom hash function, allows using GMP big integers as keys in hash maps by converting them to strings first.
struct MpzHash { 
    size_t operator()(const mpz_class& key) const {
        return hash<string>()(key.get_str());
    }
};

// ---------------------------------------------------------------------------
// --------------- Parallel processing function helpers ----------------------
// ---------------------------------------------------------------------------

// A wrapper around unordered_map. Adds thread safety through mutex locking.
template<typename K, typename V, typename H = hash<K>>
class ConcurrentMap { 
private:
    unordered_map<K, V, H> map_;
    mutable mutex mutex_;

public:
    void insert(const K& key, const V& value) {
        lock_guard<mutex> lock(mutex_);
        map_[key] = value;
    }

    bool find(const K& key, V& value) const {
        lock_guard<mutex> lock(mutex_);
        auto it = map_.find(key);
        if (it != map_.end()) {
            value = it->second;
            return true;
        }
        return false;
    }

    bool contains(const K& key) const {
        lock_guard<mutex> lock(mutex_);
        return map_.find(key) != map_.end();
    }

    unordered_map<K, V, H> get_copy() const {
        lock_guard<mutex> lock(mutex_);
        return map_;
    }

    size_t size() const {
        lock_guard<mutex> lock(mutex_);
        return map_.size();
    }

    void reserve(size_t n) {
        lock_guard<mutex> lock(mutex_);
        map_.reserve(n);
    }
};

// Splits a range of work across multiple CPU threads for faster processing. (Simple parallelization) 
template<typename Func>
void parallel_for_simple(size_t start, size_t end, Func f) {
    size_t num_threads = thread::hardware_concurrency();
    if (num_threads == 0) 
        num_threads = 4; // Default to 4 if detection fails

    vector<thread> threads;
    size_t chunk_size = (end - start) / num_threads;

    for (size_t t = 0; t < num_threads; ++t) {
        size_t thread_start = start + t * chunk_size;
        size_t thread_end = (t == num_threads - 1) ? end : thread_start + chunk_size;

        threads.emplace_back([thread_start, thread_end, &f]() {
            for (size_t i = thread_start; i < thread_end; ++i) {
                f(i);
            }
            });
    }

    for (auto& t : threads) {
        t.join();
    }
}

// ---------------------------------------------------------------------------
// ------------------ Man(Meet)-in-the-Middle attack -------------------------
// ---------------------------------------------------------------------------

// TODO

// ---------------------------------------------------------------------------
// ---------------------------- Main function --------------------------------
// ---------------------------------------------------------------------------

int main() {
    try {
        std::cout << "\n--------------------------------------------\n" << std::endl;

        auto SE_time = Small_Exponent_attack_test();
        std::cout << "'Small exponent' execution time: " << SE_time.count() << " mu s" << std::endl;

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
