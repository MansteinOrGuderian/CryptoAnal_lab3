#pragma once
#define NOMINMAX // Disable Windows min/max macros that conflict with GMP
#ifdef _WIN32
#include <windows.h>
#endif

#include <iostream>      // For console input/output (std::cout, std::cin)
#include <fstream>       // For file operations
#include <sstream>       // For string stream operations
#include <string>        // For string class
#include <unordered_map> // For hash table data structure
#include <vector>        // For dynamic arrays
#include <chrono>        // For time measurement
#include <random>        // For random number generation
#include <algorithm>     // For algorithms like swap

#include <thread>        // For multi-threading
#include <numeric>       // For numeric operations
#include <mutex>         // For thread synchronization
#include <atomic>        // For atomic operations in threading
#include <iomanip>       // For setprecision

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
const string SE_path_to_txt = "../var_data/SE_RSA_1024_5_hard_var15.txt"; // ../var_data/test_SE_RSA_256_var15.txt 
const uint32_t SE_COUNT = 5; // 3 - for test data

// Man(Meet)-in-the-Middle attack
const string MitM_path_to_txt = "../var_data/MitM_RSA_2048_20_regular_hard_var15.txt"; // ../var_data/test_MitM_RSA_256_var15.txt
const uint32_t E_CONST = 65537; // ../var_data/bonus_MitM_RSA_2048_56_var15.txt"
const uint32_t L_CONST = 56;

// For concurrency calculation in MitM
const uint32_t BLOCK_POWER = 24;
const size_t SIZE_OF_BLOCK = 1ULL << BLOCK_POWER;
const size_t USER_BLOCK_SHIFT = 0;

// ---------------------------------------------------------------------------
// ------------------ Parsing input hex to Integers --------------------------
// ---------------------------------------------------------------------------

unordered_map<string, mpz_class> read_txt(const string& path);

// ---------------------------------------------------------------------------
// ------------------------- Small Exponent attack ---------------------------
// ---------------------------------------------------------------------------

// Extended Euclidean Algorithm for Bezout coefficients
std::pair<mpz_class, mpz_class> Extended_Euclidean_algorithm(const mpz_class& a, const mpz_class& b);

// Chinese Remainder Theorem solver. System of x ≡ Cs[i] (mod Ns[i])
mpz_class Chinese_Remainder_Theorem_solver(const std::vector<mpz_class>& Ns, const std::vector<mpz_class>& Cs);

// Custom compute of n-th root using Newton's method (tangent method)
void nth_root(mpz_class& result, const mpz_class& n, unsigned long k);

// Small Exponent attack
duration<double, micro> Small_Exponent_attack();

// ---------------------------------------------------------------------------
// -------------------------- Utility functions ------------------------------
// ---------------------------------------------------------------------------

// (a * x) mod n = 1
mpz_class modinv(const mpz_class& a, const mpz_class& n);

// base^exp mod mod
mpz_class modpow(const mpz_class& base, const mpz_class& exp, const mpz_class& mod);

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

duration<double, micro> Pure_Meet_in_the_Middle_attack();

duration<double, micro> Optimized_Space_Meet_in_the_Middle_attack();

// Result structure for thread-safe result sharing
struct MitMResult {
    std::atomic<bool> found{ false };
    mpz_class M;
    mpz_class T;
    mpz_class S;
    std::mutex result_mutex;
};

// Optimized parallel streaming Meet-in-the-Middle attack. For cases ONLY М^е < N
duration<double, micro> Meet_in_the_Middle_attack_streaming_parallel();

// ---------------------------------------------------------------------------
// ------------------- Bruteforce attack time measure  -----------------------
// ---------------------------------------------------------------------------

duration<double, micro> bruteforce_try_once_time();
