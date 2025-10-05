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
duration<double, micro> Small_Exponent_attack() {
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

duration<double, micro> Pure_Meet_in_the_Middle_attack() {
    mpz_class e(static_cast<unsigned long>(E_CONST));

    std::cout << "Getting values from: '" << MitM_test_path << "'" << std::endl;
    auto test_values = read_txt(MitM_test_path);
    mpz_class N = test_values["N"];
    mpz_class C = test_values["C"];

    std::cout << "> Simple MitM attack started for:" << std::endl;
    std::cout << "l = " << L_CONST << std::endl;
    std::cout << "N: " << N << std::endl;
    std::cout << "C: " << C << std::endl;

    auto timer = steady_clock::now();

    size_t size = 1ULL << (L_CONST / 2);

    std::cout << "> MitM: Started pushing at " <<
        duration_cast<microseconds>(steady_clock::now() - timer).count() << std::endl;

    unordered_map<mpz_class, mpz_class, MpzHash> X;
    X.reserve(size);

    // Parallel generation of X table
    vector<pair<mpz_class, mpz_class>> temp_pairs(size);

    parallel_for_simple(1, size + 1, [&](size_t a) {
        mpz_class num(static_cast<unsigned long>(a));
        mpz_class num_e = modpow(num, e, N);
        temp_pairs[a - 1] = make_pair(num_e, num);
        });

    for (const auto& p : temp_pairs) {
        X[p.first] = p.second;
    }

    std::cout << "> MitM: Pushing finished at " <<
        duration_cast<microseconds>(steady_clock::now() - timer).count() << std::endl;

    // Find S
    atomic<bool> found(false);
    mpz_class M;
    mutex result_mutex;

    // Create a vector of keys to iterate over
    vector<mpz_class> keys;
    keys.reserve(X.size());
    for (const auto& item : X) {
        keys.push_back(item.first);
    }

    parallel_for_simple(0, keys.size(), [&](size_t i) {
        if (found.load()) 
            return;

        mpz_class C_S = (C * modinv(keys[i], N)) % N;
        auto find_it = X.find(C_S);

        if (find_it != X.end()) {
            lock_guard<mutex> lock(result_mutex);
            if (!found.load()) {
                found.store(true);
                M = X[keys[i]] * find_it->second;
                std::cout << "> MitM message: " << M << '\n' << "MitM message in hex : ";
                gmp_printf("%Zx\n", M); // output in hex

                // Save to file
                ofstream result_file("Simple_MitM_result.txt");
                result_file << "MitM found message: " << M << "\n";
                result_file << "T: " << find_it->second << "\n";
                result_file << "S: " << X[keys[i]] << std::endl;
            }
        }
        });

    if (!found) {
        std::cout << "> Message for MitM not found! Damn it!" << std::endl;
    }

    return steady_clock::now() - timer;
}

duration<double, micro> Optimized_Space_Meet_in_the_Middle_attack() {
    mpz_class e(static_cast<unsigned long>(E_CONST));

    std::cout << "Getting values from: '" << MitM_test_path << "'" << std::endl;
    auto test_values = read_txt(MitM_test_path);
    mpz_class N = test_values["N"];
    mpz_class C = test_values["C"];

    std::cout << "> MitM (optimized) attack started for:" << std::endl;
    std::cout << "l = " << L_CONST << std::endl;
    std::cout << "N: " << N << std::endl;
    std::cout << "C: " << C << std::endl;

    auto timer = steady_clock::now();
    size_t blocks = 1ULL << ((L_CONST / 2) - BLOCK_POWER);

    for (size_t bn_t = USER_BLOCK_SHIFT; bn_t < blocks; ++bn_t) {
        std::cout << "Start: bn_t = " << bn_t << ", bn_s = " << bn_t << " : " <<
            duration_cast<microseconds>(steady_clock::now() - timer).count() << std::endl;

        // Symmetrical case
        size_t shift_t_start = 1 + bn_t * BLOCK_SIZE;
        size_t shift_t_end = (bn_t + 1) * BLOCK_SIZE;

        ConcurrentMap<mpz_class, mpz_class, MpzHash> T_block;
        T_block.reserve(BLOCK_SIZE);

        parallel_for_simple(shift_t_start, shift_t_end + 1, [&](size_t a) {
            mpz_class num(static_cast<unsigned long>(a));
            mpz_class num_e = modpow(num, e, N);
            T_block.insert(num_e, num);
            });

        std::cout << "Generated: bn_t = " << bn_t << ", bn_s = " << bn_t << " : " <<
            duration_cast<microseconds>(steady_clock::now() - timer).count() << std::endl;

        // Check for solution in symmetric case
        auto T_map_copy = T_block.get_copy();
        for (auto it = T_map_copy.begin(); it != T_map_copy.end(); ++it) {
            mpz_class C_S = (C * modinv(it->first, N)) % N;
            mpz_class found_val;
            if (T_block.find(C_S, found_val)) {
                mpz_class M = it->second * found_val;
                std::cout << "MitM message: " << M << std::endl;

                ofstream result_file("result.txt");
                result_file << "MitM found message: " << M << std::endl;
                return steady_clock::now() - timer;
            }
        }

        std::cout << "End: bn_t = " << bn_t << ", bn_s = " << bn_t << " : " <<
            duration_cast<microseconds>(steady_clock::now() - timer).count() << std::endl;

        // Asymmetrical case
        for (size_t bn_s = bn_t + 1; bn_s < blocks; ++bn_s) {
            std::cout << "Start: bn_t = " << bn_t << ", bn_s = " << bn_s << " : " <<
                duration_cast<microseconds>(steady_clock::now() - timer).count() << std::endl;

            size_t shift_s_start = 1 + bn_s * BLOCK_SIZE;
            size_t shift_s_end = (bn_s + 1) * BLOCK_SIZE;

            ConcurrentMap<mpz_class, mpz_class, MpzHash> S_block;
            S_block.reserve(BLOCK_SIZE);

            parallel_for_simple(shift_s_start, shift_s_end + 1, [&](size_t a) {
                mpz_class num(static_cast<unsigned long>(a));
                mpz_class num_e = modpow(num, e, N);
                S_block.insert(num_e, num);
                });

            // Compare blocks
            auto S_map_copy = S_block.get_copy();
            for (auto it = S_map_copy.begin(); it != S_map_copy.end(); ++it) {
                mpz_class C_S = (C * modinv(it->first, N)) % N;
                mpz_class found_val;
                if (T_block.find(C_S, found_val)) {
                    mpz_class M = it->second * found_val;
                    std::cout << "MitM message: " << M << std::endl;

                    ofstream result_file("result.txt");
                    result_file << "MitM found message: " << M << std::endl;
                    return steady_clock::now() - timer;
                }
            }

            for (auto it = T_map_copy.begin(); it != T_map_copy.end(); ++it) {
                mpz_class C_S = (C * modinv(it->first, N)) % N;
                mpz_class found_val;
                if (S_block.find(C_S, found_val)) {
                    mpz_class M = it->second * found_val;
                    std::cout << "MitM message: " << M << std::endl;

                    ofstream result_file("Optimized_MitM_result.txt");
                    result_file << "MitM found message: " << M << std::endl;
                    return steady_clock::now() - timer;
                }
            }

            std::cout << "End: bn_t = " << bn_t << ", bn_s = " << bn_s << " : " <<
                duration_cast<microseconds>(steady_clock::now() - timer).count() << std::endl;
        }
    }

    std::cout << "Message for MitM not found! Damn it!" << std::endl;
    return steady_clock::now() - timer;
}

// Result structure for thread-safe result sharing
struct MitMResult {
    std::atomic<bool> found{ false };
    mpz_class M;
    mpz_class T;
    mpz_class S;
    std::mutex result_mutex;
};

// Optimized parallel streaming Meet-in-the-Middle attack. For cases ONLY М^е < N
duration<double, micro> Meet_in_the_Middle_attack_streaming_parallel() {
    mpz_class e = static_cast<unsigned long>(E_CONST);

    std::cout << "> Getting values from " << MitM_test_path << std::endl;
    auto test_values = read_txt(MitM_test_path);
    mpz_class N = test_values["N"];
    mpz_class C = test_values["C"];

    std::cout << "> Parallel Streaming Meet-in-the-Middle attack started for:" << std::endl;
    std::cout << "l = " << L_CONST << std::endl;
    std::cout << "N=" << N << std::endl;
    std::cout << "C=" << C << std::endl;

    auto timer = steady_clock::now();

    const size_t half_space = 1ULL << (L_CONST / 2);
    const size_t num_threads = std::thread::hardware_concurrency();
    const size_t chunk_size = half_space / num_threads;  

    MitMResult result;

    // Progress tracking
    std::atomic<size_t> processed{ 0 };

    std::cout << "Using " << num_threads << " threads, chunk size: " << chunk_size << std::endl;
    std::cout << "Search space: [1, " << half_space << "]" << std::endl;

    // Lambda for worker threads
    auto worker = [&](size_t thread_id) {
        const size_t start = thread_id * chunk_size + 1;
        const size_t end = (thread_id == num_threads - 1) ? half_space : ((thread_id + 1) * chunk_size);

        // Thread-local progress counter
        size_t local_processed = 0;

        for (size_t s = start; s <= end && !result.found.load(std::memory_order_acquire); ++s) {
            mpz_class S(static_cast<unsigned long>(s));  // Convert size_t to mpz_class needed for GMP operations
            mpz_class S_e = modpow(S, e, N);
            mpz_class C_over_S_e = (C * modinv(S_e, N)) % N;

            // Check if C/S^e is a perfect e-th root in range [1, half_space]
            mpz_class T_candidate;

            // Use mpz_root for efficient e-th root extraction
            if (mpz_root(T_candidate.get_mpz_t(), C_over_S_e.get_mpz_t(),
                mpz_get_ui(e.get_mpz_t()))) {

                if (T_candidate >= 1 && mpz_fits_ulong_p(T_candidate.get_mpz_t())) {
                    size_t T_as_size_t = mpz_get_ui(T_candidate.get_mpz_t());
                    if (T_as_size_t <= half_space) {
                        // Double-check: T^e should equal C_over_S_e
                        mpz_class verify = modpow(T_candidate, e, N);
                        if (verify == C_over_S_e) {
                            // Update shared result atomically, if found solution 
                            std::lock_guard<std::mutex> lock(result.result_mutex);
                            if (!result.found.load()) {
                                result.M = T_candidate * S;
                                result.T = T_candidate;
                                result.S = S;
                                result.found.store(true, std::memory_order_release);

                                std::cout << "\nSolution found by thread " << thread_id << "!" << std::endl;
                                std::cout << "M = " << result.M << std::endl;
                                std::cout << "T = " << result.T << ", S = " << result.S << std::endl;
                            }
                            return; // Exit thread early
                        }
                    }
                }
            }

            // Progress reporting (every 100K iterations per thread)
            if (++local_processed % 100000 == 0) {
                size_t total_processed = processed.fetch_add(100000, std::memory_order_relaxed) + 100000;
                if (thread_id == 0) { // Only thread 0 reports about progress
                    double progress = (double)total_processed / (double)half_space * 100.0;
                    auto elapsed = duration_cast<seconds>(steady_clock::now() - timer).count();
                    std::cout << "\rProgress: " << std::fixed << std::setprecision(2)
                        << progress << "% (" << total_processed << "/" << half_space
                        << "), Time: " << elapsed << "s" << std::flush;
                }
            }
        }
    };

    // Launch worker threads
    std::vector<std::thread> threads;
    threads.reserve(num_threads);

    for (size_t i = 0; i < num_threads; ++i)
        threads.emplace_back(worker, i);

    // Wait for all threads to complete
    for (auto& thread : threads)
        thread.join();

    std::cout << std::endl;

    if (result.found.load()) {
        std::cout << "> MitM message: " << result.M << std::endl;

        // Verify the solution
        mpz_class verify_M = modpow(result.M, e, N);
        if (verify_M == C)
            std::cout << "> Solution verified: M^e equiv C (mod N)" << std::endl;
        else
            std::cout << "> Warning: Solution verification failed!" << std::endl;

        // Save result to file
        std::ofstream result_file("mitm_streaming_result.txt");
        result_file << "MitM found message: " << result.M << std::endl;
        result_file << "T = " << result.T << std::endl;
        result_file << "S = " << result.S << std::endl;
        result_file << "Verification: M^e mod N = " << verify_M << std::endl;
        result_file << "Expected C = " << C << std::endl;
    }
    else
        std::cout << "> No message found for MitM in streaming search!" << std::endl;
    return steady_clock::now() - timer;
}

// ---------------------------------------------------------------------------
// ------------------- Bruteforce attack time measure  -----------------------
// ---------------------------------------------------------------------------

duration<double, micro> bruteforce_try_once_time() {
    mpz_class e(static_cast<unsigned long>(E_CONST));

    std::cout << "Getting values from: '" << MitM_test_path << "'" << std::endl;
    auto test_values = read_txt(MitM_test_path);
    mpz_class N = test_values["N"];
    mpz_class C = test_values["C"];

    random_device rd;
    gmp_randstate_t state;
    gmp_randinit_mt(state);
    gmp_randseed_ui(state, rd());

    mpz_class M;
    mpz_urandomb(M.get_mpz_t(), state, L_CONST);

    auto timer = steady_clock::now();

    M = M + 1;
    mpz_class C_maybe = modpow(M, e, N);
    bool cmp = (C_maybe == C);

    gmp_randclear(state);

    return steady_clock::now() - timer;
}

// ---------------------------------------------------------------------------
// ---------------------------- Main function --------------------------------
// ---------------------------------------------------------------------------

int main() {
    try {
        std::cout << "\n--------------------------------------------\n" << std::endl;

        auto SE_time = Small_Exponent_attack();
        std::cout << "'Small exponent' execution time: " << SE_time.count() << " mu s" <<std::endl;

        // auto simple_MitM_time = Pure_Meet_in_the_Middle_attack();
        // std::cout << "Simple 'Meet in the middle' execution time: " << simple_MitM_time.count() << " mu s" << std::endl;

        // auto optimized_MitM_time = Optimized_Space_Meet_in_the_Middle_attack();
        // std::cout << "Optimized 'Meet in the middle' execution time: " << optimized_MitM_time.count() << " mu s" << std::endl;

        // auto brutforce_time = bruteforce_try_once_time();
        // std::cout << "Bruteforce once try (find next M & get C & make a compare) time:" << brutforce_time.count() << " mu s" << std::endl;

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
