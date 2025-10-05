#include "Header.h"

// ---------------------------------------------------------------------------
// ---------------------------- Main function --------------------------------
// ---------------------------------------------------------------------------

int main() {
    try {
        std::cout << "\n--------------------------------------------\n" << std::endl;

        // auto SE_time = Small_Exponent_attack();
        // std::cout << "'Small exponent' execution time: " << SE_time.count() << " mu s" <<std::endl;

        auto simple_MitM_time = Pure_Meet_in_the_Middle_attack();
        std::cout << "Simple 'Meet in the middle' execution time: " << simple_MitM_time.count() << " mu s" << std::endl;

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
