#ifndef STOCHASTIC_HPP
#define STOCHASTIC_HPP

#include <cstdlib>
#include <algorithm>
#include <functional>
#include <vector>

#include <limits>
#include <type_traits>

namespace sto {

    template <typename T>
    using Random_Generator = std::function<auto () -> T>;

    template <typename T = double>
    auto random() -> T {
        return static_cast<T>(std::rand());
    };

    template <typename T>
    auto random_between(T init, T end) -> T {
        return random<T>() % (end - init) + init;
    }

    template<typename T = double>
    auto rng() -> T {
        return static_cast<T>(std::rand() / static_cast<double>(std::numeric_limits<std::invoke_result_t<decltype(std::rand)>>::max()));
    }

    template <typename T>
    auto n_unique_random(std::size_t n, Random_Generator<T> gen = random<T>) -> std::vector<T> {
        std::vector<T> randoms;
        randoms.reserve(n);

        std::generate_n(std::back_inserter(randoms), n, [gen, begin = std::begin(randoms), last = std::begin(randoms)]() mutable {
            auto guess = gen();

            while (std::find(begin, last, guess) != last) {
                guess = gen();
            }

            ++last;

            return guess;
        });

        return randoms;
    }

}


#endif // STOCHASTIC_HPP
