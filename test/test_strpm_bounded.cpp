#include <cassert>
#include <fstream>
#include <sstream>

#include "oink/oink.hpp"
#include "oink/pgparser.hpp"
#include "solvers/strpm.hpp"
#include "solvers/strpm_simd.hpp"
#include "verifier.hpp"

// Exercise the bounded kernel directly: the adaptive driver may select (2,1)
// and finish before reaching the formerly incomplete (2,2) host.
template <typename Kernel>
class FixedSTRPM : public Kernel
{
public:
    using Kernel::Kernel;

    void run_pair(int t)
    {
        this->Q.resize(this->nodecount());
        this->dirty.resize(this->nodecount());
        Kernel::run(t, 2, 2, 1);
        Kernel::run(t, 2, 2, 0);
    }
};

template <typename Kernel>
static void check(const char* path, int t)
{
    std::ifstream input(path);
    assert(input.good());
    auto game = pg::PGParser::parse_pgsolver(input, false);
    std::ostringstream log;
    pg::Oink oink(game, log);
    FixedSTRPM<Kernel> solver(oink, game);
    solver.run_pair(t);
    assert(game.vertexcount() == 32 && game.count_unsolved() == 0);
    for (int v = 0; v < 32; ++v) assert(game.getWinner(v) == (v == 31));
    pg::Verifier verifier(game, log);
    verifier.verify(true, true, true);
}

int main(int argc, char** argv)
{
    assert(argc == 3);
    const std::string backend = argv[2];
    assert(backend == "scalar" || backend == "simd");
    for (int t : {1, 2}) {
        if (backend == "scalar") check<pg::STRPMSolver>(argv[1], t);
        else check<pg::STRPM_SIMDSolver>(argv[1], t);
    }
}
