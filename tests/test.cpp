#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "sobolseqgenerator.hpp"

TEST_CASE( "Sobol sequence generator basic test (pass)", "[single-file]")
{
    SobolSeqGenerator net_gen;
    net_gen.Init(10, 3, "new-joe-kuo-6.21201.txt");
    REQUIRE(net_gen.GeneratePoint().coordinate == std::vector<Real>({0, 0, 0}));
    REQUIRE(net_gen.GeneratePoint().coordinate == std::vector<Real>({0.5, 0.5, 0.5}));
    REQUIRE(net_gen.GeneratePoint().coordinate == std::vector<Real>({0.75, 0.25, 0.25}));
    REQUIRE(net_gen.GeneratePoint().coordinate == std::vector<Real>({0.25, 0.75, 0.75}));
    REQUIRE(net_gen.GeneratePoint().coordinate == std::vector<Real>({0.375, 0.375, 0.625}));
    REQUIRE(net_gen.GeneratePoint().coordinate == std::vector<Real>({0.875, 0.875, 0.125}));
    REQUIRE(net_gen.GeneratePoint().coordinate == std::vector<Real>({0.625, 0.125, 0.875}));
    REQUIRE(net_gen.GeneratePoint().coordinate == std::vector<Real>({0.125, 0.625, 0.375}));
    REQUIRE(net_gen.GeneratePoint().coordinate == std::vector<Real>({0.1875, 0.3125, 0.9375}));
    REQUIRE(net_gen.GeneratePoint().coordinate == std::vector<Real>({0.6875, 0.8125, 0.4375}));
}
