#include "catch.hpp"
#include "../mcmf.hxx"

using namespace CS2_CPP;


TEST_CASE( "gte", "[gte test problems evaluation]" ) {
   std::vector<std::string> test_files = {
      {"instances/gte/gte_bad.20"},
      {"instances/gte/gte_bad.40"},
      {"instances/gte/gte_bad.60"},
      {"instances/gte/gte_bad.160"},
      {"instances/gte/gte_bad.200"},
      {"instances/gte/gte_bad.460"},
      {"instances/gte/gte_bad.510"},
      {"instances/gte/gte_bad.1160"},
      {"instances/gte/gte_bad.1700"},
      {"instances/gte/gte_bad.6410"},
      {"instances/gte/gte_bad.6830"},
      {"instances/gte/gte_bad.15100"},
      {"instances/gte/gte_bad.15710"},
      {"instances/gte/gte_bad.35620"},
      {"instances/gte/gte_bad.49320"},
      {"instances/gte/gte_bad.60090"},
      {"instances/gte/gte_bad.65330"},
      {"instances/gte/gte_bad.176280"},
      {"instances/gte/gte_bad.298300"},
      {"instances/gte/gte_bad.451760"},
      {"instances/gte/gte_bad.469010"},
      {"instances/gte/gte_bad.508829"}
   };
   for(auto& file : test_files) {
      std::cout << " testing file " << file << "\n";
      MCMF_CS2_STAT<> mcf(file);
      mcf.run_cs2();
   }
}

