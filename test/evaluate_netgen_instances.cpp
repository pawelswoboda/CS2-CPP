#include "catch.hpp"
#include "../mcmf.hxx"

using namespace CS2_CPP;

TEST_CASE("netgen", "[netgen test problems evaluation]") {
   std::vector<std::string> test_files = {
      {"instances/netgen/big1.net"},   
      {"instances/netgen/big2.net"},   
      {"instances/netgen/big3.net"},   
      {"instances/netgen/big4.net"},   
      {"instances/netgen/big5.net"},   
      {"instances/netgen/big6.net"},   
      {"instances/netgen/big7.net"},   
      {"instances/netgen/big8.net"},   
      {"instances/netgen/cap1.net"},   
      {"instances/netgen/cap10.net"},   
      {"instances/netgen/cap11.net"},   
      {"instances/netgen/cap12.net"},   
      {"instances/netgen/cap13.net"},   
      {"instances/netgen/cap14.net"},   
      {"instances/netgen/cap15.net"},   
      {"instances/netgen/cap16.net"},   
      {"instances/netgen/cap17.net"},   
      {"instances/netgen/cap18.net"},   
      {"instances/netgen/cap19.net"},   
      {"instances/netgen/cap2.net"},   
      {"instances/netgen/cap20.net"},   
      {"instances/netgen/cap21.net"},   
      {"instances/netgen/cap22.net"},   
      {"instances/netgen/cap23.net"},   
      {"instances/netgen/cap24.net"},   
      {"instances/netgen/cap25.net"},   
      {"instances/netgen/cap26.net"},   
      {"instances/netgen/cap27.net"},   
      {"instances/netgen/cap28.net"},   
      {"instances/netgen/cap29.net"},   
      {"instances/netgen/cap3.net"},   
      {"instances/netgen/cap30.net"},   
      {"instances/netgen/cap31.net"},   
      {"instances/netgen/cap32.net"},   
      {"instances/netgen/cap33.net"},   
      {"instances/netgen/cap34.net"},   
      {"instances/netgen/cap35.net"},   
      {"instances/netgen/cap36.net"},   
      {"instances/netgen/cap37.net"},   
      {"instances/netgen/cap38.net"},   
      {"instances/netgen/cap39.net"},   
      {"instances/netgen/cap4.net"},   
      {"instances/netgen/cap40.net"},   
      {"instances/netgen/cap41.net"},   
      {"instances/netgen/cap5.net"},   
//      {"instances/netgen/cap6.net"},    // file is faulty
      {"instances/netgen/cap7.net"},   
      {"instances/netgen/cap8.net"},   
      {"instances/netgen/cap9.net"},   
      {"instances/netgen/stndrd1.net"},   
      {"instances/netgen/stndrd10.net"},   
      {"instances/netgen/stndrd16.net"},   
      {"instances/netgen/stndrd17.net"},   
      {"instances/netgen/stndrd18.net"},   
      {"instances/netgen/stndrd19.net"},   
      {"instances/netgen/stndrd2.net"},   
      {"instances/netgen/stndrd20.net"},   
      {"instances/netgen/stndrd21.net"},   
      {"instances/netgen/stndrd22.net"},   
      {"instances/netgen/stndrd23.net"},   
      {"instances/netgen/stndrd24.net"},   
      {"instances/netgen/stndrd25.net"},   
      {"instances/netgen/stndrd26.net"},   
      {"instances/netgen/stndrd27.net"},   
      {"instances/netgen/stndrd28.net"},   
      {"instances/netgen/stndrd29.net"},   
      {"instances/netgen/stndrd3.net"},   
      {"instances/netgen/stndrd30.net"},   
      {"instances/netgen/stndrd31.net"},   
      {"instances/netgen/stndrd32.net"},   
      {"instances/netgen/stndrd33.net"},   
      {"instances/netgen/stndrd34.net"},   
      {"instances/netgen/stndrd35.net"},   
      {"instances/netgen/stndrd36.net"},   
      {"instances/netgen/stndrd37.net"},   
      {"instances/netgen/stndrd38.net"},   
      {"instances/netgen/stndrd39.net"},   
      {"instances/netgen/stndrd4.net"},   
      {"instances/netgen/stndrd40.net"},   
      {"instances/netgen/stndrd41.net"},   
      {"instances/netgen/stndrd42.net"},   
      {"instances/netgen/stndrd43.net"},   
      {"instances/netgen/stndrd44.net"},   
      {"instances/netgen/stndrd45.net"},   
      {"instances/netgen/stndrd46.net"},   
      {"instances/netgen/stndrd47.net"},   
      {"instances/netgen/stndrd48.net"},   
      {"instances/netgen/stndrd5.net"},   
      {"instances/netgen/stndrd50.net"},   
      {"instances/netgen/stndrd51.net"},   
      {"instances/netgen/stndrd52.net"},   
      {"instances/netgen/stndrd53.net"},   
      {"instances/netgen/stndrd54.net"},   
      {"instances/netgen/stndrd6.net"},   
      {"instances/netgen/stndrd7.net"},   
      {"instances/netgen/stndrd8.net"},   
      {"instances/netgen/stndrd9.net"},   
      {"instances/netgen/transp1.net"},   
      {"instances/netgen/transp10.net"},   
      {"instances/netgen/transp11.net"},   
      {"instances/netgen/transp12.net"},   
      {"instances/netgen/transp13.net"},   
      {"instances/netgen/transp14.net"},   
      {"instances/netgen/transp2.net"},   
      {"instances/netgen/transp3.net"},   
      {"instances/netgen/transp4.net"},   
      {"instances/netgen/transp5.net"},   
      {"instances/netgen/transp6.net"},   
      {"instances/netgen/transp7.net"},   
      {"instances/netgen/transp8.net"},   
      {"instances/netgen/transp9.net"}
   };
   for(auto& file : test_files) {
      std::cout << " testing file " << file << "\n";
      MCMF_CS2_STAT<> mcf(file);
      mcf.run_cs2();
   }
}
