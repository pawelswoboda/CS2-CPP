#include "catch.hpp"
#include "../mcmf.hxx"

using namespace CS2_CPP;

#include <random>

TEST_CASE( "assignment problem", "[assignment]" ) {
   std::random_device rd;     // only used once to initialise (seed) engine
   std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
   std::uniform_int_distribution<long> uni(0,10); // guaranteed unbiased

   for(int run =0; run<100; ++run) {
      const int num_nodes = 6;
      const int num_arcs = 9;
      MCMF_CS2<> mcf( num_nodes, num_arcs);

      for(int i=0; i<3; ++i) {
         for(int j=0; j<3; ++j) {
            mcf.set_arc(i,3+j,0,1,uni(rng));
         }
      }
      for(int i=0; i<3; ++i) {
         mcf.set_supply_demand_of_node(i,1);
         mcf.set_supply_demand_of_node(3+i,-1);
      }
      mcf.run_cs2();

      // read out solutions
      std::vector<int> flow(9);
      for(int i=0; i<9; ++i) {
         flow[i] = mcf.get_flow(i);
      }
      // construct reverse problems. It should have flow values < num_nodes
      MCMF_CS2<> mcf_reverse( num_nodes, num_arcs);
      for(int i=0; i<3; ++i) {
         for(int j=0; j<3; ++j) {
            REQUIRE(flow[i*3 + j] >= 0);
            REQUIRE(flow[i*3 + j] <= 1);
            if(flow[i*3 + j] == 0)
               mcf_reverse.set_arc(i,3+j,1,1000,mcf.get_cost(i*3+j));
            else
               mcf_reverse.set_arc(3+j,i,1,1000,mcf.get_cost(i*3+j));
         }
      }
      for(int i=0; i<3; ++i) {
         mcf_reverse.set_supply_demand_of_node(i,0);
         mcf_reverse.set_supply_demand_of_node(3+i,0);
      }
      mcf_reverse.run_cs2();
      for(int i=0; i<num_arcs; ++i) {
         REQUIRE(mcf_reverse.get_flow(i) < 3);
      }
   }
}

