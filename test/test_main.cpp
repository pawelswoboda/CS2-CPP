#define CATCH_CONFIG_MAIN 
#include "catch.hpp"
#include "../mcmf.h"

using namespace CS2_CPP;
using namespace std;

TEST_CASE( "test problem", "[test problem]" ) {
	const int num_nodes = 6;
	const int num_arcs = 8;
	MCMF_CS2 mcf( num_nodes, num_arcs);

   mcf.set_arc( 0, 1, 0, 4, 1);
   mcf.set_arc( 0, 2, 0, 8, 5);
   mcf.set_arc( 1, 2, 0, 5, 0);
   mcf.set_arc( 2, 4, 0, 10, 1);
   mcf.set_arc( 3, 1, 0, 8, 1);
   mcf.set_arc( 3, 5, 0, 8, 1);
   mcf.set_arc( 4, 3, 0, 8, 0);
   mcf.set_arc( 4, 5, 0, 8, 9);
   mcf.set_supply_demand_of_node( 0, 10);
   mcf.set_supply_demand_of_node( 5, -10);

   SECTION("solve mcf") {
      long obj = mcf.run_cs2();
      REQUIRE(obj == 70);
      // after solving, arcs are reordered lexicographically. Note that reverse arc is implicitly added by set_arc as well.
      for(int e=0; e<16; ++e) {
         std::cout << mcf.get_arc_tail(e) << " -> " << mcf.get_arc_head(e) << "; cost = " << mcf.get_cost(e) << "\n";
      }
      //REQUIRE(mcf.compute_objective_cost() == 70);

      // correct arc ordering
      REQUIRE(mcf.get_arc_tail(0) == 0); REQUIRE(mcf.get_arc_head(0) == 1);
      REQUIRE(mcf.get_arc_tail(1) == 0); REQUIRE(mcf.get_arc_head(1) == 2);
      REQUIRE(mcf.get_arc_tail(2) == 1); REQUIRE(mcf.get_arc_head(2) == 0);
      REQUIRE(mcf.get_arc_tail(3) == 1); REQUIRE(mcf.get_arc_head(3) == 2);
      REQUIRE(mcf.get_arc_tail(4) == 1); REQUIRE(mcf.get_arc_head(4) == 3);
      REQUIRE(mcf.get_arc_tail(5) == 2); REQUIRE(mcf.get_arc_head(5) == 0);
      REQUIRE(mcf.get_arc_tail(6) == 2); REQUIRE(mcf.get_arc_head(6) == 1);
      REQUIRE(mcf.get_arc_tail(7) == 2); REQUIRE(mcf.get_arc_head(7) == 4);
      REQUIRE(mcf.get_arc_tail(8) == 3); REQUIRE(mcf.get_arc_head(8) == 1);
      REQUIRE(mcf.get_arc_tail(9) == 3); REQUIRE(mcf.get_arc_head(9) == 4);
      REQUIRE(mcf.get_arc_tail(10) == 3); REQUIRE(mcf.get_arc_head(10) == 5);
      REQUIRE(mcf.get_arc_tail(11) == 4); REQUIRE(mcf.get_arc_head(11) == 2);
      REQUIRE(mcf.get_arc_tail(12) == 4); REQUIRE(mcf.get_arc_head(12) == 3);
      REQUIRE(mcf.get_arc_tail(13) == 4); REQUIRE(mcf.get_arc_head(13) == 5);
      REQUIRE(mcf.get_arc_tail(14) == 5); REQUIRE(mcf.get_arc_head(14) == 3);
      REQUIRE(mcf.get_arc_tail(15) == 5); REQUIRE(mcf.get_arc_head(15) == 4);

      // correct flow values
      // here I must know the original capacities
      //REQUIRE(mcf.get_flow(0) == 4);
      //REQUIRE(mcf.get_flow(1) == 6);
      //REQUIRE(mcf.get_flow(2) == 4);
      //REQUIRE(mcf.get_flow(3) == 10);
      //REQUIRE(mcf.get_flow(4) == 8);
      //REQUIRE(mcf.get_flow(5) == 0);
      //REQUIRE(mcf.get_flow(6) == 8);
      //REQUIRE(mcf.get_flow(7) == 2);

      REQUIRE(4-mcf.get_residual_flow(0) == 4);
      REQUIRE(8-mcf.get_residual_flow(1) == 6);
      REQUIRE(5-mcf.get_residual_flow(3) == 4);
      REQUIRE(10-mcf.get_residual_flow(7) == 10);
      REQUIRE(8-mcf.get_residual_flow(8) == 0);
      REQUIRE(8-mcf.get_residual_flow(10) == 8);
      REQUIRE(8-mcf.get_residual_flow(12) == 8);
      REQUIRE(8-mcf.get_residual_flow(13) == 2);

      REQUIRE(mcf.get_flow(0) == 4);
      REQUIRE(mcf.get_flow(1) == 6);
      REQUIRE(mcf.get_flow(3) == 4);
      REQUIRE(mcf.get_flow(7) == 10);
      REQUIRE(mcf.get_flow(8) == 0);
      REQUIRE(mcf.get_flow(10) == 8);
      REQUIRE(mcf.get_flow(12) == 8);
      REQUIRE(mcf.get_flow(13) == 2);

      // complementary slackness
      REQUIRE(mcf.get_reduced_cost(0) <= 0);
      REQUIRE(mcf.get_reduced_cost(1) == 0);
      REQUIRE(mcf.get_reduced_cost(3) == 0);
      REQUIRE(mcf.get_reduced_cost(7) <= 0);
      REQUIRE(mcf.get_reduced_cost(8) >= 0);
      REQUIRE(mcf.get_reduced_cost(10) <= 0);
      REQUIRE(mcf.get_reduced_cost(12) <= 0);
      REQUIRE(mcf.get_reduced_cost(13) == 0);
   }
}

#include <random>

std::random_device rd;     // only used once to initialise (seed) engine
std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
std::uniform_int_distribution<long> uni(0,10); // guaranteed unbiased

TEST_CASE( "assignment problem", "[assignment]" ) {
   const int num_nodes = 6;
   const int num_arcs = 9;
   MCMF_CS2 mcf( num_nodes, num_arcs);

   for(int i=0; i<3; ++i) {
      for(int j=0; j<3; ++j) {
         mcf.set_arc(i,3+j,0,1,uni(rng));
      }
   }
   for(int i=0; i<3; ++i) {
      mcf.set_supply_demand_of_node(i,1);
      mcf.set_supply_demand_of_node(3+i,-1);
   }

   // read out solutions
   std::vector<int> flow(9);
   for(int i=0; i<9; ++i) {
      flow[i] = mcf.get_flow(i);
   }
   // construct reverse problems. It should have flow values < num_nodes
   MCMF_CS2 mcf_reverse( num_nodes, num_arcs);
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

