#ifndef _MCMF_H_
#define _MCMF_H_

#include <cmath>
#include <assert.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <memory>
#include <chrono>

namespace CS2_CPP {

////////////////////////////////////////////////////////////////////////////////
//
// MCMF_CS2
//
////////////////////////////////////////////////////////////////////////////////

template<typename COST_TYPE = long long int, typename CAPACITY_TYPE = long int>
class MCMF_CS2
{
protected:
   constexpr static long MAX_64  = 0x7fffffffffffffffLL;
   constexpr static long MAX_32 = 0x7fffffff;
   constexpr static long PRICE_MAX = MAX_64;

   // parameters
   constexpr static double UPDT_FREQ = 0.4;
   constexpr static long UPDT_FREQ_S = 30;
   constexpr static double SCALE_DEFAULT = 12.0;
   // PRICE_OUT_START may not be less than 1
   constexpr static long PRICE_OUT_START = 1;
   constexpr static double CUT_OFF_POWER = 0.44;
   constexpr static double CUT_OFF_COEF = 1.5;
   constexpr static double CUT_OFF_POWER2 = 0.75;
   constexpr static double CUT_OFF_COEF2 = 1;
   constexpr static double CUT_OFF_GAP = 0.8;
   constexpr static double CUT_OFF_MIN = 12;
   constexpr static double CUT_OFF_INCREASE = 4;

   constexpr static long TIME_FOR_PRICE_IN1 = 2;
   constexpr static long TIME_FOR_PRICE_IN2 = 4;
   constexpr static long TIME_FOR_PRICE_IN3 = 6;

   constexpr static long MAX_CYCLES_CANCELLED = 0;
   constexpr static long START_CYCLE_CANCEL = 100;

 public:
	typedef long long int excess_t;
	typedef long long int price_t;

	class NODE;
	
	class ARC {
	public:
		long _rez_capacity; // residual capacity;
		price_t _cost; // cost of arc;
		NODE *_head; // head node;
		ARC *_sister; // opposite arc;
	public:
		ARC() {}
		~ARC() {}

		void set_rez_capacity( long rez_capacity) { _rez_capacity = rez_capacity; }
		void dec_rez_capacity( long delta) { _rez_capacity -= delta; }
		void inc_rez_capacity( long delta) { _rez_capacity += delta; }
		void set_cost( price_t cost) { _cost = cost; }
		void multiply_cost( price_t mult) { _cost *= mult; }
		void set_head( NODE *head) { _head = head; }
		void set_sister( ARC *sister) { _sister = sister; }
		long rez_capacity() { return _rez_capacity; }
		price_t cost() { return _cost; }
		NODE *head() { return _head; }
		ARC *sister() { return _sister; }
	};

	class NODE {
	public:
		excess_t _excess; // excess of the node;
		price_t _price; // distance from a sink;
		ARC *_first; // first outgoing arc;
		ARC *_current; // current outgoing arc;
		ARC *_suspended;
		NODE *_q_next; // next node in push queue
		NODE *_b_next; // next node in bucket-list;
		NODE *_b_prev; // previous node in bucket-list;
		long _rank; // bucket number;
		long _inp; // auxilary field;
	public:
		NODE() {}
      ~NODE() {}

		void set_excess( excess_t excess) { _excess = excess; }
		void dec_excess( long delta) { _excess -= delta; }
		void inc_excess( long delta) { _excess += delta; }
		void set_price( price_t price) { _price = price; }
		void dec_price( long delta) { _price -= delta; }
		void set_first( ARC *first) { _first = first; }
		void set_current( ARC *current) { _current = current; }
		void inc_current() { _current ++; }
		void set_suspended( ARC *suspended) { _suspended = suspended; }
		void set_q_next( NODE *q_next) { _q_next = q_next; }
		void set_b_next( NODE *b_next) { _b_next = b_next; }
		void set_b_prev( NODE *b_prev) { _b_prev = b_prev; }
		void set_rank( long rank) { _rank = rank; }
		void set_inp( long inp) { _inp = inp; }
		excess_t excess() { return _excess; }
		price_t price() { return _price; }
		ARC *first() { return _first; }
		void dec_first() { _first --; }
		void inc_first() { _first ++; }
		ARC *current() { return _current; }
		ARC *suspended() { return _suspended; }
		NODE *q_next() { return _q_next; }
		NODE *b_next() { return _b_next; }
		NODE *b_prev() { return _b_prev; }
		long rank() { return _rank; }
		long inp() { return _inp; }
	};
 
	class BUCKET {
	private:
		// 1st node with positive excess or simply 1st node in the buket;
		NODE *_p_first;
	public:
	BUCKET( NODE *p_first) : _p_first(p_first) {}
		BUCKET() {}
		~BUCKET() {}

		void set_p_first( NODE *p_first) { _p_first = p_first; }
		NODE *p_first() { return _p_first; }
	};

 public:
	long _n; // number of nodes
	long _m; // number of arcs

   // counters for operations.
   long _n_rel; // number of relabels from last price update
   long _n_ref; // current number of refines
   long _n_src; // current number of nodes with excess
   long _n_bad_pricein;
   long _n_bad_relabel;

   std::unique_ptr<long[]> _cap;
	//long *_cap; // array containig capacities
	//NODE *_nodes; // array of nodes
   std::unique_ptr<NODE[]> _nodes; // array of nodes
	NODE *_sentinel_node; // next after last
	NODE *_excq_first; // first node in push-queue
	NODE *_excq_last; // last node in push-queue
	//ARC *_arcs; // array of arcs
   std::unique_ptr<ARC[]> _arcs;
	ARC *_sentinel_arc; // next after last

	//BUCKET *_buckets; // array of buckets
   std::unique_ptr<BUCKET[]> _buckets; // array of buckets
	BUCKET *_l_bucket; // last bucket
	long _linf; // number of l_bucket + 1
	int _time_for_price_in;

	price_t _epsilon; // quality measure
	price_t _dn; // cost multiplier = number of nodes + 1
	price_t _price_min; // lowest bound for prices
	price_t _mmc; // multiplied maximal cost
	double _f_scale; // scale factor
	double _cut_off_factor; // multiplier to produce cut_on and cut_off from n and epsilon
	double _cut_on; // the bound for returning suspended arcs
	double _cut_off; // the bound for suspending arcs
	excess_t _total_excess; // total excess

	// if = 1 - signal to start price-in ASAP - 
	// maybe there is infeasibility because of susoended arcs 
	int _flag_price;
	// if = 1 - update failed some sources are unreachable: either the 
	// problem is unfeasible or you have to return suspended arcs
	int _flag_updt;
	// maximal number of cycles cancelled during price refine 
	int _snc_max;

	// dummy variables;
	ARC _d_arc; // dummy arc - for technical reasons
	NODE _d_node; // dummy node - for technical reasons
	NODE *_dummy_node; // the address of d_node
   std::unique_ptr<NODE> _dnode;

	// sketch variables used during reading in arcs;
	long _node_min; // minimal no of nodes
	long _node_max; // maximal no of nodes
   std::unique_ptr<long[]> _arc_first; // internal array for holding
                     // - node degree
                     // - position of the first outgoing arc
	std::unique_ptr<long[]> _arc_tail; // internal array: tails of the arcs
	long _pos_current;
	ARC *_arc_current;
	ARC *_arc_new;
	ARC *_arc_tmp;
	price_t _max_cost; // maximum cost
	excess_t _total_p; // total supply
	excess_t _total_n; // total demand
	// pointers to the node structure
	NODE *_i_node;
	NODE *_j_node;


 public:
	MCMF_CS2( long num_nodes, long num_arcs) {
		_n = num_nodes;
		_m = num_arcs;

      _n_bad_pricein = 0;
      _n_bad_relabel = 0;

		_flag_price = 0;
		_flag_updt = 0;
      // allocate arrays and prepare for "receiving" arcs;
		// will also reset _pos_current, etc.;
		allocate_arrays();
	}

   MCMF_CS2(const std::string& file);

	~MCMF_CS2() {
      deallocate_arrays();
   }

	void err_end( int cc);
	void allocate_arrays();
	void deallocate_arrays();
	void set_arc( long tail_node_id, long head_node_id,
				  long low_bound, long up_bound, price_t cost);
	void set_supply_demand_of_node( long id, long excess);
	void pre_processing();
	void cs2_initialize();
	void up_node_scan( NODE *i);
	void price_update();
	int relabel( NODE *i);
	void discharge( NODE *i);
	int price_in();
	void refine();
	int price_refine();
	void compute_prices();
	void price_out();
	int update_epsilon();
	void init_solution();
	void cs_cost_reinit();
	long long int cs2_cost_restart();
	void finishup( double *objective_cost);
	void cs2( double *objective_cost);
	price_t run_cs2();

   // information functions
   long no_nodes() const { return _n; }
   long no_arcs() const { return 2*_m; }
   long no_arcs(const long i) const { return _nodes[i].first() - _nodes[i+1].first(); }
   long starting_arc(const long i) const { return N_ARC(_nodes[i].first()); } // index of the first arc emanating out from node i
   price_t compute_objective_cost() const;
   long get_arc_tail(const long arc_id) const { return N_NODE(_arcs[arc_id].sister()->head()); }
   long get_arc_head(const long arc_id) const { return N_NODE(_arcs[arc_id].head()); }
   long get_flow(const long arc_id) const { return _cap[arc_id] - _arcs[arc_id].rez_capacity(); }
   long get_residual_flow(const long arc_id) const { return _arcs[arc_id].rez_capacity(); } // original flow i capacity - residual flow
   long get_cost(const long arc_id) const { return _arcs[arc_id].cost(); }
   long get_reduced_cost(const long arc_id) const { return get_cost(arc_id) + get_price(get_arc_tail(arc_id)) - get_price(get_arc_head(arc_id)); }
   long get_price(const long node_id) const { return _nodes[node_id].price(); }
   void update_cost(const long arc_id, const price_t cost) { 
      _arcs[arc_id].set_cost(cost);
      _arcs[arc_id].sister()->set_cost(-cost);
   }
   void set_cap(const long arc_id, const long cap) { assert(false); assert(cap >= 0); _cap[arc_id] = cap; } // me must adjust residual capacity as well

   //#define N_NODE( i ) ( ( (i) == NULL ) ? -1 : ( (i) - _nodes + _node_min ) )
   long N_NODE(NODE* n) const { return n - _nodes.get() + _node_min; }
   //#define N_ARC( a ) ( ( (a) == NULL )? -1 : (a) - _arcs )
   long N_ARC(ARC* a) const { return a - _arcs.get(); }

 
   // loop utilities
   template<typename LAMBDA>
   void for_each_node(LAMBDA f) {
      for(NODE* i= _nodes.get(); i<_nodes.get() + _n; ++i) {
         f(i);
      }
   }
   template<typename LAMBDA>
   void for_each_arc(NODE* i, LAMBDA f) {
      for(ARC* a = i->first(); a<i->suspended(); ++a) {
         f(a);
      }
   }
   template<typename LAMBDA>
   void for_each_arc(long i, LAMBDA f) {
      for(ARC* a = _nodes[i].first(); a<_nodes[i+1].first() ; ++a) {
         f(a);
      }
   }

	// shared utils;
	void increase_flow( NODE *i, NODE *j, ARC *a, long df) {
		i->dec_excess( df);
		j->inc_excess( df);
		a->dec_rez_capacity( df);
		a->sister()->inc_rez_capacity( df);
	}
	bool time_for_update() { 
		return ( _n_rel > _n * UPDT_FREQ + _n_src * UPDT_FREQ_S); 
	}
	// utils for excess queue;
	void reset_excess_q() {
		for ( ; _excq_first != NULL; _excq_first = _excq_last ) {
			_excq_last = _excq_first->q_next();
			_excq_first->set_q_next( _sentinel_node);
		}
	}
	bool out_of_excess_q( NODE *i) { return ( i->q_next() == _sentinel_node); }
	bool empty_excess_q() { return ( _excq_first == NULL); }
	bool nonempty_excess_q() { return ( _excq_first != NULL); }
	void insert_to_excess_q( NODE *i) {
		if ( nonempty_excess_q() ) {
			_excq_last->set_q_next( i);
		} else {
			_excq_first = i;
		}
		i->set_q_next( NULL);
		_excq_last = i;
	}
	void insert_to_front_excess_q( NODE *i) {
		if ( empty_excess_q() ) {
			_excq_last = i;
		}
		i->set_q_next( _excq_first);
		_excq_first = i;
	}
	void remove_from_excess_q( NODE *i) {
		i = _excq_first;
		_excq_first = i->q_next();
		i->set_q_next( _sentinel_node);
	}
	// utils for excess queue as a stack;
	bool empty_stackq() { return empty_excess_q(); }
	bool nonempty_stackq() { return nonempty_excess_q(); }
	void reset_stackq() { reset_excess_q(); }
	void stackq_push( NODE *i) {
		i->set_q_next( _excq_first);
		_excq_first = i;
	}
	void stackq_pop( NODE *i) {
		remove_from_excess_q( i);
	}
	// utils for buckets;
	void reset_bucket( BUCKET *b) { b->set_p_first( _dnode.get()); }
	bool nonempty_bucket( BUCKET *b) { return ( (b->p_first()) != _dnode.get()); }
	void insert_to_bucket( NODE *i, BUCKET *b) {
		i->set_b_next( b->p_first() );
		b->p_first()->set_b_prev( i);
		b->set_p_first( i);
	}
	void get_from_bucket( NODE *i, BUCKET *b) {
		i = b->p_first();
		b->set_p_first( i->b_next());
	}
	void remove_from_bucket( NODE *i, BUCKET *b) {
		if ( i == b->p_first() ) {
			b->set_p_first( i->b_next());
		} else {
			i->b_prev()->set_b_next( i->b_next());
			i->b_next()->set_b_prev( i->b_prev());
		}
	}
	// misc utils;
	void update_cut_off() {
		if ( _n_bad_pricein + _n_bad_relabel == 0) {
			_cut_off_factor = CUT_OFF_COEF2 * std::pow( _n, CUT_OFF_POWER2 );
			_cut_off_factor = std::max ( _cut_off_factor, double(CUT_OFF_MIN) );
			_cut_off = _cut_off_factor * _epsilon;
			_cut_on = _cut_off * CUT_OFF_GAP;
		} else {
			_cut_off_factor *= CUT_OFF_INCREASE;
			_cut_off = _cut_off_factor * _epsilon;
			_cut_on = _cut_off * CUT_OFF_GAP;
		}
	}
	void exchange( ARC *a, ARC *b) {
		if ( a != b) {
			ARC *sa = a->sister();
			ARC *sb = b->sister();
			long d_cap;						

			_d_arc.set_rez_capacity( a->rez_capacity());
			_d_arc.set_cost( a->cost());
			_d_arc.set_head( a->head());

			a->set_rez_capacity( b->rez_capacity());
			a->set_cost( b->cost());
			a->set_head( b->head());

			b->set_rez_capacity( _d_arc.rez_capacity());
			b->set_cost( _d_arc.cost());
			b->set_head( _d_arc.head());

			if ( a != sb) {			
				b->set_sister( sa);
				a->set_sister( sb);
				sa->set_sister( b);
				sb->set_sister( a);
			}
						
			d_cap = _cap[N_ARC(a)];
			_cap[N_ARC(a)] = _cap[N_ARC(b)];
			_cap[N_ARC(b)] = d_cap;	
		}
	}
};

// class 
// - counting how many operations were performed 
// - checking feasibility of solutions after optimization 
// - printing diagnostics and statistics after optimization
// - printing graph and obtained flow
template<typename COST_TYPE = long long int, typename CAPACITY_TYPE = long int>
class MCMF_CS2_STAT : public MCMF_CS2<COST_TYPE, CAPACITY_TYPE> {
public:
   MCMF_CS2_STAT( long num_nodes, long num_arcs) 
      : MCMF_CS2<COST_TYPE,CAPACITY_TYPE>(num_nodes, num_arcs)
        /*
        ,
      _n_push(0),
      _n_relabel(0),
      _n_discharge(0),
      _n_refine(0),
      _n_update(0),
      _n_scan(0),
      _n_prscan(0),
      _n_prscan1(0),
      _n_prscan2(0),
      _n_prefine(0),
      _no_zero_cycles(false),
      _check_solution(false),
      _comp_duals(false),
      _cost_restart(false),
      _print_ans(true)
      */
   {}
   MCMF_CS2_STAT(const std::string& file) : MCMF_CS2<COST_TYPE,CAPACITY_TYPE>(file) 
   {}

   long long int run_cs2();

	void print_solution();
	void print_graph();

	int check_feas();
	int check_cs();
	int check_eps_opt();

private:
   /*
   long _n_push;
   long _n_relabel;
   long _n_discharge;
   long _n_refine;
   long _n_update;
   long _n_scan;
   long _n_prscan;
   long _n_prscan1;
   long _n_prscan2;
   long _n_prefine;
   */

   /*
   bool _no_zero_cycles; // finds an optimal flow with no zero-cost cycles
   bool _check_solution; // check feasibility/optimality. HIGH OVERHEAD!
   bool _comp_duals; // compute prices?
   bool _cost_restart; // to be able to restart after a cost function change -> should go to main class or be deleted
	bool _print_ans;
	long long int *_node_balance;
   */
};


////////////////////////////////////////////////////////////////////////////////
//
// Implementation
//
////////////////////////////////////////////////////////////////////////////////


#define WHITE 0
#define GREY  1
#define BLACK 2
#define OPEN( a ) ( a->rez_capacity() > 0 )
#define CLOSED( a ) ( a->rez_capacity() <= 0 )
#define REDUCED_COST( i, j, a ) ( i->price() + a->cost() - j->price() )
#define FEASIBLE( i, j, a ) ( i->price() + a->cost() < j->price() )
#define ADMISSIBLE( i, j, a ) ( OPEN( a ) && FEASIBLE( i, j, a ) )
#define SUSPENDED( i, a ) ( a < i->first() )


#define REMOVE_FROM_EXCESS_Q( i )				\
	{											\
		i           = _excq_first;				\
		_excq_first = i -> q_next();			\
		i ->set_q_next( _sentinel_node );		\
	}

#define STACKQ_POP( i )							\
	{											\
		i           = _excq_first;				\
		_excq_first = i -> q_next();			\
		i ->set_q_next( _sentinel_node );		\
	}

#define GET_FROM_BUCKET( i, b )					\
	{											\
		i    = ( b -> p_first() );				\
		b ->set_p_first( i -> b_next() );		\
	}
#define REMOVE_FROM_BUCKET( i, b )								\
	{															\
		if ( i == ( b -> p_first() ) )							\
			b ->set_p_first( i -> b_next() );					\
		else													\
			{													\
				( i -> b_prev() )->set_b_next( i -> b_next() );	\
				( i -> b_next() )->set_b_prev( i -> b_prev() );	\
			}													\
	}

template<typename COST_TYPE, typename CAPACITY_TYPE>
MCMF_CS2<COST_TYPE,CAPACITY_TYPE>::MCMF_CS2(const std::string& file)
{
   _n_bad_pricein = 0;
   _n_bad_relabel = 0;

   _flag_price = 0;
   _flag_updt = 0;

   std::ifstream instance;
   instance.open(file);
   if(!instance.is_open()) {
      throw std::runtime_error("could not open file " + file);
   }
   std::string line;
   bool read_problem_def = false;
   while(std::getline(instance, line))
   {
      if(line.empty()) continue;
      std::istringstream iss(line);
      char id;
      if (!(iss >> id )) { throw std::runtime_error("in file " + file + ": cannot read line " + line); } 
      switch(id) {
         case 'c':
            break;
         case 'p': 
            {
               if(read_problem_def == true) { throw std::runtime_error("in file " + file + ": not more than one line beginning with 'p' allowed"); }
               read_problem_def = true;
               std::string min;
               if( !(iss >> min >> _n >> _m)) { throw std::runtime_error("in file " + file + ": cannot read number of nodes and arcs from line:\n " + line); } 
               if("min" != min) { throw std::runtime_error("in file " + file + ": min must come after 'p' in line:\n " + line); } 
               allocate_arrays();
               break;
            }
         case 'n': 
            {
               long id, flow;
               if( !(iss >> id >> flow)) { throw std::runtime_error("in file " + file + ": cannot read number node id and external flow from line:\n " + line); } 
               assert(id >= 1);
               set_supply_demand_of_node(id-1,flow);
               break;
            }
         case 'a': 
            {
               long i,j,lower,upper,cost;
               if( !(iss >> i >> j >> lower >> upper >> cost)) { throw std::runtime_error("in file " + file + ": cannot read arc information from line:\n " + line); } 
               assert(i >= 1);
               assert(j >= 1);
               set_arc(i-1,j-1,lower,upper,cost);
               break;
            }
         default:
            throw std::runtime_error("in file " + file + ": unknown line identifier " + std::to_string(id));
            break;
      }
   }
   if(read_problem_def == false) { throw std::runtime_error("Exactly one line beginning with 'p' allowed in file " + file); }
}

template<typename COST_TYPE, typename CAPACITY_TYPE>
void MCMF_CS2<COST_TYPE,CAPACITY_TYPE>::allocate_arrays()
{
	// (1) allocate memory for 'nodes', 'arcs' and internal arrays;

   // TODO: make unique pointers out of that.
	//_nodes = (NODE*) calloc ( _n+1,   sizeof(NODE) );
   _nodes = std::unique_ptr<NODE[]>{ new NODE[_n+1] };
   if(_nodes.get() == nullptr) { throw std::runtime_error("could not allocate _nodes"); }
   std::fill(&_nodes[0], &_nodes[_n+1], NODE());
	//_arcs = (ARC*)  calloc ( 2*_m+1, sizeof(ARC) );
   _arcs = std::unique_ptr<ARC[]>{ new ARC[2*_m+1] };
   if(_arcs.get() == nullptr) { throw std::runtime_error("could not allocate _arcs"); }
   std::fill(_arcs.get(), _arcs.get()+2*_m+1, ARC());
	//_cap = (long*) calloc ( 2*_m,   sizeof(long) );
   _cap = std::unique_ptr<long[]>{ new long[2*_m] };
   if(_cap.get() == nullptr) { throw std::runtime_error("could not allocate _cap"); }
   std::fill(&_cap[0], &_cap[2*_m], 0);
   _dnode = std::unique_ptr<NODE>{new NODE};
   if(_dnode.get() == nullptr) { throw std::runtime_error("could not allocate _dnode"); }

	_arc_tail = std::unique_ptr<long[]>{ new long[2*_m] }; // _arc_tail = (long*) calloc ( 2*_m,   sizeof(long) );
   if(_arc_tail.get() == nullptr) { throw std::runtime_error("could not allocate _arc_tail"); }
   std::fill(&_arc_tail[0], &_arc_tail[2*_m], 0);
	_arc_first = std::unique_ptr<long[]>{ new long[_n+1] }; //(long*) calloc ( _n+1,   sizeof(long) );
   if(_arc_first.get() == nullptr) { throw std::runtime_error("could not allocate _arc_first"); }
   std::fill(&_arc_first[0], &_arc_first[_n+1], 0);
	// arc_first [ 0 .. n ] = 0 - initialized by calloc;

	for_each_node([](NODE* i) { i->set_excess(0); });

	if ( _nodes == NULL || _arcs == NULL || _arc_first == NULL || _arc_tail == NULL) {
      throw std::runtime_error("Error:  Memory allocation problem inside CS2\n");
	}

	// (2) resets;
	_pos_current = 0;
	_arc_current = _arcs.get(); // set "current" pointer to the first arc
	_node_max = 0;
	_node_min = _n;
	_max_cost = 0;
	_total_p = _total_n = 0;
	// at this moment we are ready to add arcs and build the network,
	// by using set_arc()...
}

template<typename COST_TYPE, typename CAPACITY_TYPE>
void MCMF_CS2<COST_TYPE,CAPACITY_TYPE>::deallocate_arrays()
{ 
	//if ( _arcs) free ( _arcs );
	//if ( _dnode) delete _dnode;
	//if ( _cap) free ( _cap );
	//if ( _buckets) free ( _buckets );
	//if ( _check_solution == true) free ( _node_balance );
   //if ( _nodes) {
   //   _nodes = _nodes - _node_min;
   //   free ( _nodes );
   //}
}

template<typename COST_TYPE, typename CAPACITY_TYPE>
void MCMF_CS2<COST_TYPE,CAPACITY_TYPE>::set_arc( long tail_node_id, long head_node_id,
						long low_bound, long up_bound, // up_bound is basically capacity;
						price_t cost)
{
	// DIMACS format:
	// c arc has <tail> <head> <capacity l.b.> <capacity u.b> <cost>
   assert(_pos_current < 2*_m);

	if ( tail_node_id < 0 || tail_node_id >= _n || 
		 head_node_id < 0 || head_node_id >= _n ) {
      throw std::runtime_error("Error:  Arc with head or tail out of bounds inside CS2");
   }
   if ( up_bound < 0 ) {
      up_bound = MAX_32;
      std::cout << "Warning:  Infinite capacity replaced by BIGGEST_FLOW\n";
   }
   if ( low_bound < 0 || low_bound > up_bound ) {
      throw std::runtime_error("Error:  Wrong capacity bounds inside CS2\n");
	}

	// no of arcs incident to node i is placed in _arc_first[i+1]
	_arc_first[tail_node_id + 1] ++; 
	_arc_first[head_node_id + 1] ++;
	_i_node = _nodes.get() + tail_node_id;
	_j_node = _nodes.get() + head_node_id;

	// store information about the arc
	_arc_tail[_pos_current]   = tail_node_id;
	_arc_tail[_pos_current+1] = head_node_id;
	_arc_current->set_head( _j_node );
	_arc_current->set_rez_capacity( up_bound - low_bound );
	_cap[_pos_current] = up_bound;
	_arc_current->set_cost( cost );
	_arc_current->set_sister( _arc_current + 1 );
	( _arc_current + 1 )->set_head( _nodes.get() + tail_node_id );
	( _arc_current + 1 )->set_rez_capacity( 0 );
	_cap[_pos_current+1] = 0;
	( _arc_current + 1 )->set_cost( -cost );
	( _arc_current + 1 )->set_sister( _arc_current );

	_i_node->dec_excess( low_bound );
	_j_node->inc_excess( low_bound );

	// searching for minimum and maximum node
	if ( head_node_id < _node_min ) _node_min = head_node_id;
	if ( tail_node_id < _node_min ) _node_min = tail_node_id;
	if ( head_node_id > _node_max ) _node_max = head_node_id;
	if ( tail_node_id > _node_max ) _node_max = tail_node_id;

	if ( cost < 0 ) cost = -cost;
	if ( cost > _max_cost && up_bound > 0 ) _max_cost = cost;

	// prepare for next arc to be added;
	_arc_current += 2;
	_pos_current += 2;
}

template<typename COST_TYPE, typename CAPACITY_TYPE>
void MCMF_CS2<COST_TYPE,CAPACITY_TYPE>::set_supply_demand_of_node( long id, long excess)
{
	// set supply and demand of nodes; not used for transhipment nodes;
	if ( id < 0 || id > _n ) {
      throw std::runtime_error("Error:  Unbalanced problem inside CS2\n");
	}
	(_nodes.get() + id)->set_excess( excess);
	if ( excess > 0) _total_p += excess;
	if ( excess < 0) _total_n -= excess;
}

template<typename COST_TYPE, typename CAPACITY_TYPE>
void MCMF_CS2<COST_TYPE,CAPACITY_TYPE>::pre_processing()
{
	// called after the arcs were just added and before run_cs2();
	// ordering arcs - linear time algorithm;
   // arcs are ordered so that they are ordered first by tail node and second by order of insertion. Here original implementation differs.
	long i;
	long last, arc_num, arc_new_num;;
	long tail_node_id;
	NODE *head_p;
	ARC *arc_new, *arc_tmp;	
	long up_bound;
	price_t cost; // arc cost;
	excess_t cap_out; // sum of outgoing capacities
	excess_t cap_in; // sum of incoming capacities

	if ( std::abs( _total_p - _total_n ) > 0.5 ) {
      throw std::runtime_error("Error:  Unbalanced problem inside CS2\n");
	}
	
	// first arc from the first node
	( _nodes.get() + _node_min )->set_first( _arcs.get() );

	// before below loop arc_first[i+1] is the number of arcs outgoing from i;
	// after this loop arc_first[i] is the position of the first 
	// outgoing from node i arcs after they would be ordered;
	// this value is transformed to pointer and written to node.first[i]
	for ( i = _node_min + 1; i <= _node_max + 1; i ++ ) {
		_arc_first[i] += _arc_first[i-1];
		( _nodes.get() + i )->set_first( _arcs.get() + _arc_first[i] );
	}

   std::unique_ptr<long[]> perm = std::unique_ptr<long[]>{new long[_n]}; // for holding permutations for sorting
   if(perm.get() == nullptr) { throw std::runtime_error("could not allocate _perm"); }
	// scanning all the nodes except the last
	for ( i = _node_min; i <= _node_max; i ++ ) {

		last = N_ARC( ( _nodes.get() + i + 1 )->first() );// - _arcs.get();
		// arcs outgoing from i must be cited    
		// from position arc_first[i] to the position
		// equal to initial value of arc_first[i+1]-1

		for ( arc_num = _arc_first[i]; arc_num < last; arc_num ++ ) {
			tail_node_id = _arc_tail[arc_num];

			while ( tail_node_id != i ) {
				// the arc no  arc_num  is not in place because arc cited here
				// must go out from i;
				// we'll put it to its place and continue this process
				// until an arc in this position would go out from i

				arc_new_num = _arc_first[tail_node_id];
				_arc_current = _arcs.get() + arc_num;
				arc_new = _arcs.get() + arc_new_num;
	    
				// arc_current must be cited in the position arc_new    
				// swapping these arcs:

				head_p = arc_new->head();
				arc_new->set_head( _arc_current->head() );
				_arc_current->set_head( head_p );

				up_bound          = _cap[arc_new_num];
				_cap[arc_new_num] = _cap[arc_num];
				_cap[arc_num]     = up_bound;

				up_bound = arc_new->rez_capacity();
				arc_new->set_rez_capacity( _arc_current->rez_capacity() );
				_arc_current->set_rez_capacity( up_bound) ;

				cost = arc_new->cost();
				arc_new->set_cost( _arc_current->cost() );
				_arc_current->set_cost( cost );

				if ( arc_new != _arc_current->sister() ) {
					arc_tmp = arc_new->sister();
					arc_new->set_sister( _arc_current->sister() );
					_arc_current->set_sister( arc_tmp );

					_arc_current->sister()->set_sister( _arc_current );
					arc_new->sister()->set_sister( arc_new );
				}

				_arc_tail[arc_num] = _arc_tail[arc_new_num];
				_arc_tail[arc_new_num] = tail_node_id;

				// we increase arc_first[tail_node_id]
				_arc_first[tail_node_id] ++ ;

            tail_node_id = _arc_tail[arc_num];
         }
      }
      // all arcs outgoing from  i  are in place

      // sort outgoing arcs of every node by head node id.
      // this is not done in the original implementation.
      long s = _nodes[i+1].first() - _nodes[i].first();
      assert(s <= _n);
      for(int c=0; c<s; ++c) {
         perm[c] = c;
      }
      ARC* arc = _nodes[i].first();
      long* cap = &_cap[N_ARC(arc)];
      std::sort(perm.get(), perm.get()+s, 
            [arc](int i, int j) { return arc[i].head() < arc[j].head(); });
      // now follow cycles in permutation. negate permutation entries to avoid permuting back
      for(int c=0; c<s; ++c) {
         long next_idx = perm[c];
         if(next_idx == c || next_idx < 0) {
            continue;
         }
         long cur_idx = c;
         ARC tmp_arc = arc[cur_idx];
         long tmp_cap = cap[cur_idx];
         while(next_idx >= 0) {
            std::swap(tmp_arc, arc[next_idx]); // now tmp holds next element
            std::swap(tmp_cap, cap[next_idx]); // now tmp holds next element
            perm[cur_idx] -= s; // mark as visited
            cur_idx = next_idx;
            next_idx = perm[next_idx];
         }
         std::swap(tmp_arc, arc[s+next_idx]); // now tmp holds next element
         std::swap(tmp_cap, cap[s+next_idx]); // now tmp holds next element
      }
      //std::sort(_nodes[i].first(), _nodes[i+1].first(),
      //      [] (ARC& a, ARC& b) { return a.head() < b.head(); });
      // correct sister arcs
      for_each_arc(i, [](ARC* a) { a->sister()->set_sister(a); });
      //for(ARC* a = _nodes[i].first(); a<_nodes[i+1].first() ; ++a) {
      //   a->sister()->set_sister(a);
      //}
   }
   // arcs are ordered by now!
   
	// testing network for possible excess overflow
	for ( NODE *ndp = _nodes.get() + _node_min; ndp <= _nodes.get() + _node_max; ndp ++ ) {
		cap_in  =   ( ndp->excess() );
		cap_out = - ( ndp->excess() );
		for ( _arc_current = ndp->first(); _arc_current != (ndp+1)->first(); 
			  _arc_current ++ ) {
			arc_num = N_ARC(_arc_current);
			if ( _cap[arc_num] > 0 ) cap_out += _cap[arc_num];
			if ( _cap[arc_num] == 0 ) 
				cap_in += _cap[ N_ARC(_arc_current->sister()) ];
		}
	}
	if ( _node_min != 0 ) {
      throw std::runtime_error("node ids must start from 0 inside CS2");
	}

	// adjustments due to nodes' ids being between _node_min - _node_max;
	_n = _node_max - _node_min + 1;
	//_nodes = _nodes.get() + _node_min;

	// () free internal memory, not needed anymore inside CS2;
   _arc_first.reset();
   _arc_tail.reset();
	//free ( _arc_first ); 
	//free ( _arc_tail );
}

template<typename COST_TYPE, typename CAPACITY_TYPE>
void MCMF_CS2<COST_TYPE,CAPACITY_TYPE>::cs2_initialize()
{
	// initialization; 
	// called after allocate_arrays() and all nodes and arcs have been inputed;

	NODE *i; // current node
	ARC *a; // current arc
	ARC *a_stop;
	BUCKET *b; // current bucket
	long df;

	_f_scale = (long) SCALE_DEFAULT;
	_sentinel_node = _nodes.get() + _n;
	_sentinel_arc  = _arcs.get() + 2*_m;

   for_each_node( [this](NODE* i) {
		i->set_price( 0);
		i->set_suspended( i->first());
		i->set_q_next( _sentinel_node);
         });

	_sentinel_node->set_first( _sentinel_arc);
	_sentinel_node->set_suspended( _sentinel_arc);

	// saturate negative arcs, e.g. in the circulation problem case
	for ( i = _nodes.get(); i != _sentinel_node; i ++ ) {
		for ( a = i->first(), a_stop = (i + 1)->suspended(); a != a_stop; a ++ ) {
			if ( a->cost() < 0) {
				if ( ( df = a->rez_capacity()) > 0) {
					increase_flow( i, a->head(), a, df);
				}
			}
		}
	}

	_dn = _n + 1;
	//if ( _no_zero_cycles == true) { // NO_ZERO_CYCLES
	//	_dn = 2 * _dn;
	//}

	for ( a = _arcs.get(); a != _sentinel_arc; a ++ ) {
		a->multiply_cost( _dn);
	}

	//if ( _no_zero_cycles == true) { // NO_ZERO_CYCLES
	//	for ( a = _arcs.get(); a != _sentinel_arc; a ++ ) {
	//		if ((a->cost() == 0) && (a->sister()->cost() == 0)) {
	//			a->set_cost( 1);
	//			a->sister()->set_cost( -1);
	//		}
	//	}
	//}

	if ((double) _max_cost * (double) _dn > MAX_64) {
      std::cout << "Warning:  Arc lengths too large, overflow possible\n";
	}
	_mmc = _max_cost * _dn;

	_linf = (long) (_dn * std::ceil(_f_scale) + 2);

	//_buckets = (BUCKET*) calloc ( _linf, sizeof(BUCKET));
	_buckets = std::unique_ptr<BUCKET[]>( new BUCKET[_linf] );//, sizeof(BUCKET));
   if(_buckets.get() == nullptr) { throw std::runtime_error("could not allocate _buckets"); }


	//if ( _buckets == NULL )
   //   throw std::runtime_error("Allocation fault");

	_l_bucket = _buckets.get() + _linf;

	//_dnode = new NODE; // used as reference;

	for ( b = _buckets.get(); b != _l_bucket; b ++ ) {
		reset_bucket( b);
	}

	_epsilon = _mmc;
	if ( _epsilon < 1) {
		_epsilon = 1;
	}

	_price_min = -PRICE_MAX;

	_cut_off_factor = CUT_OFF_COEF * std::pow( _n, CUT_OFF_POWER);

	_cut_off_factor = std::max( _cut_off_factor, double(CUT_OFF_MIN));

	_n_ref = 0;

	_flag_price = 0;

	_dummy_node = &_d_node;

	_excq_first = NULL;

	//print_graph(); // debug;
}

template<typename COST_TYPE, typename CAPACITY_TYPE>
void MCMF_CS2<COST_TYPE,CAPACITY_TYPE>::up_node_scan( NODE *i)
{
	NODE *j; // opposite node
	ARC *a; // (i, j)
	ARC *a_stop; // first arc from the next node
	ARC *ra; // (j, i)
	BUCKET *b_old; // old bucket contained j
	BUCKET *b_new; // new bucket for j
	long i_rank;
	long j_rank; // ranks of nodes
	long j_new_rank;             
	price_t rc; // reduced cost of (j, i)
	price_t dr; // rank difference

	i_rank = i->rank();

	// scanning arcs;
	for ( a = i->first(), a_stop = (i + 1)->suspended(); a != a_stop; a ++ ) {

		ra = a->sister();

		if ( OPEN ( ra ) ) {
			j = a->head();
			j_rank = j->rank();

			if ( j_rank > i_rank ) {
				if ( ( rc = REDUCED_COST( j, i, ra ) ) < 0 ) {
					j_new_rank = i_rank;
				} else {
					dr = rc / _epsilon;
					j_new_rank = ( dr < _linf ) ? i_rank + (long)dr + 1 : _linf;
				}

				if ( j_rank > j_new_rank ) {
					j->set_rank( j_new_rank);
					j->set_current( ra);

					if ( j_rank < _linf ) {
						b_old = _buckets.get() + j_rank;
						REMOVE_FROM_BUCKET( j, b_old );
					}

					b_new = _buckets.get() + j_new_rank;
					insert_to_bucket( j, b_new );
				}
			}
		}
	}

	i->dec_price( i_rank * _epsilon);
	i->set_rank( -1);
}

template<typename COST_TYPE, typename CAPACITY_TYPE>
void MCMF_CS2<COST_TYPE,CAPACITY_TYPE>::price_update()
{
	NODE *i;
	excess_t remain;
	// total excess of unscanned nodes with positive excess;
	BUCKET *b; // current bucket;
	price_t dp; // amount to be subtracted from prices;

	for ( i = _nodes.get(); i != _sentinel_node; i ++ ) {
		if ( i->excess() < 0 ) {
			insert_to_bucket( i, _buckets.get() );
			i->set_rank( 0);
		} else {
			i->set_rank( _linf);
		}
	}

	remain = _total_excess;
	if ( remain < 0.5 ) return;

	// scanning buckets, main loop;
	for ( b = _buckets.get(); b != _l_bucket; b ++ ) {

		while ( nonempty_bucket( b) ) {

			GET_FROM_BUCKET( i, b );
			up_node_scan( i );

			if ( i ->excess() > 0 ) {
				remain -= ( i->excess());
				if ( remain <= 0 ) break; 
			}
		}
		if ( remain <= 0 ) break; 
	}

	if ( remain > 0.5 ) _flag_updt = 1;

	// finishup 
	// changing prices for nodes which were not scanned during main loop;
	dp = ( b - _buckets.get() ) * _epsilon;

	for ( i = _nodes.get(); i != _sentinel_node; i ++ ) {

		if ( i->rank() >= 0 ) {
			if ( i->rank() < _linf ) {
				REMOVE_FROM_BUCKET( i, ( _buckets.get() + i->rank()) );
			}
			if ( i->price() > _price_min ) {
				i->dec_price( dp);
			}			
		}
	}
}

template<typename COST_TYPE, typename CAPACITY_TYPE>
int MCMF_CS2<COST_TYPE,CAPACITY_TYPE>::relabel( NODE *i)
{
	ARC *a; // current arc from i
	ARC *a_stop; // first arc from the next node
	ARC *a_max; // arc which provides maximum price
	price_t p_max; // current maximal price
	price_t i_price; // price of node  i
	price_t dp; // current arc partial residual cost

	p_max = _price_min;
	i_price = i->price();

	a_max = NULL;

	// 1/2 arcs are scanned;
	for ( a = i->current() + 1, a_stop = (i + 1)->suspended(); a != a_stop; a ++ ) {
		
		if ( OPEN(a) && ( (dp = (a->head()->price() - a->cost())) > p_max ) ) {
			if ( i_price < dp ) {
				i->set_current( a);
				return ( 1);
			}
			p_max = dp;
			a_max = a;
		}
	}

	// 2/2 arcs are scanned;
	for ( a = i->first(), a_stop = i->current() + 1; a != a_stop; a ++ ) {
		if ( OPEN( a) && ( (dp = (a->head()->price() - a->cost())) > p_max ) ) {
			if ( i_price < dp ) {
				i->set_current( a);
				return ( 1);
			}
			p_max = dp;
			a_max = a;
		}
	}

	// finishup
	if ( p_max != _price_min ) {
		i->set_price( p_max - _epsilon);
		i->set_current( a_max);
	} 
	else { // node can't be relabelled;
		if ( i->suspended() == i->first() ) {
			if ( i->excess() == 0 ) {
				i->set_price( _price_min);
			} else {
				if ( _n_ref == 1 ) {
               throw std::runtime_error("problem infeasible");
				} else {
               throw std::runtime_error("price overflow");
				}
			}
		} else { // node can't be relabelled because of suspended arcs;
			_flag_price = 1;
		}
	}

	_n_rel ++;
	return ( 0);
}

template<typename COST_TYPE, typename CAPACITY_TYPE>
void MCMF_CS2<COST_TYPE,CAPACITY_TYPE>::discharge( NODE *i)
{
	ARC *a;// an arc from i
	NODE *j; // head of a
	long df; // amoumt of flow to be pushed through a
	excess_t j_exc; // former excess of j

	a = i->current();
	j = a->head();

	if ( !ADMISSIBLE( i, j, a ) ) { 
		relabel( i );
		a = i->current();
		j = a->head();
	}

	while ( 1 ) {

		j_exc = j->excess();
		if ( j_exc >= 0 ) {

			df = std::min( long(i->excess()), a->rez_capacity() );
			if ( j_exc == 0) _n_src++;
			increase_flow( i, j, a, df ); // INCREASE_FLOW 

			if ( out_of_excess_q( j ) ) {
				insert_to_excess_q( j );
			}
		} 
		else { // j_exc < 0;

			df = std::min( long(i->excess()), a->rez_capacity() );
			increase_flow( i, j, a, df ); // INCREASE_FLOW 

			if ( j->excess() >= 0 ) {
				if ( j->excess() > 0 ) {
					_n_src ++;
					relabel( j );
					insert_to_excess_q( j );
				}
				_total_excess += j_exc;
			}
			else {
				_total_excess -= df;
			}
		}
  
		if ( i->excess() <= 0) _n_src --;
		if ( i->excess() <= 0 || _flag_price ) break;

		relabel( i );

		a = i->current();
		j = a->head();
	}

	i->set_current( a);
}

template<typename COST_TYPE, typename CAPACITY_TYPE>
int MCMF_CS2<COST_TYPE,CAPACITY_TYPE>::price_in()
{
	NODE *i; // current node
	NODE *j;
	ARC *a; // current arc from i
	ARC *a_stop; // first arc from the next node
	ARC *b; // arc to be exchanged with suspended
	ARC *ra; // opposite to a
	ARC *rb; // opposite to b
	price_t rc; // reduced cost
	int n_in_bad; // number of priced_in arcs with negative reduced cost
	int bad_found; // if 1 we are at the second scan if 0 we are at the first scan
	excess_t i_exc; // excess of i
	excess_t df; // an amount to increase flow


	bad_found = 0;
	n_in_bad = 0;

 restart:

	for ( i = _nodes.get(); i != _sentinel_node; i ++ ) {

		for ( a = i->first() - 1, a_stop = i->suspended() - 1; a != a_stop; a -- ) {

			rc = REDUCED_COST( i, a->head(), a );
			if ( ( rc < 0) && ( a->rez_capacity() > 0) ) { // bad case;
				if ( bad_found == 0 ) {
					bad_found = 1;
					update_cut_off();
					goto restart;
				}
				df = a->rez_capacity();
				increase_flow( i, a->head(), a, df );
	    
				ra = a->sister();
				j  = a->head();
	    
				i->dec_first();
				b = i->first();
				exchange( a, b );
	    
				if ( SUSPENDED( j, ra ) ) {
					j->dec_first();
					rb = j->first();
					exchange( ra, rb );
				}
	    
				n_in_bad ++; 
			}
			else {
				if ( ( rc < _cut_on ) && ( rc > -_cut_on ) ) {
					i->dec_first();
					b = i->first();
					exchange( a, b );
				}
			}
		}
	}


	if ( n_in_bad != 0 ) {

		_n_bad_pricein ++;

		// recalculating excess queue;
		_total_excess = 0;
		_n_src = 0;
		reset_excess_q();

		for ( i = _nodes.get(); i != _sentinel_node; i ++ ) {
			i->set_current( i->first());
			i_exc = i->excess();
			if ( i_exc > 0 ) { // i is a source;
				_total_excess += i_exc;
				_n_src ++;
				insert_to_excess_q( i );
			}
		}

		insert_to_excess_q( _dummy_node );
	}

	if ( _time_for_price_in == TIME_FOR_PRICE_IN2)
		_time_for_price_in = TIME_FOR_PRICE_IN3;
	if ( _time_for_price_in == TIME_FOR_PRICE_IN1)
		_time_for_price_in = TIME_FOR_PRICE_IN2;

	return ( n_in_bad);
}

template<typename COST_TYPE, typename CAPACITY_TYPE>
void MCMF_CS2<COST_TYPE,CAPACITY_TYPE>::refine()
{
	NODE *i; // current node
	excess_t i_exc; // excess of i
	long np, nr, ns; // variables for additional print
	int pr_in_int; // current number of updates between price_in

	//np = _n_push; 
	//nr = _n_relabel; 
	//ns = _n_scan;

	_n_ref ++;
	_n_rel = 0;
	pr_in_int = 0;

	// initialize;
	_total_excess = 0;
	_n_src = 0;
	reset_excess_q();

	_time_for_price_in = TIME_FOR_PRICE_IN1;

	for ( i = _nodes.get(); i != _sentinel_node; i ++ ) {
		i->set_current( i->first());
		i_exc = i->excess();
		if ( i_exc > 0 ) { // i  is a source 
			_total_excess += i_exc;
			_n_src++;
			insert_to_excess_q( i );
		}
	}

	if ( _total_excess <= 0 ) return;

	// (2) main loop

	while ( 1 ) {

		if ( empty_excess_q() ) {
			if ( _n_ref > PRICE_OUT_START ) {
				pr_in_int = 0;
				price_in();
			}
	  
			if ( empty_excess_q() ) break;
		}
		
		REMOVE_FROM_EXCESS_Q( i );

		// push all excess out of i
		if ( i->excess() > 0 ) {
			discharge( i );

			if ( time_for_update() || _flag_price ) {
				if ( i->excess() > 0 ) {
					insert_to_excess_q( i );
				}

				if ( _flag_price && ( _n_ref > PRICE_OUT_START ) ) {
					pr_in_int = 0;
					price_in();
					_flag_price = 0;
				}

				price_update();

				while ( _flag_updt ) {
					if ( _n_ref == 1 ) {
                  throw std::runtime_error("problem infeasible");
					} else {
						_flag_updt = 0;
						update_cut_off();
						_n_bad_relabel ++;
						pr_in_int = 0;
						price_in();
						price_update();
					}
				}
				_n_rel = 0;

				if ( _n_ref > PRICE_OUT_START && (pr_in_int ++ > _time_for_price_in) ) {
					pr_in_int = 0;
					price_in();
				}
			}
		}
	}

	return;
}

template<typename COST_TYPE, typename CAPACITY_TYPE>
int MCMF_CS2<COST_TYPE,CAPACITY_TYPE>::price_refine()
{
	NODE *i; // current node
	NODE *j; // opposite node
	NODE *ir; // nodes for passing over the negative cycle
	NODE *is;
	ARC *a; // arc (i,j)
	ARC *a_stop; // first arc from the next node
	ARC *ar;
	long bmax;            // number of farest nonempty bucket
	long i_rank;          // rank of node i
	long j_rank;         // rank of node j
	long j_new_rank;      // new rank of node j
	BUCKET *b;              // current bucket
	BUCKET *b_old;          // old and new buckets of current node
	BUCKET *b_new;
	price_t rc = 0; // reduced cost of a
	price_t dr; // ranks difference
	price_t dp;
	int cc;              
	// return code: 1 - flow is epsilon optimal
	// 0 - refine is needed       
	long df; // cycle capacity
	int nnc; // number of negative cycles cancelled during one iteration
	int snc; // total number of negative cycle cancelled

	cc = 1;
	snc = 0;

	_snc_max = ( _n_ref >= START_CYCLE_CANCEL) ? MAX_CYCLES_CANCELLED : 0;


	// (1) main loop
	// while negative cycle is found or eps-optimal solution is constructed
	while ( 1 ) { 

		nnc = 0;
		for ( i = _nodes.get(); i != _sentinel_node; i ++ ) {
			i->set_rank( 0);
			i->set_inp( WHITE);
			i->set_current( i->first());
		}
		reset_stackq();

		for ( i = _nodes.get(); i != _sentinel_node; i ++ ) {
			if ( i->inp() == BLACK ) continue;

			i->set_b_next( NULL);

			// deapth first search 
			while ( 1 ) {
				i->set_inp( GREY);

				// scanning arcs from node i starting from current
				for ( a = i->current(), a_stop = (i + 1)->suspended(); a != a_stop; a ++) {
					if ( OPEN( a ) ) {
						j = a->head();
						if ( REDUCED_COST ( i, j, a ) < 0 ) {
							if ( j->inp() == WHITE ) { // fresh node  - step forward 
								i->set_current( a);
								j->set_b_next( i);
								i = j;
								a = j->current();
								a_stop = (j+1)->suspended();
								break;
							}

							if ( j->inp() == GREY ) { // cycle detected 
								cc = 0;
								nnc ++;
								i->set_current( a);
								is = ir = i;
								df = MAX_32;

								while ( 1 ) {
									ar = ir->current();
									if ( ar->rez_capacity() <= df ) {
										df = ar->rez_capacity();
										is = ir;
									}
									if ( ir == j ) break;
									ir = ir->b_next();
								}

								ir = i;

								while ( 1 ) {
									ar = ir->current();
									increase_flow( ir, ar->head(), ar, df);
									if ( ir == j ) break;
									ir = ir->b_next();
								}

								if ( is != i ) {
									for ( ir = i; ir != is; ir = ir->b_next() ) {
										ir->set_inp( WHITE);
									}
									i = is;
									a = is->current() + 1;
									a_stop = (is+1)->suspended();
									break;
								}
							}
						}
						// if j-color is BLACK - continue search from i
					}
				} // all arcs from i are scanned 

				if ( a == a_stop ) {
					// step back 
					i->set_inp( BLACK);
					j = i->b_next();
					stackq_push( i );
					if ( j == NULL ) break;
					i = j;
					i->inc_current();
				}

			} // end of deapth first search
		} // all nodes are scanned


		// () no negative cycle
		// computing longest paths with eps-precision

		snc += nnc;
		if ( snc < _snc_max ) cc = 1;
		if ( cc == 0 ) break;
		bmax = 0;

		while ( nonempty_stackq() ) {

			STACKQ_POP( i );
			i_rank = i->rank();
			for ( a = i->first(), a_stop = (i + 1)->suspended(); a != a_stop; a ++) {

				if ( OPEN( a ) ) {
					j  = a->head();
					rc = REDUCED_COST( i, j, a );

					if ( rc < 0 ) { // admissible arc;
						dr = (price_t) (( - rc - 0.5 ) / _epsilon);
						if (( j_rank = dr + i_rank ) < _linf ) {
							if ( j_rank > j->rank() )
								j->set_rank( j_rank);
						}
					}
				}
			} // all arcs from i are scanned 

			if ( i_rank > 0 ) {
				if ( i_rank > bmax ) bmax = i_rank;
				b = _buckets.get() + i_rank;
				insert_to_bucket( i, b );
			}
		} // end of while-cycle: all nodes are scanned - longest distancess are computed;


		if ( bmax == 0 ) // preflow is eps-optimal;
			{ break; }


		for ( b = _buckets.get() + bmax; b != _buckets.get(); b -- ) {
			i_rank = b - _buckets.get();
			dp = i_rank * _epsilon;

			while ( nonempty_bucket( b) ) {
				GET_FROM_BUCKET( i, b );

				for ( a = i->first(), a_stop = (i + 1)->suspended(); a != a_stop; a ++) {
					if ( OPEN( a ) ) {
						j = a->head();
						j_rank = j->rank();
						if ( j_rank < i_rank ) {
							rc = REDUCED_COST( i, j, a );
							if ( rc < 0 ) {
								j_new_rank = i_rank;
							} else {
								dr = rc / _epsilon;
								j_new_rank = ( dr < _linf ) ? i_rank - ( (long)dr + 1 ) : 0;
							}
							if ( j_rank < j_new_rank ) {
								if ( cc == 1 ) {
									j->set_rank( j_new_rank);
									if ( j_rank > 0 ) {
										b_old = _buckets.get() + j_rank;
										REMOVE_FROM_BUCKET( j, b_old );
									}
									b_new = _buckets.get() + j_new_rank;
									insert_to_bucket( j, b_new );
								}
								else {
									df = a->rez_capacity();
									increase_flow( i, j, a, df );
								}
							}
						}
					} // end if opened arc 
				} // all arcs are scanned

				i->dec_price( dp);

			} // end of while-cycle: the bucket is scanned 
		} // end of for-cycle: all buckets are scanned 

		if ( cc == 0 ) break;

	} // end of main loop



	// (2) finish
	// if refine needed - saturate non-epsilon-optimal arcs;

	if ( cc == 0 ) { 
		for ( i = _nodes.get(); i != _sentinel_node; i ++) {
			for ( a = i->first(), a_stop = (i + 1)->suspended(); a != a_stop; a ++) {
				if ( REDUCED_COST( i, a->head(), a ) < - _epsilon ) {
					if ( ( df = a->rez_capacity() ) > 0 ) {
						increase_flow( i, a->head(), a, df );
					}
				}
			}
		}
	}

	return ( cc );
}

template<typename COST_TYPE, typename CAPACITY_TYPE>
void MCMF_CS2<COST_TYPE,CAPACITY_TYPE>::compute_prices()
{
   NODE *i; // current node
   NODE *j; // opposite node
   ARC *a; // arc (i,j)
   ARC *a_stop; // first arc from the next node
   long bmax; // number of farest nonempty bucket
   long i_rank; // rank of node i
   long j_rank; // rank of node j
   long j_new_rank; // new rank of node j
   BUCKET *b; // current bucket
	BUCKET *b_old; // old and new buckets of current node
   BUCKET *b_new;
   price_t rc; // reduced cost of a
   price_t dr; // ranks difference
   price_t dp;
   int cc; // return code: 1 - flow is epsilon optimal 0 - refine is needed

   cc = 1;

   // (1) main loop
   // while negative cycle is found or eps-optimal solution is constructed
   while ( 1 ) {

      for ( i = _nodes.get(); i != _sentinel_node; i ++) {
         i->set_rank( 0);
         i->set_inp( WHITE);
         i->set_current( i->first());
      }
      reset_stackq();

      for ( i = _nodes.get(); i != _sentinel_node; i ++ ) {
         if ( i->inp() == BLACK ) continue;

         i->set_b_next( NULL);
         // depth first search
         while ( 1 ) {
            i->set_inp( GREY);

            // scanning arcs from node i
            for ( a = i->suspended(), a_stop = (i + 1)->suspended(); a != a_stop; a ++) {
               if ( OPEN( a ) ) {
                  j = a->head();
                  if ( REDUCED_COST( i, j, a ) < 0 ) {
                     if ( j->inp() == WHITE ) { // fresh node  - step forward
                        i->set_current( a);
                        j->set_b_next( i);
                        i = j;
                        a = j->current();
                        a_stop = (j+1)->suspended();
                        break;
                     }

                     if ( j->inp() == GREY ) { // cycle detected; should not happen
                        cc = 0;
                     }                     
                  }
                  // if j-color is BLACK - continue search from i
               } 
            } // all arcs from i are scanned 

            if ( a == a_stop ) {
               // step back
               i->set_inp( BLACK);
               j = i->b_next();
               stackq_push( i );
               if ( j == NULL ) break;
               i = j;
               i->inc_current();
            }

         } // end of deapth first search
      } // all nodes are scanned


      // no negative cycle
      // computing longest paths

      if ( cc == 0 ) break;
      bmax = 0;

      while ( nonempty_stackq() ) {
         STACKQ_POP( i );
         i_rank = i->rank();
         for ( a = i->suspended(), a_stop = (i + 1)->suspended(); a != a_stop; a ++) {
            if ( OPEN( a ) ) {
               j  = a->head();
               rc = REDUCED_COST( i, j, a );


               if ( rc < 0 ) {// admissible arc
                  dr = - rc;
                  if (( j_rank = dr + i_rank ) < _linf ) {
                     if ( j_rank > j->rank() )
                        j->set_rank( j_rank);
                  }
               }
            }
         } // all arcs from i are scanned

         if ( i_rank > 0 ) {
            if ( i_rank > bmax ) bmax = i_rank;
            b = _buckets.get() + i_rank;
            insert_to_bucket( i, b );
         }
      } // end of while-cycle: all nodes are scanned - longest distancess are computed;

      if ( bmax == 0 )
      { break; }

      for ( b = _buckets.get() + bmax; b != _buckets.get(); b -- ) {
         i_rank = b - _buckets.get();
         dp = i_rank;

         while ( nonempty_bucket( b) ) {
            GET_FROM_BUCKET( i, b );

            for ( a = i->suspended(), a_stop = (i + 1)->suspended(); a != a_stop; a ++) {
               if ( OPEN( a ) ) {
                  j = a->head();
                  j_rank = j->rank();
                  if ( j_rank < i_rank ) {
                     rc = REDUCED_COST( i, j, a );

                     if ( rc < 0 ) {
                        j_new_rank = i_rank;
                     } else {
                        dr = rc;
                        j_new_rank = ( dr < _linf ) ? i_rank - ( (long)dr + 1 ) : 0;
                     }
                     if ( j_rank < j_new_rank ) {
                        if ( cc == 1 ) {
                           j->set_rank( j_new_rank);
                           if ( j_rank > 0 ) {
                              b_old = _buckets.get() + j_rank;
                              REMOVE_FROM_BUCKET( j, b_old );
                           }
                           b_new = _buckets.get() + j_new_rank;
                           insert_to_bucket( j, b_new );
                        }
                     }
                  }
               } // end if opened arc
            } // all arcs are scanned

            i->dec_price( dp);

         } // end of while-cycle: the bucket is scanned
      } // end of for-cycle: all buckets are scanned

      if ( cc == 0 ) break;

   } // end of main loop
} 

template<typename COST_TYPE, typename CAPACITY_TYPE>
void MCMF_CS2<COST_TYPE,CAPACITY_TYPE>::price_out()
{
	NODE *i; // current node
	ARC *a; // current arc from i 
	ARC *a_stop; // first arc from the next node 
	ARC *b; // arc to be exchanged with suspended
	double n_cut_off; // -cut_off
	double rc; // reduced cost

	n_cut_off = - _cut_off;

	for ( i = _nodes.get(); i != _sentinel_node; i ++) {
		for ( a = i->first(), a_stop = (i + 1)->suspended(); a != a_stop; a ++) {

			rc = REDUCED_COST( i, a->head(), a );		
			if ( ( rc > _cut_off && CLOSED(a->sister()) ) ||
				 ( rc < n_cut_off && CLOSED(a) ) ) { // suspend the arc

				b = i->first();
				i->inc_first();
				exchange( a, b );
			}
		}
	}
}

template<typename COST_TYPE, typename CAPACITY_TYPE>
int MCMF_CS2<COST_TYPE,CAPACITY_TYPE>::update_epsilon()
{
	// decrease epsilon after epsilon-optimal flow is constructed;
	if ( _epsilon <= 1 ) return ( 1 );

	_epsilon = (price_t) (std::ceil ( _epsilon / _f_scale ));
	_cut_off = _cut_off_factor * _epsilon;
	_cut_on = _cut_off * CUT_OFF_GAP;

	return ( 0 );
}




template<typename COST_TYPE, typename CAPACITY_TYPE>
void MCMF_CS2<COST_TYPE,CAPACITY_TYPE>::init_solution()
{
	ARC *a; // current arc (i,j)
	NODE *i; // tail of a
	NODE *j; // head of a
	long df; // residual capacity

	for ( a = _arcs.get(); a != _sentinel_arc; a ++ ) {
		if ( a->rez_capacity() > 0 && a->cost() < 0 ) {
			df = a->rez_capacity();
			i  = a->sister()->head();
			j  = a->head();
			increase_flow( i, j, a, df );
		}
	}
}

template<typename COST_TYPE, typename CAPACITY_TYPE>
void MCMF_CS2<COST_TYPE,CAPACITY_TYPE>::cs_cost_reinit()
{
	//if ( _cost_restart == false)
	//	return;
	
	NODE *i; // current node
	ARC *a;          // current arc
	ARC *a_stop;
	BUCKET *b; // current bucket
	price_t rc, minc, sum;


	for ( b = _buckets.get(); b != _l_bucket; b ++) {
		reset_bucket( b);
	}

	rc = 0;
	for ( i = _nodes.get(); i != _sentinel_node; i ++) {
		rc = std::min(rc, i->price());
		i->set_first( i->suspended());
		i->set_current( i->first());
		i->set_q_next( _sentinel_node);
	}

	// make prices nonnegative and multiply 
	for ( i = _nodes.get(); i != _sentinel_node; i ++) {
		i->set_price( (i->price() - rc) * _dn);
	}

	// multiply arc costs
	for (a = _arcs.get(); a != _sentinel_arc; a ++) {
		a->multiply_cost( _dn);
	}

	sum = 0;
	for ( i = _nodes.get(); i != _sentinel_node; i ++) {
		minc = 0;
		for ( a = i->first(), a_stop = (i + 1)->suspended(); a != a_stop; a ++) {		
			if ( (OPEN(a) && ((rc = REDUCED_COST(i, a->head(), a)) < 0)) )
				minc = std::max( _epsilon, -rc);
		}
		sum += minc;
	}

	_epsilon = std::ceil(sum / _dn);

	_cut_off_factor = CUT_OFF_COEF * std::pow(_n, CUT_OFF_POWER);

	_cut_off_factor = std::max( _cut_off_factor, double(CUT_OFF_MIN));

	_n_ref = 0;

   _n_bad_pricein = _n_bad_relabel = 0;

	_flag_price = 0;

	_excq_first = NULL;
}

template<typename COST_TYPE, typename CAPACITY_TYPE>
long long int MCMF_CS2<COST_TYPE,CAPACITY_TYPE>::cs2_cost_restart()
{
	// restart after a cost update;
	//if ( _cost_restart == false)
	//	return;

	int cc; // for storing return code;
   double objective_cost;

   std::cout << "c \nc ******************************\n";
   std::cout << "c Restarting after a cost update\n";
   std::cout << "c ******************************\nc\n";

	cs_cost_reinit();
  
   std::cout << "c Init. epsilon = %6.0f" << _epsilon << "\n";
	cc = update_epsilon();
  
	if (cc != 0) {
      std::cout << "c Old solution is optimal\n";
	} 
	else {
		do { // scaling loop
			while ( 1 ) {
				if ( ! price_refine() ) 
					break;

				if ( _n_ref >= PRICE_OUT_START ) {
					if ( price_in() ) 
						break;
				}
				if ((cc = update_epsilon ())) 
					break;
			}
			if (cc) break;
			refine();
			if ( _n_ref >= PRICE_OUT_START ) {
				price_out();
			}
			if ( update_epsilon() )
				break;
		} while ( cc == 0 );
	}

	finishup( &objective_cost );
   return objective_cost;
}

template<typename COST_TYPE, typename CAPACITY_TYPE>
long long int MCMF_CS2<COST_TYPE,CAPACITY_TYPE>::compute_objective_cost() const
{
   double obj = 0;
   int na;
   ARC* a;
   long flow;
   for (a = _arcs.get(), na = 0; a != _sentinel_arc ; a ++, na ++ ) {
      double cs = a->cost();
      if ( _cap[na]  > 0 && (flow = _cap[na] - a->rez_capacity()) != 0 )
         obj += (double) cs * (double) flow;
   }

   return obj;
}



template<typename COST_TYPE, typename CAPACITY_TYPE>
void MCMF_CS2<COST_TYPE,CAPACITY_TYPE>::finishup( double *objective_cost)
{
	ARC *a; // current arc
	long na; // corresponding position in capacity array
	double obj_internal = 0; // objective
	price_t cs; // actual arc cost
	long flow; // flow through an arc
	NODE *i;

	// (1) NO_ZERO_CYCLES?
	//if ( _no_zero_cycles == true) {
	//	for ( a = _arcs.get(); a != _sentinel_arc; a ++ ) {
	//		if ( a->cost() == 1) {
	//			assert( a->sister()->cost() == -1);
	//			a->set_cost( 0);
	//			a->sister()->set_cost( 0);
	//		}
	//	}
	//}

	// (2)
	for ( a = _arcs.get(), na = 0; a != _sentinel_arc ; a ++, na ++ ) {
		cs = a->cost() / _dn;
		if ( _cap[na]  > 0 && (flow = _cap[na] - a->rez_capacity()) != 0 )
			obj_internal += (double) cs * (double) flow;
		a->set_cost( cs);
	}

	for ( i = _nodes.get(); i != _sentinel_node; i ++) {
		i->set_price( (i->price() / _dn));
	}

	// (3) COMP_DUALS?
	//if ( _comp_duals == true) {
	//	compute_prices(); // sometimes goes into infinite loop
	//}

	*objective_cost = obj_internal;
}

template<typename COST_TYPE, typename CAPACITY_TYPE>
void MCMF_CS2<COST_TYPE,CAPACITY_TYPE>::cs2( double *objective_cost)
{
	// the main calling function;
	int cc = 0; // for storing return code;

  
	// (1) update epsilon first;
	update_epsilon();


	// (2) scaling loop;
	do {
		refine();

		if ( _n_ref >= PRICE_OUT_START )
			price_out();

		if ( update_epsilon() ) 
			break;

		while (1) {
			if ( ! price_refine() ) 
				break;

			if ( _n_ref >= PRICE_OUT_START ) {
				if ( price_in() ) break; 
				if ( (cc = update_epsilon()) ) break;
			}
		}
	} while ( cc == 0 );


	// (3) finishup;
	finishup( objective_cost );
}

template<typename COST_TYPE, typename CAPACITY_TYPE>
long long int MCMF_CS2<COST_TYPE,CAPACITY_TYPE>::run_cs2()
{
	// in order to solve the flow problem we have to follow these steps:
	// 1. ctor of MCMF_CS2 // sets num of nodes and arcs
	//                     // it also calls allocate_arrays()
	// 2. call set_arc() for each arc
	// 3. call set_supply_demand_of_node() for non-transhipment nodes
	// 4. pre_processing()
	// 5. cs2_initialize()
	// 6. cs2()
	// 7. retreive results
	//
	// this function is basically a wrapper to implement steps 4, 5, 6;

	double objective_cost;


	// (4) ordering, etc.;
	pre_processing();


	// () CHECK_SOLUTION?
	//if ( _check_solution == true) {
	//	_node_balance = (long long int *) calloc (_n+1, sizeof(long long int));
	//	for ( NODE *i = _nodes.get(); i < _nodes.get() + _n; i ++ ) {
	//		_node_balance[N_NODE(i)] = i->excess();
	//	}
	//}


	// (5) initializations;
	//_m = 2 * _m;
	cs2_initialize(); // works already with 2*m;
	
	// (6) run CS2;
	cs2( &objective_cost );

   return objective_cost;
}

template<typename COST_TYPE, typename CAPACITY_TYPE>
long long int MCMF_CS2_STAT<COST_TYPE,CAPACITY_TYPE>::run_cs2()
{
	//print_graph(); // exit(1); // debug;

   std::cout << "\nc CS 4.3\n";
   std::cout << "c nodes: " << this->_n << " arcs: " << this->_m << "\n";
   std::cout << "c scale-factor: " << this->_f_scale << " cut-off-factor: " << this->_cut_off_factor << "\nc\n";

   auto start = std::chrono::steady_clock::now();
   long long int obj = MCMF_CS2<COST_TYPE,CAPACITY_TYPE>::run_cs2();
   auto duration = std::chrono::duration_cast<std::chrono::milliseconds> 
                            (std::chrono::steady_clock::now() - start);

   std::cout << "c time:         " << duration.count() << "    cost:       " << obj << "\n";
   //std::cout << "c refines:      " << _n_refine << "     discharges: " << _n_discharge << "\n";
   //std::cout << "c pushes:       " << _n_push << "     relabels:   " << _n_relabel << "\n";
   //std::cout << "c updates:      " << _n_update << "     u-scans:    " << _n_scan << "\n";
   //std::cout << "c p-refines:    " << _n_prefine << "     r-scans:    " << _n_prscan << "\n";
   //std::cout << "c dfs-scans:    " << _n_prscan1 << "     bad-in:     " << _n_bad_pricein << "  + " << _n_bad_relabel << "\n";
  
	// () CHECK_SOLUTION?
   std::cout << "c checking feasibility...\n"; 
   if ( check_feas() )
      std::cout << "c ...OK\n";
   else
      std::cout << "c ERROR: solution infeasible\n";
   std::cout << "c computing prices and checking CS...\n";
   //compute_prices();
   if ( check_cs() )
      std::cout << "c ...OK\n";
   else
      std::cout << "ERROR: CS violation\n";

	// () PRINT_ANS?
   //print_solution();

   return obj;
}

template<typename COST_TYPE, typename CAPACITY_TYPE>
void MCMF_CS2_STAT<COST_TYPE,CAPACITY_TYPE>::print_solution()
{
   for(auto* a=this->_arcs.get(); a<this->_arcs.get()+this->_m*2; ++a) {
      std::cout << this->N_NODE(a->sister()->head()) << "->" << this->N_NODE(a->head()) << ": flow = " << a->rez_capacity()  - this->_cap[this->N_ARC(a)] << "\n";
   }
   return;
	
   /*
	price_t cost;
   
	long ni;
	for (auto* i = _nodes.get(); i < _nodes.get() + _n; i ++ ) {
		ni = N_NODE( i );
		for (auto* a = i->suspended(); a != (i+1)->suspended(); a ++) {
         std::cout << a->rez_capacity();
			//if ( _cap[ N_ARC (a) ]  > 0 ) {
			//	printf("f %7ld %7ld %10ld\n", 
			//		   ni, N_NODE(a->head()), _cap[ N_ARC(a) ] - a->rez_capacity());
			//}
		}
    }

	// COMP_DUALS?
   cost = MAX_32;
   for (auto* i = _nodes.get(); i != _sentinel_node; i ++) {
      cost = std::min(cost, i->price());
   }
   for (auto* i = _nodes.get(); i != _sentinel_node; i ++) {
      std::cout << "p " << N_NODE(i) << " " << i->price() - cost << "\n";
   }

   std::cout << "c\n";
   */
}

template<typename COST_TYPE, typename CAPACITY_TYPE>
void MCMF_CS2_STAT<COST_TYPE,CAPACITY_TYPE>::print_graph()
{
	long ni, na;
   std::cout << "\nGraph: " << this->_n << "\n";
	for (auto* i = this->_nodes.get(); i < this->_nodes.get() + this->_n; i ++ ) {
		ni = this->N_NODE( i );
      std::cout << "\nNode " << ni << "\n";
		for (auto* a = i->suspended(); a != (i+1)->suspended(); a ++) {
			na = this->N_ARC( a );
         std::cout << "\n {" << na << "} " << ni << " -> " <<  this->N_NODE(a->head()) << " cap: " << this->_cap[this->N_ARC(a)] << " cost: " << a->cost() << "\n";
      }
    }
}

template<typename COST_TYPE, typename CAPACITY_TYPE>
int MCMF_CS2_STAT<COST_TYPE,CAPACITY_TYPE>::check_feas()
{
   // note: we must add _node_balance here as well to be able to check feasibility
   return 1;
	
   /*
	NODE *i;
	ARC *a, *a_stop;
	long fa;
	int ans = 1;

	for ( i = _nodes; i != _sentinel_node; i ++) {
		for ( a = i->suspended(), a_stop = (i + 1)->suspended(); a != a_stop; a ++) {
			if ( _cap[ N_ARC(a) ] > 0) {
				fa = _cap[ N_ARC(a) ] - a->rez_capacity();
				if ( fa < 0) {
					ans = 0;
					break;
				}
				_node_balance[ i - _nodes ] -= fa;
				_node_balance[ a->head() - _nodes ] += fa;
			}
		}
	}

	for ( i = _nodes; i != _sentinel_node; i ++) {
		if ( _node_balance[ i - _nodes ] != 0) {
			ans = 0;
			break;
		}
	}

	return ( ans);
   */
}

template<typename COST_TYPE, typename CAPACITY_TYPE>
int MCMF_CS2_STAT<COST_TYPE,CAPACITY_TYPE>::check_cs()
{
	// check complementary slackness;
	for (auto* i = this->_nodes.get(); i != this->_sentinel_node; i ++) {
      auto* a = i->suspended();
      auto* a_stop = (i+1)->suspended();
		for ( ; a != a_stop; a ++) {

			if ( OPEN(a) && (REDUCED_COST(i, a->head(), a) < 0) ) {
				return ( 0);
			}
		}
	}
	return(1);
}

template<typename COST_TYPE, typename CAPACITY_TYPE>
int MCMF_CS2_STAT<COST_TYPE,CAPACITY_TYPE>::check_eps_opt()
{
	for (auto* i = this->_nodes.get(); i != this->_sentinel_node; i ++) {
      auto* a = i->suspended();
      auto* a_stop = (i+1)->suspended();
		for ( ; a != a_stop; a ++) {
			if ( OPEN(a) && (REDUCED_COST(i, a->head(), a) < - this->_epsilon) ) {
				return ( 0);
			}
		}
	}
	return(1);
}


} // end namespace CS2_CPP

#endif
