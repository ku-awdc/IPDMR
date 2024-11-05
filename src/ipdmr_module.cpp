/*

  DO NOT MODIFY THIS FILE DIRECTLY!
  See module/ipdmr_module_template.cpp

*/ 

#include <Rcpp.h>

#include "SEIRmodel.h"

template <class RcppModuleClassName>
RcppModuleClassName* invalidate_default_constructor() {
  Rcpp::stop("Default constructor is disabled for this class");
  return 0;
}
#define DISABLE_DEFAULT_CONSTRUCTOR() .factory(invalidate_default_constructor)

RCPP_MODULE(IPDMRmodule){

	using namespace Rcpp;

	using SEIRdetN = SEIRmodel<update_type::deterministic, false, 0L, true>;
	class_<SEIRdetN>("SEIRdetN")
	  .constructor<const int, const Rcpp::CharacterVector>("Constructor")
	  .factory(invalidate_default_constructor)
	    
		.method("show", &SEIRdetN::show, "The show/print method")
	  .method("update", &SEIRdetN::update, "The update method")  
	  .method("run", &SEIRdetN::run, "The run method")  
	  .method("save", &SEIRdetN::save, "The save method")  
	  .method("reset", &SEIRdetN::reset, "The reset method")  
	    
		.property("state", &SEIRdetN::get_state, "Get current state")
	  .property("S", &SEIRdetN::get_S, &SEIRdetN::set_S, "Get/set current S")
	  .property("E", &SEIRdetN::get_E, &SEIRdetN::set_E, "Get/set current E")
	  .property("I", &SEIRdetN::get_I, &SEIRdetN::set_I, "Get/set current I")
	  .property("R", &SEIRdetN::get_R, &SEIRdetN::set_R, "Get/set current R")
	  .property("N", &SEIRdetN::get_N, "Get current N")
	  .property("beta", &SEIRdetN::get_beta, &SEIRdetN::set_beta, "Get/set current beta")
	  .property("omega", &SEIRdetN::get_omega, &SEIRdetN::set_omega, "Get/set current omega")
	  .property("gamma", &SEIRdetN::get_gamma, &SEIRdetN::set_gamma, "Get/set current gamma")
	  .property("delta", &SEIRdetN::get_delta, &SEIRdetN::set_delta, "Get/set current delta")
	  .property("repl", &SEIRdetN::get_repl, &SEIRdetN::set_repl, "Get/set current repl")
	  .property("cull", &SEIRdetN::get_cull, &SEIRdetN::set_cull, "Get/set current cull")
	  .property("vacc", &SEIRdetN::get_vacc, &SEIRdetN::set_vacc, "Get/set current vacc")
	  .property("time", &SEIRdetN::get_time, "Get current time")
	  .property("trans_external", &SEIRdetN::get_trans_external, &SEIRdetN::set_trans_external, "Get/set current trans_external")
	  .property("transmission_type", &SEIRdetN::get_trans_type, &SEIRdetN::set_trans_type, "Get/set current trans_type")
	; 

	using SEIRstocN = SEIRmodel<update_type::stochastic, false, 0L, true>;
	class_<SEIRstocN>("SEIRstocN")
	  .constructor<const int, const Rcpp::CharacterVector>("Constructor")
	  .factory(invalidate_default_constructor)
	    
		.method("show", &SEIRstocN::show, "The show/print method")
	  .method("update", &SEIRstocN::update, "The update method")  
	  .method("run", &SEIRstocN::run, "The run method")  
	  .method("save", &SEIRstocN::save, "The save method")  
	  .method("reset", &SEIRstocN::reset, "The reset method")  
	    
		.property("state", &SEIRstocN::get_state, "Get current state")
	  .property("S", &SEIRstocN::get_S, &SEIRstocN::set_S, "Get/set current S")
	  .property("E", &SEIRstocN::get_E, &SEIRstocN::set_E, "Get/set current E")
	  .property("I", &SEIRstocN::get_I, &SEIRstocN::set_I, "Get/set current I")
	  .property("R", &SEIRstocN::get_R, &SEIRstocN::set_R, "Get/set current R")
	  .property("N", &SEIRstocN::get_N, "Get current N")
	  .property("beta", &SEIRstocN::get_beta, &SEIRstocN::set_beta, "Get/set current beta")
	  .property("omega", &SEIRstocN::get_omega, &SEIRstocN::set_omega, "Get/set current omega")
	  .property("gamma", &SEIRstocN::get_gamma, &SEIRstocN::set_gamma, "Get/set current gamma")
	  .property("delta", &SEIRstocN::get_delta, &SEIRstocN::set_delta, "Get/set current delta")
	  .property("repl", &SEIRstocN::get_repl, &SEIRstocN::set_repl, "Get/set current repl")
	  .property("cull", &SEIRstocN::get_cull, &SEIRstocN::set_cull, "Get/set current cull")
	  .property("vacc", &SEIRstocN::get_vacc, &SEIRstocN::set_vacc, "Get/set current vacc")
	  .property("time", &SEIRstocN::get_time, "Get current time")
	  .property("trans_external", &SEIRstocN::get_trans_external, &SEIRstocN::set_trans_external, "Get/set current trans_external")
	  .property("transmission_type", &SEIRstocN::get_trans_type, &SEIRstocN::set_trans_type, "Get/set current trans_type")
	; 

	using SEIRdet3 = SEIRmodel<update_type::deterministic, true, 3L, true>;
	class_<SEIRdet3>("SEIRdet3")
	  .constructor<const int, const Rcpp::CharacterVector>("Constructor")
	  .factory(invalidate_default_constructor)
	    
		.method("show", &SEIRdet3::show, "The show/print method")
	  .method("update", &SEIRdet3::update, "The update method")  
	  .method("run", &SEIRdet3::run, "The run method")  
	  .method("save", &SEIRdet3::save, "The save method")  
	  .method("reset", &SEIRdet3::reset, "The reset method")  
	    
		.property("state", &SEIRdet3::get_state, "Get current state")
	  .property("S", &SEIRdet3::get_S, &SEIRdet3::set_S, "Get/set current S")
	  .property("E", &SEIRdet3::get_E, &SEIRdet3::set_E, "Get/set current E")
	  .property("I", &SEIRdet3::get_I, &SEIRdet3::set_I, "Get/set current I")
	  .property("R", &SEIRdet3::get_R, &SEIRdet3::set_R, "Get/set current R")
	  .property("N", &SEIRdet3::get_N, "Get current N")
	  .property("beta", &SEIRdet3::get_beta, &SEIRdet3::set_beta, "Get/set current beta")
	  .property("omega", &SEIRdet3::get_omega, &SEIRdet3::set_omega, "Get/set current omega")
	  .property("gamma", &SEIRdet3::get_gamma, &SEIRdet3::set_gamma, "Get/set current gamma")
	  .property("delta", &SEIRdet3::get_delta, &SEIRdet3::set_delta, "Get/set current delta")
	  .property("repl", &SEIRdet3::get_repl, &SEIRdet3::set_repl, "Get/set current repl")
	  .property("cull", &SEIRdet3::get_cull, &SEIRdet3::set_cull, "Get/set current cull")
	  .property("vacc", &SEIRdet3::get_vacc, &SEIRdet3::set_vacc, "Get/set current vacc")
	  .property("time", &SEIRdet3::get_time, "Get current time")
	  .property("trans_external", &SEIRdet3::get_trans_external, &SEIRdet3::set_trans_external, "Get/set current trans_external")
	  .property("transmission_type", &SEIRdet3::get_trans_type, &SEIRdet3::set_trans_type, "Get/set current trans_type")
	; 

	using SEIRstoc3 = SEIRmodel<update_type::stochastic, true, 3L, true>;
	class_<SEIRstoc3>("SEIRstoc3")
	  .constructor<const int, const Rcpp::CharacterVector>("Constructor")
	  .factory(invalidate_default_constructor)
	    
		.method("show", &SEIRstoc3::show, "The show/print method")
	  .method("update", &SEIRstoc3::update, "The update method")  
	  .method("run", &SEIRstoc3::run, "The run method")  
	  .method("save", &SEIRstoc3::save, "The save method")  
	  .method("reset", &SEIRstoc3::reset, "The reset method")  
	    
		.property("state", &SEIRstoc3::get_state, "Get current state")
	  .property("S", &SEIRstoc3::get_S, &SEIRstoc3::set_S, "Get/set current S")
	  .property("E", &SEIRstoc3::get_E, &SEIRstoc3::set_E, "Get/set current E")
	  .property("I", &SEIRstoc3::get_I, &SEIRstoc3::set_I, "Get/set current I")
	  .property("R", &SEIRstoc3::get_R, &SEIRstoc3::set_R, "Get/set current R")
	  .property("N", &SEIRstoc3::get_N, "Get current N")
	  .property("beta", &SEIRstoc3::get_beta, &SEIRstoc3::set_beta, "Get/set current beta")
	  .property("omega", &SEIRstoc3::get_omega, &SEIRstoc3::set_omega, "Get/set current omega")
	  .property("gamma", &SEIRstoc3::get_gamma, &SEIRstoc3::set_gamma, "Get/set current gamma")
	  .property("delta", &SEIRstoc3::get_delta, &SEIRstoc3::set_delta, "Get/set current delta")
	  .property("repl", &SEIRstoc3::get_repl, &SEIRstoc3::set_repl, "Get/set current repl")
	  .property("cull", &SEIRstoc3::get_cull, &SEIRstoc3::set_cull, "Get/set current cull")
	  .property("vacc", &SEIRstoc3::get_vacc, &SEIRstoc3::set_vacc, "Get/set current vacc")
	  .property("time", &SEIRstoc3::get_time, "Get current time")
	  .property("trans_external", &SEIRstoc3::get_trans_external, &SEIRstoc3::set_trans_external, "Get/set current trans_external")
	  .property("transmission_type", &SEIRstoc3::get_trans_type, &SEIRstoc3::set_trans_type, "Get/set current trans_type")
	; 


}

/*

  DO NOT MODIFY THIS FILE DIRECTLY!
  See module/ipdmr_module_template.cpp

*/ 

