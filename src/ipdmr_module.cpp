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

	using SEIRdb = SEIRmodel<true>;
	class_<SEIRdb>("SEIRdb")
	  .constructor<const int, const double, const Rcpp::CharacterVector>("Constructor")
	  .factory(invalidate_default_constructor)
		.method("show", &SEIRdb::show, "The show/print method")
		.property("state", &SEIRdb::get_state, "Get current state")
	; 

	using SEIR = SEIRmodel<false>;
	class_<SEIR>("SEIR")
	  .constructor<const int, const double, const Rcpp::CharacterVector>("Constructor")
	  .factory(invalidate_default_constructor)
		.method("show", &SEIR::show, "The show/print method")
		.property("state", &SEIR::get_state, "Get current state")
	; 


}

/*

  DO NOT MODIFY THIS FILE DIRECTLY!
  See module/ipdmr_module_template.cpp

*/ 

