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

