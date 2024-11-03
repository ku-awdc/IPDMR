class_<NAME>("NAME")
  .constructor<const int, const double, const Rcpp::CharacterVector>("Constructor")
  .factory(invalidate_default_constructor)
	.method("show", &NAME::show, "The show/print method")
	.property("state", &NAME::get_state, "Get current state")
;
