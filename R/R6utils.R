`$` <- function(lhs,rhs){
  name <- deparse(substitute(rhs))
  if(is.environment(lhs) &&
      !".pointer" %in% names(lhs) &&
      !name %in% c(names(lhs),"active","private","super","deep_clone") &&
      !str_sub(name, 1L, 1L)==".")
    stop(name, " not present in environment")
  do.call(base::`$`, args=list(lhs, name))
}

mlist <- function(..., .dummy=NULL){
  args <- list(...)
  if(any(names(args)=="")){
    an <- names(args)
    an[names(args)==""] <- "MISSING"
    stop("Empty names in member/method list: ", paste0(an,collapse=", "))
  }
  if(any(table(names(args))>1)){
    tn <- table(names(args))
    stop("Duplicated names in member/method list: ", paste0(names(tn)[tn>1], collapse=", "))
  }
  do.call(list, args)
}
