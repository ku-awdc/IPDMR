library("tidyverse")

cat("/*

  DO NOT MODIFY THIS FILE DIRECTLY!
  See module/ipdmr_module_template.cpp

*/ \n\n", file="src/ipdmr_module.cpp", append=FALSE)

cat(readLines("src/module/ipdmr_module_template.cpp"), sep="\n", file="src/ipdmr_module.cpp", append=TRUE)

ct <- paste0("\t", readLines("src/module/class_template.cpp")) |> paste(collapse="\n")

expand_grid(ut=c("deterministic","stochastic"), ne=c(-1,0,1,3), ni=c(-1,1), nr=c(-1,0,1)) |>
  mutate(Debug = "true") |>
  mutate(Name = str_c("SEIR",
    if_else(ut=="deterministic", "det", "stoc"),
    if_else(ne==-1, "N", as.character(ne)),
    if_else(ni==-1, "N", as.character(ni)),
    if_else(nr==-1, "N", as.character(nr))
  )) |>
  mutate(Template = str_c("SEIRmodel<update_type::", ut, ",", ne, ",", ni, ",", nr, ",", Debug, ">")) |>
  rowwise() |>
  group_split() |>
  lapply(function(x){
    n <- x$Name
    t <- x$Template
    cat("\tusing ", n, " = ", t, ";\n", sep="", file="src/ipdmr_module.cpp", append=TRUE)
    gsub("NAME", n, ct) |> cat("\n\n", file="src/ipdmr_module.cpp", append=TRUE)
  })

cat("\n}\n\n/*

  DO NOT MODIFY THIS FILE DIRECTLY!
  See module/ipdmr_module_template.cpp

*/ \n\n", file="src/ipdmr_module.cpp", append=TRUE)
