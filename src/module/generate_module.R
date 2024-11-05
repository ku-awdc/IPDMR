library("tidyverse")

cat("/*

  DO NOT MODIFY THIS FILE DIRECTLY!
  See module/ipdmr_module_template.cpp

*/ \n\n", file="src/ipdmr_module.cpp", append=FALSE)

cat(readLines("src/module/ipdmr_module_template.cpp"), sep="\n", file="src/ipdmr_module.cpp", append=TRUE)

ct <- paste0("\t", readLines("src/module/class_template.cpp")) |> paste(collapse="\n")

tribble(~Name, ~Template,
  "SEIRdetN", "SEIRmodel<update_type::deterministic, false, 0L, true>",
  "SEIRstocN", "SEIRmodel<update_type::stochastic, false, 0L, true>",
  "SEIRdet3", "SEIRmodel<update_type::deterministic, true, 3L, true>",
  "SEIRstoc3", "SEIRmodel<update_type::stochastic, true, 3L, true>",
) |>
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
