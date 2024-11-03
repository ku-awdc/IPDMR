library("tidyverse")

cat("/*

  DO NOT MODIFY THIS FILE DIRECTLY!
  See module/ipdmr_module_template.cpp

*/ \n\n", file="src/ipdmr_module.cpp", append=FALSE)

cat(readLines("src/module/ipdmr_module_template.cpp"), sep="\n", file="src/ipdmr_module.cpp", append=TRUE)

ct <- paste0("\t", readLines("src/module/class_template.cpp")) |> paste(collapse="\n")

tribble(~Name, ~Template,
  "SEIRdb", "SEIRmodel<true>",
  "SEIR", "SEIRmodel<false>",
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
