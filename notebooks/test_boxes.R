seq(1,10,by=1) |>
  lapply(\(x){
    model <- make_group(numE=0, numI=1, numR=x)
    model$S <- 0
    model$I <- 0
    model$R <- 1
    model$beta <- 0
    model$gamma <- 0.1
    model$delta <- 0.1
    model$run(100,0.1) |> mutate(B=x, New=S-lag(S))
  }) |>
  bind_rows() |>
  ggplot(aes(x=Time, y=New, col=factor(B))) +
  geom_line() +
  geom_vline(xintercept=1/0.1)


## TODO: allow more flexibility in transmission type, e.g. function of N?


set.seed(2024-11-05)
r6res <- make_group(update_type = "deterministic", implementation = "R6")$run(10, 0.1)

set.seed(2024-11-05)
cppres <- make_group(update_type = "deterministic", implementation = "C++")$run(10, 0.1)

## Should be identical:
stopifnot(unlist(r6res) - unlist(cppres) == 0L)

