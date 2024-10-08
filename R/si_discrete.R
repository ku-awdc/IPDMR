#' Title
#'
#' @param N
#' @param beta
#' @param type
#' @param init_I
#' @param time_points
#'
#' @importFrom checkmate qassert assert_number
#'
#' @example si_discrete(N=10, type="density") |> autoplot()
#' @example si_discrete(N=10, type="frequency") |> autoplot()
#'
#' @export
si_discrete <- function(N=10, beta=0.05, type=c("frequency","density"), init_I=1, time_step=1/24, max_time=21){

  qassert(N, "N1(0,)")
  qassert(beta, "N1(0,)")
  qassert(init_I, "N1(0,)")
  assert_number(init_I, upper=N)
  qassert(N, "N1(0,)")
  qassert(time_step, "N1[0,)")
  qassert(max_time, "N1[0,)")
  assert_number(max_time, lower=time_step)

  ntime <- ceiling(max_time/time_step)+1L
  beta <- beta*time_step
  type <- match.arg(type)

  i_dens <- function(i) 1 - exp(-beta * i)
  i_freq <- function(i) 1 - exp(-beta * i/N)
  i_fun <- if(type=="density") i_dens else i_freq

  output <- matrix(NA_real_, nrow=ntime+1L, ncol=3L, dimnames=list(NULL, c("Time","S","I")))
  I <- init_I
  S <- N - I
  time <- 0
  output[1L,"Time"] <- time
  output[1L,"S"] <- S
  output[1L,"I"] <- I
  for(row in seq_len(ntime)){
    time <- time + time_step
    newI <- S*i_fun(I)
    S <- S - newI
    I <- I + newI
    output[row+1L,"Time"] <- time
    output[row+1L,"S"] <- S
    output[row+1L,"I"] <- I
  }

  output |>
    as_tibble() |>
    select(Time, everything()) |>
    filter(Time <= max_time) ->
    output

  class(output) <- c("ipdmr_dt", class(output))
  return(output)
}


