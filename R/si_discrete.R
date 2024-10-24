#' Title
#'
#' @param S
#' @param I
#' @param beta
#' @param type
#' @param d_time
#' @param max_time
#'
#' @importFrom checkmate qassert assert_number
#'
#' @examples
#' si_discrete(N=10, type="density") |> ggplot2::autoplot()
#' si_discrete(N=10, type="frequency") |> ggplot2::autoplot()
#'
#' @export
si_discrete <- function(S=9, I=1, beta=0.05, type=c("frequency","density"), d_time=1/24, max_time=21){

  qassert(S, "N1[0,)")
  qassert(I, "N1[0,)")
  qassert(S+I, "N1(0,)")
  N <- S+I
  qassert(beta, "N1(0,)")
  qassert(d_time, "N1[0,)")
  qassert(max_time, "N1[0,)")
  assert_number(max_time, lower=d_time)

  ntime <- ceiling(max_time/d_time)+1L
  type <- match.arg(type)

  i_dens <- function(i) 1 - exp(-beta * i)
  i_freq <- function(i) 1 - exp(-beta * i/N)
  i_fun <- if(type=="density") i_dens else i_freq

  output <- matrix(NA_real_, nrow=ntime+1L, ncol=3L, dimnames=list(NULL, c("Time","S","I")))
  time <- 0
  output[1L,"Time"] <- time
  output[1L,"S"] <- S
  output[1L,"I"] <- I
  for(row in seq_len(ntime)){
    time <- time + d_time
    newI <- S*i_fun(I) * d_time
    S <- S - newI
    I <- I + newI
    output[row+1L,"Time"] <- time
    output[row+1L,"S"] <- S
    output[row+1L,"I"] <- I
  }

  output |>
    as_tibble() |>
    select(.data$Time, everything()) |>
    filter(.data$Time <= max_time) ->
    output

  class(output) <- c("ipdmr_dt", class(output))
  attr(output, "plot_caption") <- str_c("discrete; ", type, "; d_time=", round(d_time,3))
  return(output)
}


