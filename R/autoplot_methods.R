#' Title
#'
#' @param object
#' @param ...
#'
#' @import ggplot2
#' @importFrom tidyr pivot_longer
#' @importFrom rlang .data
#'
#' @export
autoplot.ipdmr_dt <- function(object, ...){

  object |>
    pivot_longer(cols=-"Time", names_to="Compartment", values_to="Number") |>
    ggplot(aes(x=.data$Time, y=.data$Number, col=.data$Compartment)) +
    geom_line() +
    geom_point() +
    theme_light() +
    theme(legend.position="bottom", legend.title = element_blank()) +
    ylim(c(0,object |> slice(1L) |> select(-"Time") |> sum())) +
    labs(caption = attr(object, "plot_caption"))

}
