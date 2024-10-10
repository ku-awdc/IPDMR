#' Title
#'
#' @param object
#' @param ...
#'
#' @importFrom ggplot2 autoplot
#'
#' @export
autoplot.ipdmr_dt <- function(object, ...){

  object |>
    pivot_longer(cols=-Time, names_to="Compartment", values_to="Number") |>
    ggplot(aes(x=Time, y=Number, col=Compartment)) +
    geom_line() +
    geom_point() +
    theme_light() +
    theme(legend.pos="bottom", legend.title = element_blank()) +
    ylim(c(0,object |> slice(1L) |> select(-Time) |> sum()))

}
