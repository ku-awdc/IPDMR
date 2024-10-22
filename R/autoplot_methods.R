#' Title
#'
#' @param object
#' @param ...
#'
#' @import ggplot2
#' @importFrom tidyr pivot_longer
#' @importFrom rlang .data

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

#' @export
autoplot.ipdmr_st <- function(object, ...){

  if(attr(object,'iterations') == 1L){

    object |>
      select(-"Iteration") |>
      pivot_longer(cols=-"Time", names_to="Compartment", values_to="Number") |>
      ggplot(aes(x=.data$Time, y=.data$Number, col=.data$Compartment)) +
      geom_step() +
      theme_light() +
      theme(legend.position="bottom", legend.title = element_blank()) +
      ylim(c(0,object |> slice(1L) |> select(-"Time") |> sum())) +
      labs(caption = attr(object, "plot_caption"))

  }else if(attr(object,'iterations') < 5L){

    object |>
      pivot_longer(cols=!"Time" & !"Iteration", names_to="Compartment", values_to="Number") |>
      ggplot(aes(x=.data$Time, y=.data$Number, col=.data$Compartment)) +
      geom_step() +
      theme_light() +
      theme(legend.position="bottom", legend.title = element_blank()) +
      ylim(c(0,object |> slice(1L) |> select(-"Time") |> sum())) +
      labs(caption = attr(object, "plot_caption")) +
      facet_wrap(~str_c("Iteration: ", Iteration))

  }else{

    object |>
      pivot_longer(cols=!"Time" & !"Iteration", names_to="Compartment", values_to="Number") |>
      group_by(Time, Compartment) |>
      summarise(Mean = mean(Number), LCI = quantile(Number, 0.025), UCI = quantile(Number, 0.975), .groups="drop") |>
      ggplot(aes(x=.data$Time, y=.data$Mean, ymin=.data$LCI, ymax=.data$UCI, col=.data$Compartment, fill=.data$Compartment)) +
      geom_ribbon(alpha=0.25, col="transparent") +
      geom_step() +
      theme_light() +
      theme(legend.position="bottom", legend.title = element_blank()) +
      ylim(c(0,object |> slice(1L) |> select(-"Time") |> sum())) +
      labs(caption = attr(object, "plot_caption"))

  }

}
