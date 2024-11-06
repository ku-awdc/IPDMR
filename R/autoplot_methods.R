#' Autoplot methods associated with IPDMR model functions
#'
#' @param object an object returned by an IPDMR model function
#' @param ... additional arguments (currently ignored)
#'
#' @import ggplot2
#' @importFrom tidyr pivot_longer replace_na
#' @importFrom rlang .data set_names

#' @export
autoplot.ipdmr_ct <- function(object, ...){

  object |>
    pivot_longer(cols=-"Time", names_to="Compartment", values_to="Number") |>
    mutate(Compartment = fctcomp(Compartment)) |>
    ggplot(aes(x=.data$Time, y=.data$Number, col=.data$Compartment)) +
    geom_line(alpha=0.25) +
    geom_point(size=0.5) +
    theme_light() +
    scale_color_manual(values=compcol) + scale_fill_manual(values=compcol) +
    theme(legend.position="bottom", legend.title = element_blank()) +
    ylim(c(0,object |> slice(1L) |> select(-"Time") |> sum())) +
    labs(caption = attr(object, "plot_caption"))

}

#' @export
autoplot.ipdmr_dt <- function(object, ...){

  object |>
    pivot_longer(cols=-"Time", names_to="Compartment", values_to="Number") |>
    mutate(Compartment = fctcomp(Compartment)) |>
    ggplot(aes(x=.data$Time, y=.data$Number, col=.data$Compartment)) +
    geom_step(alpha=0.25) +
    geom_point(size=0.5) +
    theme_light() +
    scale_color_manual(values=compcol) + scale_fill_manual(values=compcol) +
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
      mutate(Compartment = fctcomp(Compartment)) |>
      ggplot(aes(x=.data$Time, y=.data$Number, col=.data$Compartment)) +
      geom_step() +
      theme_light() +
      scale_color_manual(values=compcol) + scale_fill_manual(values=compcol) +
      theme(legend.position="bottom", legend.title = element_blank()) +
      ylim(c(0,object |> slice(1L) |> select(-"Time") |> sum())) +
      labs(caption = attr(object, "plot_caption"))

  }else if(attr(object,'iterations') < 5L){

    object |>
      pivot_longer(cols=!"Time" & !"Iteration", names_to="Compartment", values_to="Number") |>
      mutate(Compartment = fctcomp(Compartment)) |>
      ggplot(aes(x=.data$Time, y=.data$Number, col=.data$Compartment)) +
      geom_step() +
      theme_light() +
      scale_color_manual(values=compcol) + scale_fill_manual(values=compcol) +
      theme(legend.position="bottom", legend.title = element_blank()) +
      ylim(c(0,object |> slice(1L) |> select(-"Time") |> sum())) +
      labs(caption = attr(object, "plot_caption")) +
      facet_wrap(~str_c("Iteration: ", Iteration))

  }else{

    object |>
      pivot_longer(cols=!"Time" & !"Iteration", names_to="Compartment", values_to="Number") |>
      mutate(Compartment = fctcomp(Compartment)) |>
      group_by(Time, Compartment) |>
      summarise(Mean = mean(Number), LCI = quantile(Number, 0.025), UCI = quantile(Number, 0.975), .groups="drop") |>
      ggplot(aes(x=.data$Time, y=.data$Mean, ymin=.data$LCI, ymax=.data$UCI, col=.data$Compartment, fill=.data$Compartment)) +
      geom_ribbon(alpha=0.25, col="transparent") +
      geom_step() +
      theme_light() +
      scale_color_manual(values=compcol) + scale_fill_manual(values=compcol) +
      theme(legend.position="bottom", legend.title = element_blank()) +
      ylim(c(0,object |> slice(1L) |> select(-"Time") |> sum())) +
      labs(caption = attr(object, "plot_caption"))

  }

}


#' @export
autoplot.ipdmr_dm <- function(object, ...){

  object$GroupName <- NULL
  object$GroupIndex <- NULL

  if(attr(object,'ngroups') <= 5L){

    object |>
      pivot_longer(cols=c(-"Time",-"Group"), names_to="Compartment", values_to="Number") |>
      mutate(Compartment = fctcomp(Compartment)) |>
      ggplot(aes(x=.data$Time, y=.data$Number, col=.data$Compartment)) +
      geom_step(alpha=0.25) +
      geom_point(size=0.5) +
      facet_wrap(~Group, scales="free_y", ncol=1) +
      theme_light() +
      scale_color_manual(values=compcol) + scale_fill_manual(values=compcol) +
      theme(legend.position="bottom", legend.title = element_blank()) +
      labs(caption = attr(object, "plot_caption"))

  }else{

    object |>
      pivot_longer(cols=c(-"Time",-"Group"), names_to="Compartment", values_to="Number") |>
      mutate(Compartment = fctcomp(Compartment)) |>
      group_by(Time, Compartment) |>
      summarise(TotalAnimals = sum(Number), InfectedGroups = sum(Number>0), .groups="drop") |>
      pivot_longer(TotalAnimals:InfectedGroups, names_to="Metric", values_to="Number") |>
      filter(Metric=="TotalAnimals" | Compartment=="I") |>
      ggplot(aes(x=.data$Time, y=.data$Number, col=.data$Compartment, fill=.data$Compartment)) +
      geom_step(alpha=0.25) +
      geom_point(size=0.5) +
      facet_wrap(~Metric, scales="free_y", ncol=2) +
      theme_light() +
      scale_color_manual(values=compcol) + scale_fill_manual(values=compcol) +
      theme(legend.position="bottom", legend.title = element_blank()) +
      labs(caption = attr(object, "plot_caption"))

  }

}


#' @export
autoplot.ipdmr_sm <- function(object, ...){

  object$GroupName <- NULL
  object$GroupIndex <- NULL

  if(attr(object,'iterations') == 1L){

    if(attr(object,'ngroups') <= 5L){

      object |>
        pivot_longer(cols=c(-"Time",-"Group"), names_to="Compartment", values_to="Number") |>
        mutate(Compartment = fctcomp(Compartment)) |>
        ggplot(aes(x=.data$Time, y=.data$Number, col=.data$Compartment)) +
        geom_step(alpha=0.25) +
        geom_point(size=0.5) +
        facet_wrap(~Group, scales="free_y", ncol=1) +
        theme_light() +
        scale_color_manual(values=compcol) + scale_fill_manual(values=compcol) +
        theme(legend.position="bottom", legend.title = element_blank()) +
        labs(caption = attr(object, "plot_caption"))

    }else{

      object |>
        pivot_longer(cols=c(-"Time",-"Group"), names_to="Compartment", values_to="Number") |>
        mutate(Compartment = fctcomp(Compartment)) |>
        group_by(Time, Compartment) |>
        summarise(TotalAnimals = sum(Number), InfectedGroups = sum(Number>0), .groups="drop") |>
        pivot_longer(TotalAnimals:InfectedGroups, names_to="Metric", values_to="Number") |>
        filter(Metric=="TotalAnimals" | Compartment=="I") |>
        ggplot(aes(x=.data$Time, y=.data$Number, col=.data$Compartment, fill=.data$Compartment)) +
        geom_step(alpha=0.25) +
        geom_point(size=0.5) +
        facet_wrap(~Metric, scales="free_y", ncol=2) +
        theme_light() +
        scale_color_manual(values=compcol) + scale_fill_manual(values=compcol) +
        theme(legend.position="bottom", legend.title = element_blank()) +
        labs(caption = attr(object, "plot_caption"))

    }

  }else{

    if(attr(object,'ngroups') <= 5L){

      object |>
        pivot_longer(cols=c(-"Time",-"Group",-"Iteration"), names_to="Compartment", values_to="Number") |>
        mutate(Compartment = fctcomp(Compartment)) |>
        group_by(Time, Compartment, Group) |>
        summarise(Mean = mean(Number), LCI = quantile(Number, 0.025), UCI = quantile(Number, 0.975), .groups="drop") |>
        ggplot(aes(x=.data$Time, y=.data$Mean, ymin=.data$LCI, ymax=.data$UCI, col=.data$Compartment, fill=.data$Compartment)) +
        geom_ribbon(alpha=0.25, col="transparent") +
        geom_step() +
        facet_wrap(~Group, scales="free_y", ncol=1) +
        theme_light() +
        scale_color_manual(values=compcol) + scale_fill_manual(values=compcol) +
        theme(legend.position="bottom", legend.title = element_blank()) +
        labs(caption = attr(object, "plot_caption"))

    }else{

      object |>
        pivot_longer(cols=c(-"Time",-"Group"), names_to="Compartment", values_to="Number") |>
        group_by(Time, Compartment) |>
        summarise(TotalAnimals = sum(Number), InfectedGroups = sum(Number>0), .groups="drop") |>
        pivot_longer(TotalAnimals:InfectedGroups, names_to="Metric", values_to="Number") |>
        filter(Metric=="TotalAnimals" | Compartment=="I") |>
        ggplot(aes(x=.data$Time, y=.data$Number, col=.data$Compartment, fill=.data$Compartment)) +
        geom_step(alpha=0.25) +
        geom_point(size=0.5) +
        facet_wrap(~Metric, scales="free_y", ncol=2) +
        theme_light() +
        theme(legend.position="bottom", legend.title = element_blank()) +
        scale_color_manual(values=compcol) + scale_fill_manual(values=compcol) +
        labs(caption = attr(object, "plot_caption"))

    }

  }
}


fctcomp <- function(Compartment){
  poss <- c("S","E","I","R")
  factor(Compartment, levels=poss[poss %in% Compartment])
}
compcol <- c("S"="#7CAE00", "E"="#C77CFF", "I"="#F8766D", "R"="#00BFC4")
