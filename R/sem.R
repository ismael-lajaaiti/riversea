#' Build a SEM.
#'
#' @param observable list of observable variables
#' @param mediator list of mediator (intemediate) variables
#' @param response list of response variable that we want to predict
#' @param random random effects
#' @param data data frame
#'
#' @return psem output
#' @export
build_sem <- function(observable, mediator, response, random, data) {
  # Observables to Mediators.
  mediator_model <- purrr::map(mediator, function(med) {
    glmmTMB::glmmTMB(make_formula(med, observable, random), data = data)
  })
  # Mediators to Responses.
  response_model <- purrr::map(response, function(resp) {
    glmmTMB::glmmTMB(make_formula(resp, mediator, random), data = data)
  })
  do.call(
    piecewiseSEM::psem,
    c(mediator_model, response_model, list(data = data))
  )
}

make_formula <- function(response, predictor, random) {
  rhs <- paste(c(predictor, random), collapse = "+")
  stats::as.formula(paste(response, "~", rhs))
}

build_sem_graph <- function(model) {
  edges <- build_sem_edges(model)
  nodes <- tibble(
    name = unique(c(edges$from, edges$to)),
    raw_label = name |>
      stringr::str_replace_all("_", " ") |>
      stringr::str_to_title(),
    label = case_when(
      raw_label == "Log Species Richness" ~ "Species\nRichness",
      raw_label == "Frac Piscivorous" ~ "%\nPiscivorous",
      raw_label == "Max Fish Size" ~ "Max\nFish Size",
      TRUE ~ raw_label
    )
  )
  # Plot SEM graph.
  tidygraph::tbl_graph(nodes = nodes, edges = edges, directed = TRUE)
}

build_sem_edges <- function(model) {
  coef_table <- piecewiseSEM::coefs(model, standardize = "scale")
  column_name <- names(coef_table)
  column_name[length(column_name)] <- "Significant"
  names(coef_table) <- column_name
  edges <- coef_table |>
    select(
      from = Predictor,
      to = Response,
      estimate = Std.Estimate,
      p = P.Value
    ) |>
    mutate(
      significant = p < 0.05,
      sign = if_else(estimate > 0, "positive", "negative")
    )
}

#' Compute net effects of observables on responses.
#'
#' @param model piecewiseSEM model
#'
#' @return data frame
#' @export
get_sem_net_effects <- function(model) {
  edges <- build_sem_edges(model) |>
    mutate(estimate = if_else(significant, estimate, 0))
  direct <- edges |>
    filter(
      from %in% observable,
      to %in% response
    ) |>
    select(from, to, direct = estimate)
  indirect <- edges |>
    filter(from %in% observable) |>
    rename(mediator = to, path1 = estimate) |>
    inner_join(
      edges |>
        filter(to %in% response) |>
        rename(mediator = from, path2 = estimate),
      by = "mediator"
    ) |>
    mutate(indirect = path1 * path2) |>
    group_by(from, to) |>
    summarise(indirect = sum(indirect), .groups = "drop")
  indirect |>
    left_join(direct, by = c("from", "to")) |>
    mutate(
      direct = tidyr::replace_na(direct, 0),
      total = indirect + direct,
      sign = if_else(total > 0, "positive", "negative"),
      from_label = from |>
        stringr::str_replace_all("_", " ") |>
        stringr::str_to_title(),
      to_label = to |>
        stringr::str_replace_all("_", " ") |>
        stringr::str_to_title()
    )
}

#' Plot SEM graph
#'
#' @param model piecewiseSEM model
#'
#' @return plot
#' @export
plot_sem <- function(model) {
  graph <- build_sem_graph(model)
  p_graph <- ggraph::ggraph(graph, layout = "sugiyama") +
    ggraph::geom_edge_link(
      ggplot2::aes(
        edge_color = sign,
        edge_linetype = if_else(significant, "solid", "dashed"),
        label = if_else(significant, as.character(round(estimate, 2)), "")
      ),
      arrow = grid::arrow(length = grid::unit(4, "mm"), type = "closed"),
      end_cap = ggraph::rectangle(23, 11, "mm"),
      angle_calc = "along",
      label_dodge = grid::unit(2.5, "mm"),
      label_push = grid::unit(20, "mm")
    ) +
    ggraph::geom_node_label(ggplot2::aes(label = label)) +
    ggraph::scale_edge_color_manual(values = c(
      "positive" = "springgreen3",
      "negative" = "firebrick3"
    )) +
    ggraph::theme_graph() +
    ggplot2::theme(legend.position = "none")
  net_effects <- get_sem_net_effects(model)
  p_effects <- ggplot2::ggplot(
    net_effects,
    ggplot2::aes(x = total, y = from_label, fill = sign)
  ) +
    ggplot2::geom_col() +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    ggplot2::facet_wrap(~to_label, ncol = 1, scales = "free_x") +
    ggplot2::scale_fill_manual(values = c(
      "positive" = "springgreen3",
      "negative" = "firebrick3"
    )) +
    ggplot2::labs(x = "Net effect", y = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")
  patchwork::wrap_plots(p_graph, p_effects)
}
