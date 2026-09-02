plot_targets <- list(
  # Figures.
  tar_target(
    plot_dag,
    create_plot_dag()
  ),
  tar_target(
    foodweb_fig,
    {
      dot_file <- "figures/foodweb.dot"
      out_file <- "figures/foodweb.png"
      system2("dot", c("-Tpng", dot_file, "-o", out_file))
      out_file
    },
    format = "file"
  ),
  tar_target(
    venn_diagram,
    {
      fname <- "figures/venn-diagram.png"
      ggsave(fname, plot_venn_diagram(sea_data_tidy$size))
    },
    format = "file"
  ),
  tar_target(
    plot_station,
    plot_sampling(sea_data_tidy, station_file)
  ),
  tar_target(
    plot_stations_overview,
    plot_station_overview(
      operation_location, params$district_kept, basin, lang = "fr"
    )
  ),
  tar_target(
    plot_net_groups,
    plot_network_groups(diet_resource)
  ),
  tar_target(
    plot_diet_categories_overview,
    plot_diet_categories(diet, size_year_filtered)
  ),
  tar_target(
    plot_foodweb_overview_report,
    plot_foodweb_overview(
      web_list$metaweb, resource_list, diet, size_year_filtered,
      foodweb_structure,
      ops = c(river = "83544", nurse = "2022_1_GIR12_14_S", pomet = "61832307"),
      titles = c(
        river = "Rivière (Dordogne)",
        nurse = "Estuaire (Gironde)",
        pomet = "Côte (Gironde)"
      ),
      lang = "fr"
    )
  ),
  tar_target(
    plot_foodweb_structure_map_overview,
    plot_foodweb_structure_map(
      foodweb_structure, operation_location, params$district_kept, basin
    )
  ),
  tar_target(
    plot_rephy_predictions_map_overview,
    plot_rephy_predictions_map(
      rephy_operation_predictions, operation_location, params$district_kept, basin
    )
  ),
  tar_target(
    plot_distance_to_mouth_overview_report,
    plot_distance_to_mouth_overview(
      river_network_geo, river_mouth_district, inland_distance_to_mouth,
      sea_graph, sea_mesh, sea_offshore_distance_to_mouth,
      distance_to_mouth, params$district_kept, basin,
      river_operation_id = "61833954",
      sea_operation_id = "2019_2_GIR33_1_S",
      lang = "fr"
    )
  )
)
