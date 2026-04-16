plot_targets <- list(
  # Figures.
  tar_target(
    plot_dag,
    create_plot_dag()
  ),
  tar_target(
    foodweb_fig,
    {
      dot_file <- here("figures", "foodweb.dot")
      out_file <- here("figures", "foodweb.png")
      system2("dot", c("-Tpng", dot_file, "-o", out_file))
      out_file
    },
    format = "file"
  ),
  tar_target(
    venn_diagram,
    {
      fname <- here("figures", "venn-diagram.png")
      ggsave(fname, plot_venn_diagram(sea_data_tidy$size))
    },
    format = "file"
  ),
  tar_target(
    plot_station,
    plot_sampling(sea_data_tidy, station_file)
  ),
  tar_target(
    plot_metaweb_connectance,
    plot_sizeclass_connectance(metaweb_table)
  )
)
