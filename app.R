# Optional single-file launcher. You can keep ui.R/server.R for shinyapps.io;
# app.R is useful for local development with shiny::runApp().
source("R/packages.R")
load_rnaseviz_packages()
source("R/data_config.R")
source("R/helpers_data.R")
source("R/helpers_positions.R")
source("R/plot_alignment.R")
source("R/ui_components.R")
source("R/server.R")

shiny::shinyApp(ui = rnaseviz_ui(), server = rnaseviz_server)
