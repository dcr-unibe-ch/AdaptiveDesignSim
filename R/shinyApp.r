#' AdaptiveDesignSim shiny app
#'
#' The shiny app to AdaptiveDesignSim provides an intuitive way to 
#' 1) check the results of individual simulations (via forest plots) and 
#' 2) of operating characteristics over a range of conditions.
#'
#' @importFrom shiny runApp
#'
#' @usage
#' launch_ADSim_app()
#' @export
#' @examples
#' # launch the app
#' \dontrun{
#' launch_ADSim_app()
#' }
#'

launch_ADSim_app <- function(){
  shiny::runApp(system.file("shinyApp", package = "AdaptiveDesignSim"))
}
