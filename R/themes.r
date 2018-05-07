#' Blackbox theme
#' 
#' ggplot2 theme that places a black margin on every
#' panel
#'
#' @param border.size Thickness of the margin.
#'
#' @return A \code{theme gg} object.
#' @export
#' @author Sur Herrera Paredes
theme_blackbox <- function(border.size = 3){
  theme(axis.text = element_text(color="black"),
        axis.title = element_text(face="bold"),
        panel.background = element_rect(color = "black",
                                        size = border.size,
                                        fill = NA),
        panel.grid = element_blank())
}
  
  