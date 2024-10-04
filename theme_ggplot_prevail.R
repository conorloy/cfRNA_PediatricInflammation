theme_prevail <- function(){
  theme_bw() %+replace%
    theme(#panel.grid.major = element_line(color='black'),
          #panel.grid.minor = element_line(color='black'),
          panel.border = element_rect(color='black', fill=NA),
          axis.title=element_text(family="Helvetica", size=8, color='black'),
          axis.ticks = element_line(color='black'),
          axis.text=element_text(family="Helvetica", size=6, color='black'),
          legend.position = "none")
}
