library("ggplot2")

#font
t <- list(
  family = "Helvetica",
  size = 18,
  face = "bold",
  color = 'black')

fontfamily <- "Helvetica"

#xaxis plotly layout
xaxis <- list(
  showgrid = FALSE,
  showline = TRUE,
  linewidth = 2,
  mirror = TRUE,
  ticks = "outside",
  zeroline = FALSE,
  font = t
)

#yaxis plotly layout
yaxis <- list(
  mirror = TRUE,
  linewidth = 2,
  tickwidth = 2,
  tickformat = "digits",
  showgrid = FALSE,
  zeroline = FALSE,
  showline = TRUE,
  font = t,
  ticks = "outside",
  hoverformat = '.2f')

#ggplot theme
ggplottheme <- theme(axis.text.x=element_text(colour="black", face = "bold", family = fontfamily, size = 4.76 * .pt),
                    axis.text.y=element_text(colour="black", face = "bold", family = fontfamily),
                    axis.title.x = element_text(family = fontfamily),
                    axis.title.y = element_text(family = fontfamily),
                    axis.ticks = element_line(size = 2, color="black"),
                    axis.ticks.length=unit(.2, "cm"),
                    legend.title=element_text(face = "bold", family = fontfamily, size=5.81 * .pt),
                    legend.text=element_text(face = "bold", family = fontfamily, size=4.76 * .pt),
                    legend.margin = margin(0, 0, 0, 0),
                    legend.spacing.y = unit(0, "mm")
                    )
