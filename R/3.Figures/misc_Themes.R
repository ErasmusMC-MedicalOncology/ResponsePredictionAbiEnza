# Helper functions --------------------------------------------------------

themeAnno_Job <- ggplot2::theme(
  legend.position = 'right',
  axis.ticks = ggplot2::element_blank(),
  axis.title.y = ggtext::element_textbox_simple(size = 8, family='Nimbus Sans', face = 'bold', orientation = 'left-rotated', width = NULL, halign = .5),
  axis.text.x = ggplot2::element_blank(),
  text = ggplot2::element_text(size = 8, family='Nimbus Sans', face = 'bold'),
  panel.background = ggplot2::element_rect(fill = NA),
  panel.grid.major = ggplot2::element_line(NULL),
  panel.grid.major.x = ggplot2::element_line(NULL),
  panel.grid.major.y = ggplot2::element_line(NULL),
  plot.margin = ggplot2::unit(c(0, 0, 0, 0), "cm"),
  legend.text = ggtext::element_markdown(family='Nimbus Sans', face = 'bold'),
  legend.title = ggtext::element_markdown(family='Nimbus Sans', face = 'bold'),
  legend.key = ggplot2::element_blank()
)

themeTrack_Job <- ggplot2::theme(
  legend.position = 'right',
  legend.direction = 'horizontal',
  text = ggplot2::element_text(size = 7, family = 'Nimbus Sans', face = 'bold'),
  axis.title.y = ggtext::element_textbox_simple(size = 6, family='Nimbus Sans', face = 'bold', orientation = 'left-rotated', width = NULL, halign = .5),
  axis.text.x = ggplot2::element_blank(),
  axis.ticks.x = ggplot2::element_blank(),
  panel.grid.major.x = ggplot2::element_blank(),
  panel.grid.major.y = ggplot2::element_line(colour = '#E5E5E5', linetype = 'dotted'),
  panel.grid.minor.y = ggplot2::element_blank(),
  axis.title.x = ggplot2::element_blank(),
  panel.background = ggplot2::element_rect(fill = NA, colour = 'black'),
  panel.border = ggplot2::element_rect(fill = NA, colour = NA),
  strip.background = ggplot2::element_rect(colour = 'grey20', fill = 'white'),
  plot.margin = ggplot2::unit(c(0, 0, 0, 0), "cm"),
  legend.text = ggtext::element_markdown(family='Nimbus Sans', face = 'bold'),
  legend.title = ggtext::element_markdown(family='Nimbus Sans', face = 'bold'),
  legend.key = ggplot2::element_blank()
)

theme_Job <- ggplot2::theme(
  legend.position = 'bottom',
  legend.direction = 'horizontal',
  text = ggplot2::element_text(size=9, family='Nimbus Sans', face = 'bold'),
  axis.title.x = ggtext::element_textbox_simple(family='Nimbus Sans', face = 'bold', width = NULL, halign = .5),
  axis.title.y = ggtext::element_textbox_simple(family='Nimbus Sans', face = 'bold', size = 8, orientation = 'left-rotated', width = NULL, halign = .5),
  strip.text = ggtext::element_textbox_simple(width = NULL, halign = .5),
  panel.grid.major.x = ggplot2::element_line(colour = 'grey90', linetype = 'dotted'),
  panel.grid.major.y = ggplot2::element_line(colour = '#E5E5E5', linetype = 'dotted'),
  panel.grid.minor.y = ggplot2::element_blank(),
  panel.background = ggplot2::element_rect(fill = NA, colour = 'black'),
  panel.border = ggplot2::element_rect(fill = NA, colour = NA),
  strip.background = ggplot2::element_rect(colour = 'black', fill = 'white'),
  legend.text = ggtext::element_markdown(family='Nimbus Sans', face = 'bold'),
  legend.title = ggtext::element_markdown(family='Nimbus Sans', face = 'bold'),
  legend.key = ggplot2::element_blank()
)

# Half-Half Boxplot(s)
geomJob_HalfHalfBox <- list(
  gghalves::geom_half_boxplot(side = 'l', outlier.shape = NA, color = 'black'),
  gghalves::geom_half_point_panel(side = 'r', color = 'black', shape = 21, transformation = ggbeeswarm::position_quasirandom(width = .15, groupOnX = T)),
  ggplot2::stat_summary(fun = median, colour = 'black', geom = 'text', size = 3, show.legend = FALSE, hjust = .5, vjust = -2.5, angle = 90)
)

geomJob_HalfHalfBox_Small <- list(
  gghalves::geom_half_boxplot(side = 'l', outlier.shape = NA, color = 'black'),
  gghalves::geom_half_point(side = 'r', color = 'black', alpha = .66, size = .75, shape = 21, transformation = ggbeeswarm::position_quasirandom(width = .075)),
  ggplot2::stat_summary(fun = min, colour = 'black', geom = 'text', size = 3, show.legend = FALSE, hjust = 1.25, angle = 90)
)

# Color palette(s) --------------------------------------------------------

colorPalette <- c(

  # Classes.
  'Poor Responder (≤100 days)' = '#8B0000', 
  'Ambiguous Responder (101-179 days)' = '#FFFFE0',
  'Good Responder (≥180 days)' = '#008080',
  'Ambiguous Prediction' = '#FFFFE0',
  'Predicted - Poor Responder' = '#8B0000',
  'Predicted - Good Responder' = '#008080',
  
  # Treatments
  'Enzalutamide' = '#2c7fb8', 
  'Abiraterone' = '#D93A28',
  
  # SV Classes.
  'Translocations' = '#375D96', 
  'Deletions' = '#ff8c00', 
  'Tandem Duplications' ='#fc6769', 
  'Insertions' = 'yellow', 
  'Inversions' = 'skyblue', 
  'Singles' = 'grey75',
  
  # Biopsy sites
  'Liver' = '#FF3500', 
  'Bone' = '#FEFEFE', 
  'Other' = '#4CA947', 
  'Lung' = '#9E4CD7', 
  'Lymph node' = '#0A6C94', 
  'Soft tissue' = '#EDAEAE'
  
)
