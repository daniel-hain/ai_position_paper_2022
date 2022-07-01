####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
##	Plotting
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

##################################################################
##	BEGIN: gg_color_hue():
##################################################################
# Little helper for teh standard Ggplot palette
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

##################################################################
##	BEGIN: gg_color_select():
##################################################################
# Little to make sure the colors of one or more plots are always aligned with a certain palette

gg_color_select <- function(x, pal = NULL) {
  require(RColorBrewer)

  n_col <- x %>% n_distinct() 
  if(is.null(pal)){mycol <- gg_color_hue(n_col) } 
  if(!is.null(pal)){mycol <- brewer.pal(n_col, pal)}
  return(mycol)
}
  
##################################################################
##	BEGIN: plot_timeline():
##################################################################
# Nicely plotting a timeline...
plot_summary_timeline <- function(x, y1, y2, t, by, 
                                  y1_text = "Number", 
                                  y2_text = "Share", 
                                  pal = NULL, 
                                  select_cat = NULL,  
                                  t_min = NULL, 
                                  t_max = NULL, 
                                  label = FALSE){
  require(RColorBrewer)
  require(directlabels)
  require(patchwork)
  
  x %>% arrange({{by}}, {{t}}) #%>% complete({{by}}, {{t}}, fill = list({{y1}} = 0, {{y2}} = 0)) %>%
  
  # select subcategories if necessary
  if(!is.null(select_cat)){ x %<>% filter({{by}}  %in% select_cat)}
  if(!is.null(t_min)){ x %<>% filter({{t}} >= t_min)}
  if(!is.null(t_max)){ x %<>% filter({{t}} <= PY_max)}
  
  # colors
  mycol <- x %>% pull({{by}}) %>% gg_color_select(pal = pal)
  #mycol <- x %>% gg_color_select(cat = {{by}}, pal = pal)
  
  # generate the plots
  x %<>% select({{t}}, {{y1}}, {{y2}}, {{by}}) 
  
  p1 <- x %>% ggplot(aes(x = {{t}},  y = {{y1}}, col = factor({{by}}) )) +
    geom_line(size = 1, alpha = 0.8, show.legend = FALSE) +
    labs(x = "Year", y = y1_text)  +
    scale_colour_manual(name = "Legend", values = mycol)
  
  p2 <- x %>% ggplot(aes(x = {{t}}, y = {{y2}}, fill = factor({{by}}))) +
    geom_area(position = "stack") +
    labs(x = "Year", y = y2_text)  +
    scale_fill_manual(name = "Legend", values = mycol)
  
  if(label == TRUE){
    p1 <- p1 + geom_dl(aes(label = {{by}}), method = list("last.bumpup", cex =0.75, hjust = 0.5, colour = "black"))
    p2 <- p2 +  geom_text(data = x %>% 
                            filter({{t}} == max({{t}})) %>% 
                            arrange(desc({{by}})) %>% 
                            mutate(pos = cumsum({{y2}}) - ({{y2}} / 2)),
                          aes(x = max({{t}}), 
                              y = pos, label =  {{by}}, 
                              size = 3) )
  }
  
  out <- p1 + p2 + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
  
  return(out)
}
