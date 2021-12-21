library(dplyr)
library(data.table)
library(progress)
library(data.table)
library(ranger)
library(pdp)
library(stringr)
library(igraph)
library(viridis)
library(ggplot2)
library(ggpubr)
#library(ggwordcloud)
library(ggheatmap) # CRAN, not the version on https://rdrr.io/github/tleonardi/ggheatmap/
library(aplot)

load('prepared datasets.Rdata')
if (!file.exists('outputs'))
{
  dir.create('outputs')
}


### HEATMAPS - states that are the least costly to reach
draw_heatmap <- function(matrix_list, i)
{
  cat('.')
  
  df = matrix_list[[i]]
  row.names(df) = 1:nrow(df)
  
  g <- ggheatmap(df, 
                 color=c('antiquewhite','seagreen'),
                 border='black',
                 cluster_rows  = TRUE,
                 cluster_cols = TRUE) %>%
    ggheatmap_theme(1,theme =list(
      theme(axis.text.x = element_text(angle = 90),
            axis.text.y = element_text(colour = "gray", size=0),
            legend.position='none')
    ))
  
  g = g %>% insert_top(
    ggplot() + theme_void() + ggtitle(names(matrix_list)[i]), 
    height = 0.05
  )
  
  return(g)
}

draw_heatmap_all = function(matrix_list)
{
  heatmaps = lapply(1:length(matrix_list), draw_heatmap, matrix_list=matrix_list)
  
  heatmaps_gridded = plot_list(gglist=heatmaps, nrow=2, ncol=ceiling(length(matrix_list)/2))
  
  return(heatmaps_gridded)
}

df_top_to = data_processed %>% 
  group_by(name) %>% 
  arrange(proportional_cost_improvement) %>%
  slice_tail(n=50) %>% select(name,to,n)

matrix_list_top_to = by(df_top_to, df_top_to$name, function(x) {
  m = matrix(data=0,nrow=nrow(x),ncol=max(x$n))
  
  for (i in 1:nrow(x))
  {
    states_this = as.numeric(strsplit(get_state_env(x$to[i])$state[1],split='\\*')[[1]])
    m[i,states_this] = 1
  }
  
  m = data.frame(m)
  names(m) = df_taxa$taxon[df_taxa$name %in% x$name] # this only works because the ordering is assumed
  return(m)
})

g_heatmap_top_to = draw_heatmap_all(matrix_list_top_to)
ggsave(g_heatmap_top_to, file='outputs/g_heatmap_top_to.png', width=10,height=5)







# PREDICTIVE MODELING



xvars <- c("richness_from", "richness_to","jaccard_vals", "net_species_gain", "name", "n", "m")

do_rf <- function(df)
{
  ranger(formula(sprintf("proportional_cost_improvement ~ %s", paste(xvars, collapse=" + "))), 
         data=df,
         importance='permutation',
         respect.unordered.factors = TRUE
  )
}

do_pdps <- function(m_rf, df_train)
{
  pdps = rbindlist(lapply(setdiff(xvars,c("name")), function(xvar) { 
    pdp_this = partial(m_rf, pred.var=c(xvar,"name"), 
                       train=df_train %>% sample_n(1000),
                       progress='text',
                       plot=FALSE,
                       chull=TRUE)
    pdp_this$predictor_name=xvar
    names(pdp_this)[1] = "xval"
    pdp_this = pdp_this %>% 
      as_tibble %>%
      mutate(xval = as.numeric(as.character(xval)), 
             name = as.character(name), 
             predictor_name = as.character(predictor_name))
    return(pdp_this)
  }))
  
  return(pdps)
}


m_rf_empirical = do_rf(data_empirical_processed)
pdps_empirical = do_pdps(m_rf_empirical, data_empirical_processed)



g_pdps_empirical = ggplot(pdps_empirical,aes(x=xval,y=yhat,col=name,group=name,shape=name)) + 
  geom_line() + 
  geom_point(size=2) +
  facet_wrap(~predictor_name,scales='free_x', labeller = labeller(predictor_name=
                                                                    c(jaccard_vals="Net Jaccard similarity", 
                                                                      n="Number of species in dataset (n)", 
                                                                      m="Number of environments in dataset (m)",
                                                                      richness_from="Starting state richness",
                                                                      richness_to="Target state richness",
                                                                      net_species_gain="Net richness change"))) +
  theme_bw() +
  ylab("Proportional cost improvement") +
  xlab("Predictor value") +
  labs(color='Dataset',shape='Dataset') +
  theme(legend.position='bottom') +
  scale_color_brewer(palette='Set1')

g_improvement_empirical = ggplot(data_empirical_processed, aes(x=proportional_cost_improvement,fill=name)) +
  geom_histogram(binwidth = 0.05) +
  theme_bw() + 
  facet_wrap(~name, scales='free_y') +
  theme(legend.position='none') + 
  xlab("Proportional cost improvement") +
  ylab("Number of state pairs") +
  scale_fill_brewer(palette='Set1') +
  scale_y_sqrt()

ggsave(ggarrange(g_improvement_empirical, g_pdps_empirical,nrow=2,ncol=1,labels=c('(a)','(b)'),heights=c(1,1.5)), 
       file='g_results_empirical.png',
       width=8,height=8)





# summarize findings
wilcox.test(data_empirical_processed$proportional_cost_improvement, 
            mu=0, 
            alternative='greater',
            paired=FALSE, 
            conf.int=TRUE)

pci_quantiles_empirical = data_empirical_processed %>% 
  group_by(name) %>%
  summarize(pci.25 = quantile(proportional_cost_improvement, 0.25), 
            pci.50 = quantile(proportional_cost_improvement, 0.50), 
            pci.75 = quantile(proportional_cost_improvement, 0.75), 
            pci.mean=mean(proportional_cost_improvement), 
            pci.sd=sd(proportional_cost_improvement),
            prob.improvement = length(which(proportional_cost_improvement > 0.1)) / length(proportional_cost_improvement))










plot_path <- function(df, n_env_max, n_actions_max, n_species_max)
{
  nr = nrow(df)
  el = as.matrix(data.frame(from=1:(nr-1),to=2:nr))
  g = graph_from_edgelist(el)
  V(g)$state = df$state
  V(g)$env = df$env
  E(g)$name = edge_labels[head(df$transition,-1)]
  
  colors = inlmisc::GetColors(n_env_max, start = 0.2, end = 0.9)
  
  plot(g, 
       layout=df %>% 
         select(id, richness) %>% 
         rename(x = id, y = richness) %>% 
         mutate(x = x/n_actions_max, y = y/n_species_max) %>% 
         as.matrix, 
       asp=1,
       rescale=FALSE,
       frame=TRUE,
       vertex.label=paste(df$state,"\n",df$env,sep=""),
       edge.arrow.size=0.1,
       edge.arrow.color='black',
       edge.label=E(g)$name,
       vertex.shape='none',
       vertex.label.color=colors[as.numeric(df$env)],
       vertex.label.family='sans') 
}

paths_cilate_ss = data_empirical_processed %>% 
  filter(name=="Ciliate (large)") %>%
  filter(proportional_cost_improvement > 0) %>%
  sample_n(10)

pdf(file='paths_ciliate.pdf',width=5,height=5)
for (i in 1:nrow(paths_cilate_ss))
{
  print(i)
  plot_path(get_path_df(paths_cilate_ss$path_astar[i]), n_env_max = max(data_empirical_processed$m), n_species_max = max(data_empirical_processed$n), n_actions_max = 10)
}
dev.off()



paths_human_gut_ss = data_empirical_processed %>% 
  mutate(num_natural = (net_cost_astar - floor(net_cost_astar))/0.1) %>%
  filter(name=="Human gut (n=12, m=1)") %>% 
  filter(num_natural >= 3)

pdf(file='best_paths_human_gut.pdf',width=5,height=5)
for (i in 1:nrow(paths_human_gut_ss))
{
  plot_path(get_path_df(paths_human_gut_ss$path_astar[i]))
}
dev.off()

