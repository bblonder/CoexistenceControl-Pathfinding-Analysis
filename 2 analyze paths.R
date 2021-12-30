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
library(ggheatmap) # CRAN, not the version on https://rdrr.io/github/tleonardi/ggheatmap/
library(aplot)
library(pals)
library(ggrepel)

# load in datasets 
load('prepared datasets.Rdata')
if (!file.exists('outputs'))
{
  dir.create('outputs')
}


int_breaks_rounded <- function(x, n = 5)  pretty(x, n)[round(pretty(x, n),1) %% 1 == 0]

# write out taxon names
write.csv(df_taxa %>% select(name, id, taxon),file='outputs/df_taxa_raw.csv',row.names=F)


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
  
  heatmaps_gridded = plot_list(gglist=heatmaps, nrow=3, ncol=ceiling(length(matrix_list)/3))
  
  return(heatmaps_gridded)
}

df_top_to = data_processed %>% 
  group_by(name) %>% 
  arrange(proportional_cost_improvement) %>%
  slice_tail(n=50) %>% select(name,to,n) %>%
  mutate(name = as.character(name))

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
ggsave(g_heatmap_top_to, file='outputs/g_heatmap_top_to.png', width=10,height=10)







# PREDICTIVE MODELING
trim_scale_df <- function(df, xvars, scale=TRUE, quantile_alpha=NULL, val_max=10)
{
  df_processed = df %>% 
                  group_by(name) %>%
                  mutate(across(all_of(xvars), function(x) {
                    if (!is.null(quantile_alpha))
                    {
                      q_lo = quantile(x,quantile_alpha/2,na.rm=T)
                      q_hi = quantile(x,1-quantile_alpha/2,na.rm=T)
                      x[x>q_hi] = NA
                      x[x<q_lo] = NA
                    }
                    else
                    {
                      x[x > val_max] = NA
                      x[x < -1*val_max] = NA
                    }
                    
                    if (scale==TRUE)
                    {
                      x = scale(x)
                    }
                    
                    return(as.numeric(x))
                    })) %>%
                  ungroup

  return(df_processed)
}

do_rf <- function(df, yvar, xvars)
{
  ranger(formula(sprintf("%s ~ %s", yvar, paste(xvars, collapse=" + "))), 
         data=df,
         #importance='permutation',
         verbose=TRUE,
         seed=1,
         respect.unordered.factors = TRUE,
         # tuning parameters
         mtry=3,
         num.trees=500
  )
}

do_pdps <- function(m_rf, df_train, xvars)
{
  df_train_this = df_train %>% sample_n(5000)
  pdps = rbindlist(lapply(setdiff(xvars,c("name")), function(xvar) {
    message(xvar)
    pdp_this = partial(m_rf, pred.var=c(xvar,"name"), 
                       train=df_train_this,
                       progress='text',
                       grid.resolution=10,
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

# trim and scale the dataset
data_processed = data_processed %>%
  mutate(delta_A = A_to_mean - A_from_mean) %>%
  mutate(delta_r = r_to_mean - r_from_mean)# %>%
# calculate a proportional length improvement (can be > 1 if the a* path is more circuitous, don't include naturals as they don't take time)
#  mutate(proportional_length_improvement = 1 - (path_length_astar - num_transitions_natural) / (path_length))

xvars <- c(#"richness_from", 
           #"richness_to",
           "jaccard_vals", 
           "net_species_gain", 
           "name", 
           "n", 
           "m", 
           #"A_from_mean", 
           #"r_from_mean", 
           #"A_to_mean", 
           #"r_to_mean",
           "delta_A", 
           "delta_r",
           "proportion_states")# 
xvars_to_trim = c(#"A_from_mean", "r_from_mean", "A_to_mean", "r_to_mean",
                  "delta_A", "delta_r")#c("A_from_mean", "r_from_mean", "A_to_mean", "r_to_mean") # #

data_processed_trimmed = trim_scale_df(data_processed, xvars_to_trim, scale = FALSE, val_max = 5)
data_processed_trimmed_melted = reshape2::melt(data_processed_trimmed[,xvars],id.vars='name')

# plot distributions
g_predictor_distribution = ggplot(data_processed_trimmed_melted, aes(x=value,color=name)) + 
  geom_density() + 
  facet_wrap(~variable,scales='free', labeller = labeller(variable=
                                                            c(jaccard_vals="Net Jaccard similarity", 
                                                              n="# species in dataset (n)", 
                                                              m="# environments in dataset (m)",
                                                              richness_from="Starting richness",
                                                              richness_to="Target richness",
                                                              net_species_gain="Net richness change",
                                                              A_from_mean="Starting A",
                                                              r_from_mean="Starting r",
                                                              A_to_mean="Target A",
                                                              r_to_mean="Target r",
                                                              delta_A = "Net A change",
                                                              delta_r = "Net r change",
                                                              proportion_states = "Frac. feasible+stable"))) +
  theme_bw() +
  scale_color_brewer(palette='Set1',name='Dataset') +
  xlab("Value") +
  ylab("Probability density") +
  scale_x_continuous(breaks = int_breaks_rounded, limits=function(x) { c(floor(min(x)), ceiling(max(x)))  })
ggsave(g_predictor_distribution, file='outputs/g_predictor_distribution.png',width=10,height=7)

# find rows with NA predictors
rows_na = which(apply(is.na(data_processed_trimmed[,xvars]),1,sum)>0)
data_processed_trimmed = data_processed_trimmed[-rows_na,]

# run random forest on the NA-dropped rows
# cost improvement
do_models = TRUE

if (do_models==TRUE)
{
  m_rf_pci = do_rf(data_processed_trimmed, 
               yvar = "proportional_cost_improvement", 
               xvars = xvars)
  
  pdps_pci = do_pdps(m_rf_pci, data_processed_trimmed, xvars)
  
  # length improvement
  # m_rf_pli = do_rf(data_processed_trimmed, 
  #                  yvar = "proportional_length_improvement", 
  #                  xvars = xvars)
  # 
  # pdps_pli = do_pdps(m_rf_pli, data_processed_trimmed, xvars)
  
  ## save PDPS and rf
  save(m_rf_pci, pdps_pci, 
       #m_rf_pci, pdps_pli, 
       file='models pdps.Rdata')
} else
{
  load('outputs/models pdps.Rdata')
}

plot_pdp <- function(pdps, ylab)
{
  g_pdps = ggplot(pdps,aes(x=xval,y=yhat,col=name,group=name,shape=name)) + 
    geom_line() + 
    geom_point(size=2) +
    facet_wrap(~predictor_name,scales='free_x', labeller = labeller(predictor_name=
                                                                      c(jaccard_vals="Net Jaccard similarity", 
                                                                        n="# species in dataset (n)", 
                                                                        m="# environments in dataset (m)",
                                                                        richness_from="Starting richness",
                                                                        richness_to="Target richness",
                                                                        net_species_gain="Net richness change",
                                                                        A_from_mean="Starting A",
                                                                        r_from_mean="Starting r",
                                                                        A_to_mean="Target A",
                                                                        r_to_mean="Target r",
                                                                        delta_A = "Net A change",
                                                                        delta_r = "Net r change",
                                                                        proportion_states = "Frac. feasible+stable"))) +
    theme_bw() +
    ylab(ylab) +
    xlab("Predictor value") +
    labs(color='Dataset',shape='Dataset') +
    theme(legend.position='bottom') +
    scale_color_brewer(palette='Set1') +
    #scale_color_manual(values = unname(alphabet())) +
    scale_x_continuous(breaks = int_breaks_rounded, limits=function(x) { c(floor(min(x)), ceiling(max(x))) } )
  
  return(g_pdps)
}

g_pdps_pci = plot_pdp(pdps_pci, "Proportional cost improvement")
ggsave(g_pdps_pci, file='outputs/g_pdps_pci.png',width=7.5,height=8)

# g_pdps_pli = plot_pdp(pdps_pli, "Proportional length improvement")
# ggsave(g_pdps_pli, file='outputs/g_pdps_pli.png',width=8,height=5)
# 

 





plot_histogram <- function(df, yvar, ylab)
{
  g_hist = ggplot(data_processed, aes_string(x=yvar,fill="name")) +
    geom_histogram(binwidth = 0.05) +
    theme_bw() + 
    facet_wrap(~name, scales='free_y') +
    theme(legend.position='none') + 
    xlab("Proportional cost improvement") +
    ylab("Number of state pairs") +
    scale_fill_brewer(palette='Set1') +
    scale_y_sqrt()
  
  return(g_hist)
}

ggsave(plot_histogram(data_processed, yvar="proportional_cost_improvement", ylab="Proportional cost improvement"), 
       file='outputs/g_histogram_pci.png',
       width=8,height=6)

# ggsave(plot_histogram(data_processed, yvar="proportional_length_improvement", ylab="Proportional length improvement"), 
#        file='outputs/g_histogram_pli.png',
#        width=8,height=6)



table_summaries = data_processed %>% 
  group_by(name) %>% 
  summarize(prob.pci.0.1=length(which(proportional_cost_improvement>0.1))/length(proportional_cost_improvement),
            pci.mean = mean(proportional_cost_improvement),
            pci.sd = sd(proportional_cost_improvement),
            #pli.mean = mean(proportional_length_improvement),
            #num_pairs_candidate = unique(num_states_candidate),
            proportion_fs = mean(proportion_states),
            n = unique(n),
            m = unique(m),
            )
table_summaries_to_print = table_summaries %>%
  select(name, 
         n,
         m,
         proportion_fs,
         prob.pci.0.1,
         pci.mean,
         #pli.mean,
         ) %>%
  mutate(across(proportion_fs:pci.mean, round, digits=2))
write.csv(table_summaries_to_print, file='outputs/table_summaries_to_print.csv',row.names=F)













plot_paths <- function(path_arg, name_arg, cost_max, richness_max, action_max, show_y_axis, title_text)
{
  costs = c(`=`=5,
            `>`=0.1,
            `+`=1,
            `-`=3)
  
  path_arg = as.character(path_arg)
  name_arg = as.character(name_arg)
  
  path_df_this = get_path_df(path_arg) %>%
    mutate(cost=costs[transition]) %>%
    mutate(total_cost = c(0,head(cumsum(cost),-1)))
  
  taxa_this = df_taxa %>% filter(name==name_arg)
  
  path_df_this$path_taxa = sapply(strsplit(path_df_this$state,"\\*"), function(taxon_id) { 
    taxa_names = taxa_this$taxon[as.numeric(as.character(taxon_id))]
    return(paste(taxa_names,collapse="\n"))
    })
  
  actions_this = data.frame(x=(head(path_df_this$id,-1)+tail(path_df_this$id,-1))/2,
                            y=(head(path_df_this$richness,-1)+tail(path_df_this$richness,-1))/2,
                            label=head(path_df_this$transition,-1))
  
  g_path_this = ggplot(path_df_this, aes(x=id, y=richness, label=path_taxa)) +
    geom_line(arrow = arrow(length=unit(0.30,"cm"), ends="last", type = "closed"),alpha=0.5,color='black') + 
    theme_bw() +
    geom_label(fill=gray(0.9,alpha = 0.5),size=1.5) +
    geom_label(data=actions_this,aes(x=x,y=y,label=label,fill=label),size=1.5) +
    xlab("Action order") +
    ylab("Richness") +
    scale_fill_manual(values=c(`+`='orange',`-`='red',`=`='lightgreen',`>`='lightblue')) +
    theme(legend.position='none') +
    scale_y_continuous(breaks = 0:(1+richness_max),limits=c(0,1+richness_max)) +
    scale_x_continuous(breaks = 1:(action_max+1),limits=c(1,(action_max+1)))
  
  if (!show_y_axis)
  {
    g_path_this = g_path_this + 
      theme(axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank())
  }
  else
  {
    g_path_this = g_path_this + 
                    geom_text(label=paste("\n",name_arg,collapse=""),x=1,y=Inf,hjust=0)
  }
  
  
  g_cost_this = ggplot(path_df_this, aes(x=id, y=total_cost, ymax=total_cost,ymin=0)) +
    geom_ribbon(alpha=0.25,fill='magenta') +
    geom_line(color='purple4') +
    theme_bw() +
    xlab("Action order") +
    ylab("Net cost") +
    ylim(0, cost_max) +
    scale_x_continuous(breaks = 1:(action_max+1),limits=c(1,(action_max+1))) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    ggtitle(title_text)
  
  if (!show_y_axis)
  {
    g_cost_this = g_cost_this + 
      theme(axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank())
  }
  
  return(plotlist=list(cost=g_cost_this, path=g_path_this))
}



plot_paths_comparison <- function(paths_all, i, save=TRUE)
{
  cost_max = max(paths_all$net_cost_astar[i],paths_all$net_cost[i]) 
  action_max = max(paths_all$path_length[i],paths_all$path_length_astar[i])
  richness_max = max(c(get_path_df(paths_all$path_astar[i])$richness,get_path_df(paths_all$path_nominal[i])$richness))
  
  g_astar = plot_paths(path_arg = paths_all$path_astar[i], 
               name_arg = paths_all$name[i],
               cost_max = cost_max,
               action_max = action_max,
               richness_max = richness_max,
               show_y_axis = TRUE,
               title_text="Optimal path")
  
  g_nominal = plot_paths(path_arg = paths_all$path_nominal[i], 
                         name_arg = paths_all$name[i],
                         cost_max = cost_max,
                         action_max = action_max,
                         richness_max = richness_max,
                         show_y_axis = FALSE,
                         title_text="Nominal path")
  
  g_final = ggarrange(g_astar$cost, g_nominal$cost, g_astar$path, g_nominal$path, 
            nrow=2,ncol=2,heights=c(1,3),
            align='v')
  
  if (save==TRUE)
  {
    ggsave(g_final, file=sprintf('outputs/paths_comparison_%d_%s.png',i,paths_all$name[i]),
         width=7,height=5)
  }
  
  cat('.')
  
  return(g_final)
}

set.seed(2)
high_pci_ids = data_processed %>% 
  mutate(row = 1:nrow(.)) %>%
  group_by(name) %>%
  arrange(proportional_cost_improvement) %>%
  filter(path_length_astar > 1 & proportional_cost_improvement > 0.1) %>%
  sample_n(15) %>%
  pull(row)

plots_paths_example = lapply(high_pci_ids, plot_paths_comparison, paths_all = data_processed, save=FALSE)
pdf(file='outputs/g_paths_examples.pdf',width=7,height=7)
for (i in 1:length(high_pci_ids))
{
  cat('+')
  print(plots_paths_example[[i]])
}
dev.off()

# find interesting cases

# identify which is Clodif 
id_clodif = df_taxa %>% 
  filter(name=="Mouse gut" & taxon=="Clodif") %>% 
  mutate(id = as.character(id)) %>%
  pull(id)

data_lost_clodif = data_processed %>% 
  mutate(row = 1:nrow(data_processed)) %>%
  filter(name=="Mouse gut") %>% 
  mutate(has_clodif_initial = grepl(id_clodif, from)) %>% 
  mutate(has_clodif_target = grepl(id_clodif, to)) %>% 
  filter(has_clodif_initial==TRUE & has_clodif_target==FALSE)

# this next part is a little slow! lots of regex...
do_regex = FALSE
if (do_regex==TRUE)
{
  pb = progress_bar$new(total = nrow(data_lost_clodif), format = "[:bar] :current :total")
  natural_losses = sapply(1:nrow(data_lost_clodif), function(i) {
    pb$tick()
    
    yesno = get_path_df(data_lost_clodif$path_astar[i]) %>% # get the i'th path
      mutate(has_clodif = grepl(id_clodif, state)) %>% # determine if it has clodif at each transition
      mutate(lost_clodif = c(as.logical(abs(diff(has_clodif))), NA)) %>% # find out when it lost clodif
      summarize(lost_clodif_natural = any(lost_clodif==TRUE & transition==">", na.rm=T)) # whether the time of loss is also a natural transition
    
    return(yesno)
  })
  save(natural_losses,file='clodif_natural_losses.Rdata')
} else
{
  load('clodif_natural_losses.Rdata')
}

data_lost_clodif$lost_clodif_natural_transition = unlist(natural_losses)

table_lost_clodif = data_lost_clodif %>% 
  group_by(lost_clodif_natural_transition) %>% 
  summarize(pci.mean = mean(proportional_cost_improvement), count=n()) %>%
  mutate(frac = count / sum(count))

sum(table_lost_clodif$count)
table_lost_clodif

rows_clodif = data_lost_clodif %>% filter(richness_to > 2 & 
                                            lost_clodif_natural_transition==TRUE & 
                                            richness_from > 3 & 
                                            proportional_cost_improvement > 0.4 &
                                            net_length_improvement > 0) %>% pull(row)

plot_paths_comparison(data_processed,i = rows_clodif[1114])





rows_ciliate = data_processed %>% 
  mutate(row = 1:nrow(.)) %>%
  filter(name=="Ciliate (m=3)") %>%
  filter(proportional_cost_improvement > 0.2 & num_transitions_environment==1 & richness_to > 0) %>% sample_n(20) %>% pull(row)

plot_paths_comparison(data_processed, i = rows_ciliate[7]) # 2, 5



# two interesting cases
g_paths_interesting = lapply(c(rows_clodif[1114], rows_ciliate[7]), plot_paths_comparison, paths_all = data_processed, save=FALSE)

ggsave(ggarrange(plotlist=g_paths_interesting, nrow=2,ncol=1,labels='auto'),width=7,height=7,file='outputs/g_paths_interesting.png')


