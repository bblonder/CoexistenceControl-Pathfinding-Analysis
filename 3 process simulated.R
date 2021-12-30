library(dplyr)
library(ggplot2)
library(data.table)

files_simulated = dir(path='../CoexistenceControl-Private/data/results/synthetic',
                      pattern='*csv',
                      full.names = TRUE,
                      recursive = TRUE)

process_simulated <- function(fn)
{
   frac_fs = readLines(gsub("csv$","txt",fn))[17] %>% 
     strsplit(split=": ") %>% 
     unlist %>% 
     tail(1) %>% 
     as.numeric
  # 
  # num_states_total = readLines(gsub("csv$","txt",fn))[19] %>% 
  #   strsplit(split=": ") %>% 
  #   unlist %>% 
  #   tail(1) %>% 
  #   as.numeric
  
  num_species = readLines(gsub("csv$","txt",fn))[3] %>% 
    strsplit(split=": ") %>% 
    unlist %>% 
    tail(1) %>% 
    as.numeric
  
  num_env = readLines(gsub("csv$","txt",fn))[4] %>% 
    strsplit(split=": ") %>% 
    unlist %>% 
    tail(1) %>% 
    as.numeric
  
  distribution_A = readLines(gsub("csv$","txt",fn))[10] %>% 
    gsub(pattern='A parameters: Normal{Float64}(μ=',replacement="",fixed=TRUE) %>% 
    gsub(pattern=' σ=',replacement="",fixed=TRUE) %>% 
    gsub(pattern=")",replacement="",fixed=TRUE) %>%
    strsplit(split=",") %>%
    unlist %>%
    as.numeric
  names(distribution_A) = c("mu.A","sigma.A")
  
  distribution_r = readLines(gsub("csv$","txt",fn))[11] %>% 
    gsub(pattern='r parameters: Normal{Float64}(μ=',replacement="",fixed=TRUE) %>% 
    gsub(pattern=' σ=',replacement="",fixed=TRUE) %>% 
    gsub(pattern=")",replacement="",fixed=TRUE) %>%
    strsplit(split=",") %>%
    unlist %>%
    as.numeric
  names(distribution_r) = c("mu.r","sigma.r")
  # 
  # df_astar_results = read.csv(fn)
  # cost_improvement_stats = df_astar_results %>%
  #                           mutate(proportional_cost_improvement = net_cost_improvement / net_cost) %>%
  #                           summarize(mean_pci = mean(proportional_cost_improvement), mean_ci = mean(net_cost_improvement))
  # 
  
  
  df_stats = data.frame(#fraction_fs = num_states_fs / num_states_total, 
                        #cost_improvement_mean = cost_improvement_stats$mean_ci,
                        #proportional_cost_improvement_mean = cost_improvement_stats$mean_pci,
                        frac_fs,
                        num_species, 
                        num_env,
                        mu.A=distribution_A["mu.A"], sigma.A=distribution_A["sigma.A"],
                        mu.r=distribution_r["mu.r"], sigma.r=distribution_r["sigma.r"])
  
  return(df_stats)
}

stats_all = rbindlist(lapply(files_simulated,process_simulated)) %>%
  mutate(name = sprintf("n=%d m=%d", num_species, num_env))

table_glv_stats = stats_all %>% 
  group_by(name) %>%
  summarize(across(mu.A:sigma.r, mean))

write.csv(table_glv_stats, file='outputs/table_glv_stats.csv', row.names = FALSE)



stats_all_for_plotting = stats_all %>% select(frac_fs,  name)

g_fs = ggplot(stats_all_for_plotting, aes(x=frac_fs,color=name,fill=name)) +
  geom_density(alpha=0.5) +
  theme_bw() +
  theme(legend.position='bottom') +
  xlab("Fraction of states feasible+stable") +
  ylab("Probability density") +
  scale_color_brewer(name='Dataset',palette='Set1') +
  scale_fill_brewer(name='Dataset',palette='Set1')

ggsave(g_fs, file='outputs/g_fs.png',width=6,height=4)
