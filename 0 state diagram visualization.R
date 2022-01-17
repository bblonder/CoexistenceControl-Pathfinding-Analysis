library(igraph)
library(viridis)
library(ggplot2)
library(ggrepel)
library(randomForest)
library(dplyr)
library(tidyr)

if (!file.exists('outputs_r'))
{
  dir.create('outputs_r')
}

generate_assemblages <- function(n, labels=letters)
{
  df = expand.grid(replicate(n, 0:1, simplify = FALSE))
  names(df) = labels[1:n]
  
  new_states = as.data.frame(matrix(data=0,nrow=nrow(df),ncol=ncol(df)))
  names(new_states) = paste(names(df),"new",sep=".")
  
  df_final = data.frame(df, 
                        stable=NA, 
                        feasible=NA,
                        tau=NA,
                        richness=NA, 
                        abundance.mean.resident=NA, 
                        abundance.sd.resident=NA, 
                        new_states)
  df_final$richness = apply(df_final[,1:n],1,sum)
  
  return(df_final)
}

generate_params <- function(n, labels=letters, A.diag=-1,
                            distribution.A='norm',
                            A.norm.mean=-0.2, A.norm.sd=0.2,  
                            A.lnorm_neg.meanlog=log(0.25), A.lnorm_neg.sdlog=log(1.5),
                            A.lnorm_pos.meanlog=log(0.25), A.lnorm_pos.sdlog=log(1.5),
                            distribution.r='norm',
                            r.norm.mean=1, r.norm.sd=0.5, 
                            r.lnorm_pos.meanlog=log(1), 
                            r.lnorm_pos.sdlog=log(2))
{
  # interaction matrix
  A = matrix(data=NA,nrow=n,ncol=n)
  dimnames(A) = list(labels[1:n],labels[1:n])
  
  # intrinsic growth rate vector
  r = rep(NA, n)
  names(r) = labels[1:n]
  
  # fill A
  if (distribution.A=='norm')
  {
    A[] = rnorm(n=n*n,mean=A.norm.mean,sd=A.norm.sd) 
  }
  else if (distribution.A=='lnorm_neg')
  {
    A[] = -1*rlnorm(n=n*n,meanlog=A.lnorm_neg.meanlog,sdlog=A.lnorm_neg.sdlog) 
  }
  else if (distribution.A=='lnorm_pos')
  {
    A[] = rlnorm(n=n*n,meanlog=A.lnorm_pos.meanlog,sdlog=A.lnorm_pos.sdlog)
  }
  else
  {
    stop('distribution not found')
  }
  
  # fill r
  if (distribution.r=='norm')
  {
    r[] = rnorm(n=n,mean=r.norm.mean,sd=r.norm.sd)
  }
  else if (distribution.r=='lnorm_pos')
  {
    r[] = rlnorm(n=n,meanlog=r.lnorm_pos.meanlog,sdlog=r.lnorm_pos.sdlog)
  }
  else
  {
    stop('distribution not found')
  }
  
  # set diagonal for A if requested
  if (!is.null(A.diag))
  {
    diag(A) <- rep(A.diag,n)
  }
  
  # return combo
  return(list(A=A,r=r))
}


assign_params <- function(assemblage, params)
{
  stopifnot(nrow(assemblage)==1)
  # find species that are present
  species_indices_present = which(as.numeric(assemblage)==1)
  
  # pick subset of parameters (assuming that the parameters don't change when subsetting)
  A_this = params$A[species_indices_present, species_indices_present, drop=FALSE]
  r_this = params$r[species_indices_present]
  
  return(list(A=A_this,r=r_this))
}

determine_feasibility <- function(params)
{
  if (length(params$r) > 0)
  {
    x = -1 * solve(params$A) %*% params$r
    
    feasibility = all(x > 0)
    new_state = names(params$r)[which(x > 0)]
    if (length(new_state)==0)
    {
      new_state = NULL
    }
    
    abundance.mean.resident = mean(x[x>0])
    abundance.sd.resident = sd(x[x>0])
  }
  else
  {
    feasibility = TRUE
    new_state = NULL
    abundance.mean.resident = NA
    abundance.sd.resident = NA
  }
  
  return(list(feasible=feasibility, 
              new_state = new_state, 
              abundance.mean.resident = abundance.mean.resident, 
              abundance.sd.resident = abundance.mean.resident))
}

determine_stability <- function(params)
{
  if (length(params$r) > 0)
  {
    print(str(params$r))
    print(str(params$A))
    lambda = eigen( diag(params$r,nrow=length(params$r)) %*% params$A  )$values
    stability = all(Re(lambda) < 0)
    if (stability == TRUE)
    {
      tau = -1/max(Re(lambda)) # time constant is inverse of the biggest eigenvalue
    }
    else
    {
      tau = 0
    }
  }
  else
  {
    lambda = NA
    stability = TRUE
    tau = NA
  }
  
  return(list(stable=stability,lambda=lambda,tau=tau))
}

generate_assemblage_transitions <- function(assemblages, params)
{
  n = log2(nrow(assemblages))
  for (i in 1:nrow(assemblages)) # skip the no-species scenario
  {
    params_this = assign_params(assemblage = assemblages[i,1:n], params = params)
    
    f = determine_feasibility(params_this)
    s = determine_stability(params_this)
    
    assemblages[i,"stable"] = s$stable
    assemblages[i,"feasible"] = f$feasible
    assemblages[i,"tau"] = s$tau
    assemblages[i,"abundance.mean.resident"] = f$abundance.mean.resident
    assemblages[i,"abundance.sd.resident"] = f$abundance.sd.resident
    
    if (!is.null(f$new_state))
    {
      assemblages[i,paste(f$new_state,"new",sep=".")] = 1
    }
  }
  return(assemblages)
}


classify_reachable_vertices <- function(network)
{
  vertices_all = V(network)
  
  reachable = rep(NA, length(vertices_all))
  
  for (i in 1:length(vertices_all))
  {
    # if there is at least one natural link pointing in to this verte
    reachable[i] = any(incident(network,v=V(network)[i],mode='in')$type=="natural")
  }
  return(reachable)
}

make_network <- function(assemblages, include.natural = TRUE)
{
  n = log2(nrow(assemblage_transitions))
  # generate labels for natural transitions
  assemblages$from = apply(assemblages[,1:n], 1, function(x) { paste(names(x)[which(x>0)],collapse=" ") })
  # weird numbers in columns help us pick the right columns for the dataframe created in generate_assemblages()
  assemblages$to = apply(assemblages[,(n+7):(2*n+6)], 1, function(x) { paste(gsub("\\.new$","",names(x))[which(x>0)],collapse=" ") })
  
  # identify all single-species additions/deletions  
  transitions_single_matrix = (as.matrix(dist(assemblages_all[,1:n],method='manhattan'))==1)
  dimnames(transitions_single_matrix) = list(assemblages$from, assemblages$from)
  
  transitions_single_df = which(transitions_single_matrix==TRUE,arr.ind=TRUE)
  row.names(transitions_single_df) = NULL
  transitions_single_df = as.data.frame(transitions_single_df)
  # assign labels
  transitions_single_df$from = assemblages$from[ transitions_single_df$row ]
  transitions_single_df$to = assemblages$from[ transitions_single_df$col ]
  # determine if species addition or removal
  transitions_single_df$type = ifelse(sign(sapply(transitions_single_df$to, nchar) - sapply(transitions_single_df$from, nchar))>0,"addition","deletion")
  # keep necessary columns
  transitions_single_df = transitions_single_df[,c("from","to","type")]
  
  # add natural relaxations of unstable states
  transitions_natural_df = assemblages[,c("from","to")]
  transitions_natural_df$type = "natural"
  
  # put all the transition types together
  if (include.natural == TRUE)
  {
    transitions_all = rbind(transitions_single_df, transitions_natural_df)
  }
  else
  {
    transitions_all = transitions_single_df
  }
  
  # simplify the transitions
  # remove loops
  transitions_all = transitions_all[transitions_all$from != transitions_all$to, ]
  
  g = graph_from_data_frame(transitions_all,directed=TRUE,vertices=assemblages[,c("from","stable","feasible","richness")])
  
  # determine reachable vertices
  reachable = classify_reachable_vertices(g)
  V(g)$reachable = reachable
  
  # simplify the transitions, keeping only the minimum cost links when there are multiple choices, and removing loops
  #g = simplify(g, remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = 'min')
  
  return(list(graph=g,transitions=transitions_all)) # consider also returning the transition table
}


layout_network <- function(graph)
{
  vertex_names = V(graph)$name
  
  y = sapply(vertex_names, nchar)
  
  df = data.frame(name=vertex_names,x=NA,y=y)
  #df = df[order(df$y),]
  
  for (i in unique(df$y))
  {
    rows_this = which(df$y==i)
    x_out = seq(0,1,length.out=length(rows_this))
    x_out = x_out - mean(x_out)
    df$x[rows_this] = x_out
  }
  
  # select key columns
  df = as.matrix(df[,c("x","y")])
  
  return(df)
}


plot_transitions <- function(network, edge_highlight=NULL, vertex_highlight_start=NULL, vertex_highlight_end=NULL)
{
  V(network)$name = gsub(" ","",V(network)$name)
  edge_types = E(network)$type
  
  edge_colors = c(natural=rgb(0,0,1),addition=rgb(1,0.5,0,0.5),deletion=rgb(1,0,0,0.5))[edge_types]
  edge_widths = rep(0.75, length(edge_types))
  edge_arrow_sizes = rep(0.5, length(edge_types))
  if(!is.null(edge_highlight))
  {
    edge_widths[edge_highlight] = 6
    #edge_arrow_sizes[edge_highlight] = 10 # not currently implemented in igraph...
  }
  
  vertex_colors = ifelse(V(network)$feasible==1 & V(network)$stable==1,'white','gray')
  vertex_frame_colors = 'black'
  
  vertex_shapes = rep('rectangle',length(V(network)))           
  # if(!is.null(vertex_highlight_start))
  # {
  #   vertex_shapes[ which(V(network)$name==vertex_highlight_start) ] = 'circle'
  # }
  # if(!is.null(vertex_highlight_end))
  # {
  #   vertex_shapes[ which(V(network)$name==vertex_highlight_end) ] = 'circle'
  # }
  
  plot(network,
       asp=0,
       layout=layout_network(network),
       edge.arrow.size=edge_arrow_sizes, 
       edge.label=NA,
       edge.curved=TRUE,
       edge.color=edge_colors,
       edge.width=edge_widths,
       vertex.frame.color=vertex_frame_colors,
       vertex.shape=vertex_shapes,
       vertex.label.family='sans',
       vertex.size=10,
       vertex.size2=10,
       vertex.label.cex=0.7,
       vertex.label.color='black',
       vertex.color=vertex_colors)
}

# show an example network
set.seed(1) # replicability
nsp <- 6
params_all = generate_params(n = nsp, A.norm.mean = -0.7, r.norm.mean = 1.7)
assemblages_all = generate_assemblages(n = nsp)
assemblage_transitions = generate_assemblage_transitions(params = params_all, assemblages = assemblages_all)
network_model = make_network(assemblage_transitions)
network_model_no_natural = make_network(assemblage_transitions, include.natural = FALSE)






# identify subset of vertices that are stable/feasible endpoints
candidate_vertices = V(network_model$graph)$name[V(network_model$graph)$feasible & V(network_model$graph)$stable]
# count # of species in each endpoint
candidate_n_species = sapply(sapply(candidate_vertices, strsplit," ",fixed=TRUE),length)
# calculate species gain in each transition
candidate_net_species_gain = -1*outer(candidate_n_species, candidate_n_species, FUN="-")
candidate_endpoints = outer(candidate_vertices, candidate_vertices, paste, sep="->")
candidate_from = outer(candidate_vertices, candidate_vertices,function(a,b){a})
candidate_to = outer(candidate_vertices, candidate_vertices,function(a,b){b})

# calculate transition cost with/without natural paths
weight_table = c(addition=1,deletion=3,natural=0)
candidate_net_cost = distances(network_model$graph, 
                               v=candidate_vertices, 
                               to=candidate_vertices,
                               mode='out',
                               weights=weight_table[E(network_model$graph)$type])

candidate_net_cost_no_natural = distances(network_model_no_natural$graph, 
                                          v=candidate_vertices, 
                                          to=candidate_vertices,
                                          mode='out',
                                          weights=weight_table[E(network_model_no_natural$graph)$type])

candidate_net_cost_improvement = candidate_net_cost_no_natural - candidate_net_cost
# make a table
candidate_table = data.frame(paths = as.vector(candidate_endpoints), 
                             from = as.vector(candidate_from),
                             to = as.vector(candidate_to),
                             net_species_gain = as.vector(candidate_net_species_gain),
                             net_cost = as.vector(candidate_net_cost), 
                             net_cost_no_natural = as.vector(candidate_net_cost_no_natural),
                             net_cost_improvement = as.vector(candidate_net_cost_improvement))
# make space for the best edge and vertex paths
candidate_table$epath = NA
candidate_table$vpath = NA
candidate_table$epath_no_natural = NA
candidate_table$vpath_no_natural = NA

# count whether we are doing a net species addition/removal
add_path_info <- function(network, from, to)
{
  sp = shortest_paths(graph=network, from=as.character(from), to=as.character(to), mode='out',weights=NULL,output='both')
  epath = paste(as.numeric(sp$epath[[1]]),collapse="-")
  vpath = paste(V(network)[as.numeric(sp$vpath[[1]])]$name,collapse="-")
  return(c(epath=epath,vpath=vpath))
}

for (i in 1:nrow(candidate_table))
{
  path_this = add_path_info(network=network_model$graph,from=candidate_table$from[i],to=candidate_table$to[i])
  path_this_no_natural = add_path_info(network=network_model_no_natural$graph,from=candidate_table$from[i],to=candidate_table$to[i])
  
  candidate_table$epath[i] = path_this["epath"]
  candidate_table$vpath[i] = path_this["vpath"]
  candidate_table$epath_no_natural[i] = path_this_no_natural["epath"]
  candidate_table$vpath_no_natural[i] = path_this_no_natural["vpath"]
  print(i/nrow(candidate_table))
}

best_improvements  = candidate_table[candidate_table$net_cost_improvement == 6,]

for (i in 1:1)#nrow(best_improvements))
{
  # find the best path
  epath_this = as.numeric(strsplit(best_improvements[i,"epath"],"-")[[1]])
  epath_types = E(network_model$graph)$type[epath_this]
  
  epath_this_no_natural = as.numeric(strsplit(best_improvements[i,"epath_no_natural"],"-")[[1]])
  epath_types_no_natural = E(network_model_no_natural$graph)$type[epath_this_no_natural]
  
  pdf(file=sprintf('outputs_r/g_lowest_cost_%s-%s.pdf',best_improvements[i,"from"],best_improvements[i,"to"]),width=1.5*nsp,height=1*nsp)
  plot_transitions(network_model$graph, 
                   edge_highlight = NULL, 
                   vertex_highlight_start = NULL, 
                   vertex_highlight_end = NULL)
  plot_transitions(network_model$graph, 
                   edge_highlight = epath_this, 
                   vertex_highlight_start = best_improvements[i,"from"], 
                   vertex_highlight_end = best_improvements[i,"to"])
  plot_transitions(network_model_no_natural$graph, 
                   edge_highlight = epath_this_no_natural, 
                   vertex_highlight_start = best_improvements[i,"from"], 
                   vertex_highlight_end = best_improvements[i,"to"])
  dev.off()
  
  png(file=sprintf('outputs_r/g_lowest_cost_%s-%s.png',best_improvements[i,"from"],best_improvements[i,"to"]),width=1.5*nsp,height=1*nsp,units='in',res=300)
  plot_transitions(network_model$graph, 
                   edge_highlight = NULL, 
                   vertex_highlight_start = NULL, 
                   vertex_highlight_end = NULL)
  dev.off()
  
  print(i/nrow(best_improvements))
}
