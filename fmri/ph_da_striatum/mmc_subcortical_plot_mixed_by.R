##########################
## AD's adaptation of Angela's Explore script
##########################

# NB: evt_time is shifted in MMC subcortical (1 should be 0)

# rm(list=ls())

library(wesanderson)
library(pracma)
library(ggplot2)
library(tidyverse)
library(viridis)
library(grid)
library(viridis) 
library(ggrepel)
library(scales)
library(data.table)

#Choose flags:
toalign <- "response" #options are clock or response
toprocess <- "subregion" # Options I want are atlas_value (every ROI and hemisphere separately), subregion (collapse across subregions and hemispheres), network (collapse across networks and hemispheres), structure (coarser groups, combined across hemispheres), structure_hipp (with anteiror and posterior hippocampus separated)
colorby <- "group" #options are subregion, group and structure
toplot <- "estimate" #options: are estimate and statistic
model <- "model_4" # The name of the model you want to plot (name defined in "Build_mixed_by.R" script)
# structures <- c("Caudate/ant. putamen","Amygdala","Hippocampus", "VS/Post. putamen", "Thalamus")
structures <- c("VTA")
region_to_run = "subcortical" #options: cortical (vmPFC), subcortical 
flip_stat = 0 #to flip statistic for plots (i.e. if model set up to look at attempters vs. depressed but you want to flip to see depressed vs. attempters)
split_hipp_amyg = TRUE #if want to separate out hippocampus and amygdala
exclude_prefeedback = 0 #set to 1 if you want to exclude the pre-feedback results from plots

#Set the working directory based on the options
subcortical_cache_dir = '~/OneDrive - University of Pittsburgh/Documents/SCEPTIC_fMRI/mmclock_subcortical_medusa/data'  # where the decons are
subcortical_output_dir = '~/OneDrive - University of Pittsburgh/Documents/SCEPTIC_fMRI/mmclock_subcortical_medusa/'  # where the decons are
map_dir = "~/code/schaefer_wb_parcellation/labels/"
# cortical_cache_dir = '/Users/angela/Documents/RESEARCH/ANALYSES/Medusa_analysis/cortical'
# vmpfc_hc_cache_dir = '/Users/angela/Documents/RESEARCH/ANALYSES/Medusa_analysis/vmPFC_hipp_interactions'

# Define where data is located and where plots will be saved
if(region_to_run == "subcortical") {
  if (toprocess == "structure") {dir <- paste0(subcortical_output_dir,'/Structure_CombinedHemis')
  dir.create(file.path(subcortical_output_dir, '/Structure_CombinedHemis'), showWarnings = FALSE)
  } else if (toprocess == "subregion") {dir <- paste0(subcortical_output_dir,'/Subregion_CombinedHemis')
  dir.create(file.path(subcortical_output_dir, '/Subregion_CombinedHemis'), showWarnings = FALSE)
  } else if (toprocess == "network") {dir <- paste0(subcortical_output_dir,'/Network_CombinedHemis')
  dir.create(file.path(subcortical_output_dir, '/Network_CombinedHemis'), showWarnings = FALSE)
  } else if (toprocess == "atlas_value") {dir <- paste0(subcortical_output_dir,'/Network_CombinedHemis')
  dir.create(file.path(subcortical_output_dir, '/Network_CombinedHemis'), showWarnings = FALSE)}
} else if(region_to_run == "cortical") {
  if (toprocess == "structure") {dir <- paste0(cortical_output_dir,'/Structure_CombinedHemis')
  dir.create(file.path(cortical_output_dir, '/Structure_CombinedHemis'), showWarnings = FALSE)
  } else if (toprocess == "subregion") {dir <- paste0(cortical_output_dir,'/Subregion_CombinedHemis')
  dir.create(file.path(cortical_output_dir, '/Subregion_CombinedHemis'), showWarnings = FALSE)
  } else if (toprocess == "network") {dir <- paste0(cortical_output_dir,'/Network_CombinedHemis')
  dir.create(file.path(cortical_output_dir, '/Network_CombinedHemis'), showWarnings = FALSE)
  } else if (toprocess == "atlas_value") {dir <- paste0(cortical_output_dir,'/Network_CombinedHemis')
  dir.create(file.path(cortical_output_dir, '/Network_CombinedHemis'), showWarnings = FALSE)}
} else if(toprocess == "vmpfc_hc") {dir <- vmpfc_hc_output_dir}
specific_output_dir <- dir
setwd(specific_output_dir)

#Load the Rdata file with mixed_by results
load(dir(full.names=T,pattern=(paste("*-",toalign,"-",model,"_with_VTA.Rdata", sep=""))))

if(toalign == "response"){epoch_label <- 'response'
} else if(toalign == "clock"){epoch_label <- 'clock'}

if(toplot == "statistic") {ylabel <- 't-statistic'
}else if(toplot == "estimate") {ylabel <- 'Coefficient (AU)'}

#Load processed labels including VTA as 445
setwd(subcortical_output_dir)
labels <- fread("Schaefer_445_region_labels_with_VTA.csv") %>% filter(atlas_value > 400)
# labels <- read_csv("Schaefer_444_region_labels.csv") %>% mutate(roi_num = as.numeric(roi_num)) %>% 
#   inner_join(read_csv("Schaefer_444_region_lookup.csv"), by = "roi_num") %>% mutate(atlas_value = as.numeric(roi_num)) 
# setwd(specific_output_dir)
# # label missing networks as amyg_hipp_thal (amygdala, hippocampus, thalamus)
# if (split_hipp_amyg){
#   labels$network[str_detect(labels$subregion, regex("Hippocampus", ignore_case=TRUE))] <- "hippocampus"
#   labels$network[str_detect(labels$subregion, regex("BLA", ignore_case=TRUE))] <- "amygdala"
#   labels$network[str_detect(labels$subregion, regex("CMN", ignore_case=TRUE))] <- "amygdala"
# } else {labels$network[is.na(labels$network)] <- "amyg_hipp_thal"}
# # #Merge labels with subcortical medusa data
# # # cortical_df <- cortical_df %>% inner_join(labels, by="atlas_value") #add labels to subcortical_df 
# # #Code the groupings
# labels <- labels %>% filter(atlas_value > 400) %>% mutate(structure = as.factor(case_when(
#   str_detect(subregion, regex("Anterior Putamen", ignore_case=TRUE)) ~ "Caudate/ant. putamen",
#   str_detect(subregion, regex("Caudate Tail and Lateral Putamen", ignore_case=TRUE)) ~ "Caudate/ant. putamen",
#   str_detect(subregion, regex("Caudate Head", ignore_case=TRUE)) ~ "Caudate/ant. putamen",
#   str_detect(subregion, regex("Ventral Striatum", ignore_case=TRUE)) ~ "VS/Post. putamen",
#   str_detect(subregion, regex("Posterior Putamen", ignore_case=TRUE)) ~ "VS/Post. putamen",
#   str_detect(subregion, regex("BLA", ignore_case=TRUE)) ~ "Amygdala",
#   str_detect(subregion, regex("CMN", ignore_case=TRUE)) ~ "Amygdala",
#   str_detect(subregion, regex("Hippocampus", ignore_case=TRUE)) ~ "Hippocampus",
#   str_detect(subregion, regex("Thalamus", ignore_case=TRUE)) ~ "Thalamus"))) %>%
#   mutate(hippocampus_subregion = as.factor(case_when( #Code hippocampus subregions - anterior vs. posterior
#     str_detect(subregion, regex("Anterior Putamen", ignore_case=TRUE)) ~ "Striatum",
#     str_detect(subregion, regex("Caudate Tail and Lateral Putamen", ignore_case=TRUE)) ~ "Striatum",
#     str_detect(subregion, regex("Caudate Head", ignore_case=TRUE)) ~ "Striatum",
#     str_detect(subregion, regex("Ventral Striatum", ignore_case=TRUE)) ~ "Striatum",
#     str_detect(subregion, regex("Posterior Putamen", ignore_case=TRUE)) ~ "Striatum",
#     str_detect(subregion, regex("BLA", ignore_case=TRUE)) ~ "Amygdala",
#     str_detect(subregion, regex("CMN", ignore_case=TRUE)) ~ "Amygdala",
#     str_detect(subregion, regex("Thalamus", ignore_case=TRUE)) ~ "Thalamus",
#     str_detect(subregion, regex("anterior", ignore_case=TRUE)) ~ "AH",
#     str_detect(subregion, regex("posterior", ignore_case=TRUE)) ~ "PH")))



if (toprocess == "atlas_value") {plot_df <- ddf$coef_df_reml %>% inner_join(labels, by="atlas_value")
} else if (toprocess == "subregion") {plot_df <- ddf$coef_df_reml %>% filter(effect == "fixed" & subregion != "Cerebellum") %>% select(!group) %>%
  inner_join(labels %>% select(subregion, group, network) %>% unique(), by="subregion")
} else if (toprocess == "structure") {plot_df <- ddf$coef_df_reml[,-7]
} else if (toprocess == "structure_hipp") {plot_df <- ddf$coef_df_reml[,-7]
} 

# ddf is the output structure; reml and ml are two different convergence methods (use reml - contains AIC and BIC andn log likelihood)
# ddf$coeff_df_reml - has all fixed and random effects; value is the estimate and estimated error is std.error
#consolidate the ddf to just the coef_df_reml
plot_df <- plot_df %>% filter(subregion != "Cerebellum") # drop cerebellum
#ddf$atlas_value <- as.factor(ddf$atlas_value)
#ddf$subregion <- as.factor(ddf$subregion)
plot_df$network <- as.factor(plot_df$network)

if (exclude_prefeedback) {
  plot_df<- plot_df %>% filter(evt_time>-0.1)
}

# filter structure[s] of interest
# plot_df <- plot_df %>% filter(structure %in% structures)
# 
# plot_df <- plot_df %>% mutate(p_level_fdr = as.factor(case_when(
#   padj_fdr_term_subregion > .05 ~ '1',
#   padj_fdr_term < .05 & padj_fdr_term > .01 ~ '2',
#   padj_fdr_term < .01 & padj_fdr_term > .001 ~ '3',
#   padj_fdr_term <.001 ~ '4')))

plot_df <- plot_df %>% mutate(p_level_fdr = as.factor(case_when(
  padj_fdr_term_subregion > .05 ~ '1',
  padj_fdr_term_subregion < .05 & padj_fdr_term_subregion > .01 ~ '2',
  padj_fdr_term_subregion < .01 & padj_fdr_term_subregion > .001 ~ '3',
  padj_fdr_term_subregion <.001 ~ '4')))

plot_df$p_level_fdr <- factor(plot_df$p_level_fdr, levels = c('1', '2', '3', '4'), labels = c("NS","p < .05", "p < .01", "p < .001"))
plot_df$`p, FDR-corrected` = plot_df$p_level_fdr

terms <- unique(plot_df$term[plot_df$effect=="fixed"])
plot_df$t <- plot_df$evt_time
##################
plot_df$t <- as.numeric(plot_df$t) - 1 ####### !!!!! correcting error in alignment
##################
#define color schemes
if (region_to_run == "subcortical") {if (split_hipp_amyg){
  color_scheme <- c('darkgoldenrod','#c912ab','brown4','darkgreen','#644bc8', 'darkred', 'darkslategrey', 'black')
  # color_scheme <- c('#66CCEE','#c912ab','#e17014','#32af64','#644bc8')
} else {color_scheme <- c('#32af64','#c912ab','#e17014','#644bc8')}
  background_color='gray70'
  major_grid_color='gray90' 
  minor_grid_color='gray80' 
  hline_color='white' 
  vline_color='white'
} else if (region_to_run == "cortical") {color_scheme <- c('#c912ab','#e17014','#644bc8') 
background_color='white' 
major_grid_color='gray80'
minor_grid_color='white'
hline_color='gray40'
vline_color='gray40'}

#If want to flip statistic for PLOTS
if (flip_stat) {plot_df$statistic <- -plot_df$statistic
plot_df$estimate <- -plot_df$estimate}

for (fe in terms) {
  message(paste0("\nPlotting ",fe,"\n"))
  edf <- plot_df %>% filter(term == paste(fe) & plot_df$effect=='fixed')
  termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
  if (region_to_run == "subcortical") {if (flip_stat) {fname = paste0(toalign,'-',fe,'-Plot-',toprocess,'-',model,'-Colorby-',colorby,"FLIPPED.pdf")}
    else if (structures == "Thalamus") {fname = paste0(paste0(structures, collapse = "_", "_"), toalign,'-',fe,'-Plot-',toprocess,'-',model,'-Colorby-',colorby,".pdf")
    }    else {fname = paste0(toalign,'-',fe,'-Plot-',toprocess,'-',model,'-Colorby-',colorby,".pdf")}
  } else if (region_to_run == "cortical") {if (flip_stat) {fname = paste0('vmPFC-',toalign,'-',fe,'-Plot-',toprocess,'-',model,'-Colorby-',colorby,"FLIPPED.pdf")}
    else {fname = paste0('vmPFC-',toalign,'-',fe,'-Plot-',toprocess,'-',model,'-Colorby-',colorby,".pdf")}}
  if (exclude_prefeedback) {
    if (region_to_run == "cortical") {pdf(fname, width = 4.5, height = 4)
    } else {pdf(fname, width = 5, height = 4)}
  } else {pdf(fname, width = 6, height = 5)}
  pd <- position_dodge(0.2)
  gg<-ggplot(edf, aes(x=t, y=.data[[toplot]], color=.data[[colorby]], group=subregion)) + #changed to flip y for groups compared to attempters 5/16 x=t, y=statistic, color=network)) + 
    geom_vline(xintercept = 0, lty = 'solid', color = vline_color, size = 1.5)+ xlab(epoch_label) + ylab('') +
    geom_point(aes(size=`p, FDR-corrected`, alpha=`p, FDR-corrected`), position=pd) + #show.legend = FALSE
    geom_label_repel(data = subset(edf, t == 5), aes(label = subregion, alpha = NULL, colour = .data[[colorby]] , x = t, y = .data[[toplot]]),
                     arrow = arrow(length = unit(0.01, "npc"), type = "closed"), xlim  = c(7,11), point.padding = NA, force = 3, size = 2, box.padding = 0.1)    +
    geom_label_repel(data = subset(edf, t == -3), aes(label = subregion, alpha = NULL, colour = .data[[colorby]] , x = t, y = .data[[toplot]]),
                     arrow = arrow(length = unit(0.01, "npc"), type = "closed"), xlim  = c(-10, -5), point.padding = NA, force = 2, size = 2, box.padding = 0.1)    +
    scale_y_continuous(expand = c(.5, 0)) +
    scale_colour_manual(values = color_scheme, guide = 'none') + 
    geom_line(aes(alpha = `p, FDR-corrected`), size=1.5, position=pd, show.legend = FALSE) + 
    geom_linerange(aes(ymin=estimate-std.error, ymax=estimate+std.error, alpha = `p, FDR-corrected`),show.legend = FALSE, position=pd) +
    scale_alpha_discrete(range=c(0.4, 1), drop=FALSE, guide = 'none') +
    scale_size_discrete(range=c(2,6), drop=FALSE, guide = 'none') # I think adding this last bit fixed my issue of the points getting resized when all factors weren't in the data
  # scale_alpha_discrete(guide = 'none') + 
  # scale_size_discrete(guide = 'none') # I think adding this last bit fixed my issue of the points getting resized when all factors weren't in the data
  gg <- gg  + theme_bw(base_size=13) +
    xlab(paste("time with respect to",toalign," (seconds)")) +
    ylab(ylabel) + xlim(c(min(plot_df$t) - 5, max(plot_df$t) + 5)) +
    theme(legend.title = element_blank(),
          panel.grid.major = element_line(colour = major_grid_color), 
          panel.grid.minor = element_line(colour = minor_grid_color), 
          panel.background = element_rect(fill = background_color),
          axis.title.y = element_text(margin=margin(r=6)),
          axis.title.x = element_text(margin=margin(t=6))#,
          # plot.margin = unit(c(1,3,1,1), "lines")
    )
  print(gg)
  dev.off()
}
