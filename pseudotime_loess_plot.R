# 
# setwd("/Volumes/bgudenas/Proj/PB_origins/CB_PN_RT_Atlas/src")
# 
# ps_mat = readRDS(out_df)
# gs_mat = readRDS(out_gs)
# gs_mat = gs_mat[ ,c("Rb_model_TFs11", "PB_G3_RETB10")]
# gs_mat[gs_mat < 0 ] = 0
# 
# pseudotime_vec = ps_mat$slingPseudotime_2
# celltype_vec = ps_mat$Celltype
# timepoint_vec = ps_mat$Timepoint
#saveRDS(df, "../Data/DiffusionMap/Pseudotime_Loess_Pineal_testdf.rds")

Pseudotime_GS_Smooth = function(gs_mat,
                                pseudotime_vec,
                                celltype_vec,
                                timepoint_vec,
                                out_prefix,
                                col_pal = NULL,
                                data_dir
                                ){
  library(ggplot2)
  library(viridis)
  th = theme(text = element_text(size = 10, face = "bold"),
             axis.text.x = element_text(angle = 45, hjust = 1),
             plot.title = element_text(size = 14, hjust = 0.5 ),
             axis.title.x = element_text(size = 10, face = "bold"),
             axis.title.y = element_text(size = 10, face = "bold")) 
  
  stopifnot( nrow(gs_mat) == length(pseudotime_vec))
  stopifnot( length(celltype_vec) == length(pseudotime_vec))
  
  df = data.frame("Barcode" = rownames(gs_mat),
                  "Pseudotime" = pseudotime_vec,
                  "Celltype" = celltype_vec,
                  "Timepoint" = timepoint_vec)
  
  gs_mat = gs_mat[!is.na(df$Pseudotime), ]
  df = df[!is.na(df$Pseudotime), ]
 # saveRDS(df, "../Data/DiffusionMap/Pseudotime_Loess_Pineal_testdf.rds") ## for tesint pseudotime block finder
  for (i in colnames(gs_mat)){
    gs_mat[ ,i] = minmax(gs_mat[ ,i])
  }
  df = cbind(df, gs_mat)
  
  sel_col = colnames(gs_mat)[ncol(gs_mat)]
  minq = 0.025
  maxq = 0.975
  gs_vec = df[ ,sel_col]
  min_bound = quantile(df$Pseudotime[gs_vec > 0], minq )
  max_bound = quantile(df$Pseudotime[gs_vec > 0], maxq )
  
  out_res = df[(df$Pseudotime >= min_bound) & (df$Pseudotime <= max_bound),  ]
  out_fp = paste0(data_dir, "Dev_Window_", out_prefix, ".csv")
  write.csv(out_res, out_fp, quote = FALSE )
  
# Create plot just to find Ymax of first GS -------------------------------
  g1 = ggplot(df, aes(x = Pseudotime,
                      y = .data[[sel_col]], colour = "black")) +
    geom_smooth( span=0.5, method = "loess", se = FALSE, colour = "black") +
    theme_bw() +
    th
  ymax = max(ggplot_build(g1)$layout$panel_params[[1]]$y.range)*1.5
  
  g1 = ggplot(df, aes(x = Pseudotime,
                      y = .data[[sel_col]], colour = "black")) +
    annotate(geom = "rect", xmin=min_bound, xmax=max_bound, ymin=0, ymax=ymax, alpha = 0.4,  fill = "grey") +
   # geom_rect(inherit.aes = FALSE, aes(xmin=min_bound, xmax=max_bound, ymin=0, ymax=ymax), color= "grey", fill = "grey", alpha = 0.5) +
    geom_smooth( span=0.5, method = "loess", se = FALSE, colour = "black") +
    theme_bw() +
    th
  
  gs_names = colnames(gs_mat)
  colos = rainbow(length(gs_names))
  colos = c("purple3","#EBBC2E")
  names(colos) = gs_names
  
  for ( i in 1:length(gs_names)){
    gs_sel = gs_names[i]
    colo = colos[i]
    g1 = g1 +
      geom_smooth(data = df, aes(x = Pseudotime, y = .data[[gs_sel]]), lwd=2, span=0.5, method = "loess", se = FALSE, color = colo)
  }
  
  ymax = max(ggplot_build(g1)$layout$panel_params[[1]]$y.range)
  g1 = g1  + coord_cartesian(ylim = c(0, ymax)) +
    ylab("Enrichment Score") +
    theme(axis.title.x =  element_blank(),
          panel.grid = element_blank())
  
  ## dummy variable for y axis to prevent points overlapping (jitter with height based on density)
  density_vec <- density(df$Pseudotime)
  desired_length <- length(df$Pseudotime)
  density_interp <- approx(density_vec$x, density_vec$y, n = desired_length)
  density_vector <- (density_interp$y)
  ## set half to negative for random spread
  set.seed(42)
  
  density_vector[seq(1,length(density_vector), 2)] = density_vector[seq(1,length(density_vector), 2)] *-1
  df$Dummy = density_vector
  df$Dummy =1
  
  ## celltype  scatter
  g2 = ggplot(df, aes(x =Pseudotime, y = Dummy )) +
   # geom_jitter(height = 0, width = 0, aes(color = Celltype ), size = 0.4, alpha = 0.3) +
    ggbeeswarm::geom_quasirandom(groupOnX = FALSE, aes(color = Celltype ), size = 0.25, alpha = 0.66, width = 0.6) +
    theme_bw() +
    th +
    ylab("") +
    guides(color = guide_legend(override.aes = list(size = 2))) +
   # coord_cartesian(ylim = c(0.5, 1.5)) +
    theme(axis.text.y =  element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  if (!is.null(col_pal)){g2 = g2 + scale_color_brewer(palette = col_pal)}
  
  ## timepoint scatter
  # g3 = ggplot(df, aes(x = Pseudotime, y = Dummy )) +
  #   geom_jitter(height = 0.1, width = 0, aes(color = Timepoint ), size = 0.4, alpha = 0.3) +
  #  # scale_color_brewer(palette = "YlOrRd") +
  #   scale_color_viridis(discrete = TRUE) +
  #   theme_bw() +
  #   th +
  #   ylab("") +
  #   guides(color = guide_legend(override.aes = list(size = 2))) +
  #   coord_cartesian(ylim = c(0.5, 1.5)) +
  #   theme(axis.text.y =  element_blank(),
  #         axis.ticks.y = element_blank(),
  #         panel.border = element_blank(),
  #         panel.grid.major = element_blank(),
  #         panel.grid.minor = element_blank())
  
  ldf = reshape2::melt(df[ ,-c(2,4)], "Celltype")
  colos = c(colos, "black")
  names(colos)[length(colos)] = sel_col
  
  glegend = ggplot(ldf, aes(x = Celltype,
                       y = variable)) +
    geom_line(aes(color = variable)) +
    scale_color_manual(values = colos ) +
    theme_bw() +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.y = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 2)))
  
  out_plot = paste0(data_dir, "legend_", out_prefix, ".pdf")
  ggsave(glegend, device = "pdf", filename = out_plot, width = 10, height = 10, dpi = 300)
  #g4 = cowplot::plot_grid(plotlist = list(g1,g2,g3), nrow = 3, rel_heights = c(5,1,1), align = "v", axis = "lr")
  g4 = cowplot::plot_grid(plotlist = list(g1,g2), nrow = 2,ncol=1, rel_heights = c(2.5,1), align = "v", axis = "lr")
  
  out_plot = paste0(data_dir, "pseudotime_loess_", out_prefix, ".pdf")
  ggsave(g4, device = "pdf", filename = out_plot, dpi = 300, width = 6, height = 6, bg="white")
  
}



# -------------------------------------------------------------------------
minmax = function(nums){
  output = (nums - min(nums))/(max(nums) - min(nums)) 
}
