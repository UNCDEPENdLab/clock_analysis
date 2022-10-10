plot_medusa <- function(coef_obj, x="evt_time", y="estimate", ymin=NULL, ymax=NULL, color=NULL, facet_by=NULL, 
                        p.value="p.value", lty_by=NULL, pdf_by=c("term", "model_name"), panel_by=NULL, out_dir=getwd(), 
                        plot_type="line", flip=FALSE, term_filter=NULL, width=9, height=7, include_title=TRUE) {

  require(patchwork)  
  require(dplyr)
  require(ggplot2)
  require(data.table)
  
  colors <- RColorBrewer::brewer.pal(4, "Dark2") %>% setNames(c("1" = "MT+","2" = "Premotor","3" = "Rostral PPC","4" = "Caudal PPC"))
  
  
  if (is.data.frame(coef_obj)) {
    to_plot <- coef_obj
  } else if (checkmate::test_data_frame(coef_obj$coef_df_reml)) {
    to_plot <- coef_obj$coef_df_reml
  } else if (checkmate::test_data_frame(coef_obj$coef_df_ml)) {
    to_plot <- coef_obj$coef_df_ml
  }
  
  if(!is.null(lty_by)) plot_type <- "line_type"
  
  if (!checkmate::test_directory_exists(out_dir)) {
    dir.create(out_dir, recursive=TRUE)
  }
  
  stopifnot(pdf_by %in% names(to_plot))
  
  to_plot <- to_plot %>%
    dplyr::rename(..p = !!p.value) %>%
    mutate(
      p_level = case_when(
        ..p > .05 ~ '1',
        ..p < .05 & ..p > .01 ~ '2',
        ..p < .01 & ..p > .001 ~ '3',
        ..p <.001 ~ '4'),
      p_level = ordered(p_level, levels = c('1', '2', '3', '4'), labels = c("NS","p < .05", "p < .01", "p < .001"))
    )
  
  if (!is.null(term_filter)) {
    to_plot <- to_plot %>% filter(grepl(term_filter, term))
  }
  
  make_plot <- function(data, title=NULL) {
    if (plot_type=="line") {
      g <- ggplot(data, aes_string(x=x, y=y, color=color, ymin=ymin, ymax=ymax)) +
        geom_line(size=1, position=position_dodge(width=0.4)) + 
        geom_pointrange(aes(size=p_level), position=position_dodge(width=0.4)) +
        #scale_color_brewer(palette="Dark2", labels=c("1" = "MT+, control", "2" = "Caudal post. parietal", "3" = "Rostral post. parietal", "4" = "Frontal premotor")) +
        scale_color_manual(values = colors) +
        geom_hline(yintercept = 0, size=1.5, alpha=0.6) +
        geom_vline(xintercept = 0, size=1.5, alpha=0.6) +
        scale_size_manual(values=c(0.5, 0.8, 1.1, 1.4)) + theme_bw(base_size=15)
    } else if (plot_type=="line_type") {
      g <- ggplot(data, aes_string(x=x, y=y, color=color, ymin=ymin, ymax=ymax, lty = lty_by)) +
        geom_line(size=1, position=position_dodge(width=0.4)) + 
        geom_pointrange(aes(size=p_level), position=position_dodge(width=0.4)) +
        scale_color_manual(values = colors) +
        # scale_color_brewer(palette="Dark2", labels=c("1" = "MT+, control", "2" = "Caudal post. parietal", "3" = "Rostral post. parietal", "4" = "Frontal premotor")) +
        geom_hline(yintercept = 0, size=1.5, alpha=0.6) +
        geom_vline(xintercept = 0, size=1.5, alpha=0.6) +
        scale_size_manual(values=c(0.5, 0.8, 1.1, 1.4)) + theme_bw(base_size=15)
    } else if (plot_type == "heat") {
      g <-  ggplot(data, aes_string(x=x, y=y, fill=color)) +
        geom_tile() +
        scale_fill_viridis_c() +
        coord_flip() +
        theme_bw(base_size=15)
    }
    
    if (isTRUE(flip)) {
      g <- g + coord_flip()
    }
    
    if (!is.null(title)) {
      g <- g + ggtitle(title)
    }
    
    if (!is.null(facet_by)) {
      g <- g + facet_wrap({{facet_by}})
    }
    
    if (!is.null(lty_by)) {
      g <- g + facet_wrap({{facet_by}}) #?
    }
    
    return(g)
  }
  
  to_plot_list <- split(to_plot, by=pdf_by)
  for (ff in seq_along(to_plot_list)) {
    this_df <- to_plot_list[[ff]]
    this_name <- make.names(names(to_plot_list[ff]))
    
    if (!is.null(panel_by)) {
      data_to_plot <- split(this_df, by=panel_by)
    } else {
      data_to_plot <- list(this_df) # single element list
    }
    
    glist <- lapply(seq_along(data_to_plot), function(x) { make_plot(data_to_plot[[x]], title=names(data_to_plot)[x]) })
    
    pdf(file.path(out_dir, paste0(this_name, ".pdf")), width=width, height=height)
    #plot(g)
    pobj <- patchwork::wrap_plots(glist) + 
      plot_annotation(title = this_name)
    
    plot(pobj)
    dev.off()
  }
  
}

