# Generate KM-plots with p-values and median of time of strata.
plotKM.Treatment <- function(fit, ylim, palette = 'jco', hr = NULL){
    
    # Clean-up strata.
    names(fit$strata) <-  base::gsub('.*=', '', names(fit$strata))
    
    # Generate survival plot.
    x <- survminer::ggsurvplot(
        fit = fit,
        pval = F,
        size = .825,
        break.time.by = 500,
        break.y.by = .2, 
        palette = palette,
        risk.table = T,
        tables.height = .25,
        xlab = 'Treatment duration (in days)',
        axes.offset = F,
        ylim = c(0, 1.05),
        xlim = c(0, ylim), strata.labels = 'c',
        risk.table.col = 'strata', censor.shape = '+',
        fontsize = 3, 
        conf.int = F,
        surv.median.line = 'hv',
        risk.table.title = 'No. at risk',
        ggtheme = ggplot2::theme(
            legend.position = 'none',
            axis.title.x = ggtext::element_textbox_simple(width = NULL, halign = .5, family='Roboto', face = 'bold'),
            axis.title.y = ggtext::element_textbox_simple(size = 8, orientation = 'left-rotated',family='Roboto', face = 'bold', width = NULL, halign = .5),
            panel.grid.major.x = ggplot2::element_line(colour = 'grey90', linetype = 'dotted'),
            panel.grid.major.y = ggplot2::element_line(colour = '#E5E5E5', linetype = 'dotted'),
            panel.grid.minor.y = ggplot2::element_blank(),
            panel.background = ggplot2::element_rect(fill = NA, colour = 'black'),
            text = ggplot2::element_text(size = 8, family='Roboto', face = 'bold'),
            legend.text = ggtext::element_markdown(family='Roboto', face = 'bold')
        ),
        font.family = 'Roboto'
    )
    
    # Add the log-rank p-value.
    p.logrank <- survminer::surv_pvalue(fit = fit, method = 'log-rank', test.for.trend = F)
    x$plot <- x$plot + ggplot2::annotate('text', x = 1500, y = .95, size = 2.5, label = paste0('log-rank: ', p.logrank$pval.txt), family='Roboto')
    
    # Add the good vs. bad responder q-value.
    p.goodvsbad <- survminer::pairwise_survdiff(fit$call$formula, data = fit$call$data, p.adjust.method = 'none')
    if(all(dim(p.goodvsbad$p.value) == c(2,2))){
        p.goodvsbad <- p.goodvsbad$p.value[2,2]
        p.goodvsbad <- round(p.goodvsbad, 5)
        x$plot <- x$plot + ggplot2::annotate('text', x = 1500, y = .8, size = 2.5, label = paste0('Poor vs. Good Responders: log-rank q = ', p.goodvsbad), family='Roboto')
    }
    
    # Add the median + 95% CI treatment duration.
    medians <- survminer::surv_median(fit = fit) %>% dplyr::mutate(label = sprintf('%s: %s (%s-%s)', strata, median, lower, upper))
    x$plot <- x$plot + ggplot2::annotate('text', x = 1500, y = .7, size = 2.5, label = paste(medians$label, collapse = '\n'), family='Roboto')
    
    
    # Add HR (if two groups)
    if(!is.null(hr)){
        
        HR.CI <- round(summary(hr)$conf.int, 2)
        HR.p <- round(summary(hr)$waldtest[[3]], 2)
        HR.CI <- sprintf('HR (95%% CI): %s (%s - %s)', HR.CI[[1]], HR.CI[[3]], HR.CI[[4]])
        x$plot <- x$plot + ggplot2::annotate('text', x = 20, y = .95, label = HR.CI, size = 2.5, fontface = 'bold')
        
    }
    
    # Remove legends.
    x$plot <- x$plot + theme_Job + ggplot2::theme(legend.position = 'none')
    x$table <- x$table + theme_Job + ggplot2::theme(legend.position = 'none') + ggplot2::ylab(NULL)
    
    # Add perc. scaling.
    x$plot <- x$plot + ggplot2::scale_y_continuous(labels = scales::percent_format(), limits = c(0,1), expand = c(0,0.01))
    x$plot <- x$plot + ggplot2::ylab('Cum. Ongoing Treatment') + ggplot2::xlab(NULL)
    
    return(x)
}

plotKM.OS <- function(fit, ylim, palette = 'jco', hr = NULL){
    x <- plotKM.Treatment(fit, ylim, palette, hr)
    x$plot <- x$plot + ggplot2::ylab('Cum. Survival')
    x$table <- x$table + ggplot2::xlab('Overall survival')
    
    return(x)
    
}
