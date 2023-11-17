
phemirror<-function (top, bottom, phegroup, tline, bline, chroms = c(1:22, 
    "X", "Y"), log10 = TRUE, yaxis, opacity = 1, annotate_snp, 
    annotate_p, highlight_snp, highlight_p, highlighter = "red", 
    toptitle = NULL, bottomtitle = NULL, chrcolor1 = "#AAAAAA", 
    chrcolor2 = "#4D4D4D", groupcolors, freey = FALSE, background = "variegated", 
    chrblocks = TRUE, file = "phemirror", type = "png", hgtratio = 0.5, 
    hgt = 7, wi = 12, res = 800) 
{
    topn <- names(top)
    bottomn <- names(bottom)
    top$Location <- "Top"
    bottom$Location <- "Bottom"
    if (!identical(topn, bottomn)) {
        stop("Please ensure both inputs have the same metadata columns.")
    }
    d <- as.data.frame(rbind(top, bottom))
    d$POS <- as.numeric(as.character(d$POS))
    d$CHR <- droplevels(factor(d$CHR, levels = as.character(chroms)))
    d <- d[d$CHR %in% chroms, ]
    if (!missing(phegroup)) {
        print("Only phenotypes with grouping information will be plotted")
        d_phe <- merge(phegroup, d, by = "PHE")
        names(d_phe)[names(d_phe) == "Group"] <- "Color"
    }
    else {
        d_phe <- d
        names(d_phe)[names(d_phe) == "PHE"] <- "Color"
    }
    d_order <- d_phe[order(d_phe$CHR, d_phe$POS), ]
    d_order$pos_index <- seq.int(nrow(d_order))
    d_order_sub <- d_order[, c("SNP", "CHR", "POS", "pvalue", 
        "pos_index")]
    maxRows <- by(d_order_sub, d_order_sub$CHR, function(x) x[which.max(x$pos_index), 
        ])
    minRows <- by(d_order_sub, d_order_sub$CHR, function(x) x[which.min(x$pos_index), 
        ])
    milimits <- do.call(rbind, minRows)
    malimits <- do.call(rbind, maxRows)
    lims <- merge(milimits, malimits, by = "CHR")
    names(lims) <- c("Color", "snpx", "px", "posx", "posmin", 
        "snpy", "py", "posy", "posmax")
    lims$av <- (lims$posmin + lims$posmax)/2
    lims <- lims[order(lims$Color), ]
    lims$shademap <- rep(c("shade_ffffff", "shade_ebebeb"), length.out = nrow(lims), 
        each = 1)
    nchrcolors <- nlevels(factor(lims$Color))
    base_color <- c(rep(x = c(chrcolor1, chrcolor2), length.out = nchrcolors, 
        each = 1), "#FFFFFF", "#EBEBEB")
    names(base_color) <- c(levels(factor(lims$Color)), "shade_ffffff", 
        "shade_ebebeb")
    if (!missing(groupcolors)) {
        topcols <- groupcolors
        bottomcols <- groupcolors
    }
    else {
        ngroupcolors <- nlevels(factor(d_order$Color[d_order$Location == 
            "Top"]))
        if (ngroupcolors > 15) {
            topcols <- hudson:::Turbo(out.colors = ngroupcolors)
        }
        else {
            pal <- c("#009292", "#920000", "#490092", "#db6d00", 
                "#24ff24", "#ffff6d", "#000000", "#006ddb", "#004949", 
                "#924900", "#ff6db6", "#6db6ff", "#b66dff", "#ffb6db", 
                "#b6dbff")
            topcols <- pal[1:ngroupcolors]
        }
        names(topcols) <- levels(factor(d_order$Color[d_order$Location == 
            "Top"]))
        ngroupcolors <- nlevels(factor(d_order$Color[d_order$Location == 
            "Bottom"]))
        if (ngroupcolors > 15) {
            bottomcols <- hudson:::Turbo(ngroupcolors)
        }
        else {
            pal <- c("#009292", "#920000", "#490092", "#db6d00", 
                "#24ff24", "#ffff6d", "#000000", "#006ddb", "#004949", 
                "#924900", "#ff6db6", "#6db6ff", "#b66dff", "#ffb6db", 
                "#b6dbff")
            bottomcols <- pal[1:ngroupcolors]
        }
        names(bottomcols) <- levels(factor(d_order$Color[d_order$Location == 
            "Bottom"]))
    }
    if (log10 == TRUE) {
        d_order$pval <- -log10(d_order$pvalue)
        yaxislab1 <- expression(paste("-log"[10], "(p-value)", 
            sep = ""))
        yaxislab2 <- expression(paste("-log"[10], "(p-value)", 
            sep = ""))
        if (!missing(tline)) {
            tredline <- -log10(tline)
        }
        if (!missing(bline)) {
            bredline <- -log10(bline)
        }
    }
    else {
        d_order$pval <- d_order$pvalue
        yaxislab1 <- yaxis[1]
        yaxislab2 <- yaxis[2]
        if (!missing(tline)) {
            tredline <- tline
        }
        if (!missing(bline)) {
            bredline <- bline
        }
    }
    yaxismax1 <- ifelse(freey == FALSE, max(d_order$pval[which(d_order$pval < 
        Inf)]), max(d_order$pval[which(d_order$pval < Inf) & 
        d_order$Location == "Top"]))
    yaxismax2 <- ifelse(freey == FALSE, max(d_order$pval[which(d_order$pval < 
        Inf)]), max(d_order$pval[which(d_order$pval < Inf) & 
        d_order$Location == "Bottom"]))
    yaxismin1 <- ifelse(freey == FALSE, 0, min(d_order$pval[d_order$Location == 
        "Top"]))
    yaxismin2 <- ifelse(freey == FALSE, 0, min(d_order$pval[d_order$Location == 
        "Bottom"]))
    backpanel1 <- ifelse(background == "#c0bcbc", "NULL", "geom_rect(data = lims, aes(xmin = posmin-.5, xmax = posmax+.5, ymin = yaxismin1, ymax = Inf, fill=factor(shademap)), alpha = 0.5)")
    backpanel2 <- ifelse(background == "#c0bcbc", "NULL", "geom_rect(data = lims, aes(xmin = posmin-.5, xmax = posmax+.5, ymin = yaxismin2, ymax = Inf, fill=factor(shademap)), alpha = 0.5)")
    p1 <- ggplot() + eval(parse(text = backpanel1))
    if ("Shape" %in% topn) {
        p1 <- p1 + geom_point(data = d_order[d_order$Location == 
            "Top", ], aes(x = pos_index, y = pval, color = factor(Color), 
            shape = factor(Shape)), alpha = opacity,size = 0.8)
    }
    else {
        p1 <- p1 + geom_point(data = d_order[d_order$Location == 
            "Top", ], aes(x = pos_index, y = pval, color = factor(Color)), 
            alpha = opacity,size = 0.8)
    }
    p1 <- p1 + scale_x_continuous(breaks = lims$av, labels = lims$Color, 
        expand = c(0, 0))
    if (chrblocks == TRUE) {
        p1 <- p1 + geom_rect(data = lims, aes(xmin = posmin - 
            0.5, xmax = posmax + 0.5, ymin = -Inf, ymax = min(d_order$pval), 
            fill = as.factor(Color)), alpha = 1)
    }
    p1 <- p1 + scale_colour_manual(name = "Color", values = topcols) + 
        scale_fill_manual(values = base_color)
    p1 <- p1 + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), 
        axis.title.x = element_blank(), legend.position = "top", 
        legend.title = element_blank())
    p2 <- ggplot() + eval(parse(text = backpanel2))
    if ("Shape" %in% bottomn) {
        p2 <- p2 + geom_point(data = d_order[d_order$Location == 
            "Bottom", ], aes(x = pos_index, y = pval, color = factor(Color), 
            shape = factor(Shape)), alpha = opacity,size = 1.1)
    }
    else {
        p2 <- p2 + geom_point(data = d_order[d_order$Location == 
            "Bottom", ], aes(x = pos_index, y = pval, color = factor(Color)), 
            alpha = opacity,size = 1.1)
    }
    p2 <- p2 + scale_x_continuous(breaks = lims$av, labels = lims$Color, 
        expand = c(0, 0))
    if (chrblocks == TRUE) {
        p2 <- p2 + geom_rect(data = lims, aes(xmin = posmin - 
            0.5, xmax = posmax + 0.5, ymin = -Inf, ymax = min(d_order$pval), 
            fill = as.factor(Color)), alpha = 1)
    }
    p2 <- p2 + scale_colour_manual(name = "Color", values = bottomcols) + 
        scale_fill_manual(values = base_color)
    p2 <- p2 + theme(axis.text.x = element_text(angle = 90), 
        panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), 
        axis.title.x = element_blank(), legend.position = "bottom", 
        legend.title = element_blank())
    if (!missing(highlight_snp)) {
        if ("Shape" %in% topn) {
            p1 <- p1 + geom_point(data = d_order[d_order$SNP %in% 
                highlight_snp & d_order$Location == "Top", ], 
                aes(x = pos_index, y = pval, shape = Shape), 
                colour = highlighter)
            p1 <- p1 + guides(shape = guide_legend(override.aes = list(colour = "black")))
        }
        else {
            p1 <- p1 + geom_point(data = d_order[d_order$SNP %in% 
                highlight_snp & d_order$Location == "Top", ], 
                aes(x = pos_index, y = pval), colour = highlighter)
        }
        if ("Shape" %in% bottomn) {
            p2 <- p2 + geom_point(data = d_order[d_order$SNP %in% 
                highlight_snp & d_order$Location == "Bottom", 
                ], aes(x = pos_index, y = pval, shape = Shape), 
                colour = highlighter)
            p2 <- p2 + guides(shape = guide_legend(override.aes = list(colour = "black")))
        }
        else {
            p2 <- p2 + geom_point(data = d_order[d_order$SNP %in% 
                highlight_snp & d_order$Location == "Bottom", 
                ], aes(x = pos_index, y = pval), colour = highlighter)
        }
    }
    if (!missing(highlight_p)) {
        if ("Shape" %in% topn) {
            p1 <- p1 + geom_point(data = d_order[d_order$pvalue < 
                highlight_p[1] & d_order$Location == "Top", ], 
                aes(x = pos_index, y = pval, shape = Shape), 
                colour = highlighter)
            p1 <- p1 + guides(shape = guide_legend(override.aes = list(colour = "black")))
        }
        else {
            p1 <- p1 + geom_point(data = d_order[d_order$pvalue < 
                highlight_p[1] & d_order$Location == "Top", ], 
                aes(x = pos_index, y = pval), colour = highlighter)
        }
        if ("Shape" %in% bottomn) {
            p2 <- p2 + geom_point(data = d_order[d_order$pvalue < 
                highlight_p[2] & d_order$Location == "Bottom", 
                ], aes(x = pos_index, y = pval, shape = Shape), 
                colour = highlighter)
            p2 <- p2 + guides(shape = guide_legend(override.aes = list(colour = "black")))
        }
        else {
            p2 <- p2 + geom_point(data = d_order[d_order$pvalue < 
                highlight_p[2] & d_order$Location == "Bottom", 
                ], aes(x = pos_index, y = pval), colour = highlighter)
        }
    }
    if (!missing(tline)) {
        for (i in 1:length(tline)) {
            p1 <- p1 + geom_hline(yintercept = tredline[i], colour = "red")
        }
    }
    if (!missing(bline)) {
        for (i in 1:length(bline)) {
            p2 <- p2 + geom_hline(yintercept = bredline[i], colour = "red")
        }
    }
    if (!missing(annotate_p)) {
        if (!requireNamespace(c("ggrepel"), quietly = TRUE) == 
            TRUE) {
            print("Consider installing 'ggrepel' for improved text annotation")
            p1 <- p1 + geom_text(data = d_order[d_order$pvalue < 
                annotate_p[1] & !is.na(d_order$SNP) & d_order$Location == "Top", ], 
                aes(pos_index, pval, label = SNP),size = 1.4,fontface ="bold", max.overlaps = 100,color="black")
            #p1 + geom_text(aes(color=factor(SNP2)))
            p2 <- p2 + geom_text(data = d_order[d_order$pvalue < 
                annotate_p[2] & !is.na(d_order$SNP) & d_order$Location == "Bottom",],
                 aes(pos_index, pval, label = SNP),size = 1.4,fontface ="bold",max.overlaps = 100,color="black")
        }
        else {
            p1 <- p1 + ggrepel::geom_text_repel(data = d_order[d_order$pvalue < 
                annotate_p[1] & !is.na(d_order$SNP) & d_order$Location == "Top", ], 
                aes(pos_index, pval, label = SNP),size = 1.4,fontface ="bold",max.overlaps = 100)


            p2 <- p2 + ggrepel::geom_text_repel(data = d_order[d_order$pvalue < 
                annotate_p[2] & !is.na(d_order$SNP) & d_order$Location == "Bottom", 
                ], aes(pos_index, pval, label = SNP),size = 1.4,fontface ="bold",max.overlaps = 100)
        }
    }
    if (!missing(annotate_snp)) {
        if (!requireNamespace(c("ggrepel"), quietly = TRUE) == 
            TRUE) {
            print("Consider installing 'ggrepel' for improved text annotation")
            p1 <- p1 + geom_text(data = d_order[d_order$SNP %in% 
                annotate_snp & d_order$Location == "Top", ], 
                aes(pos_index, pval, label = SNP),size = 1.5,fontface ="bold")
            p2 <- p2 + geom_text(data = d_order[d_order$SNP %in% 
                annotate_snp & d_order$Location == "Bottom", 
                ], aes(pos_index, pval, label = SNP),size = 1.5,fontface ="bold")
        }
        else {
            p1 <- p1 + ggrepel::geom_text_repel(data = d_order[d_order$SNP %in% 
                annotate_snp & d_order$Location == "Top", ], 
                aes(pos_index, pval, label = SNP),size = 1.5,fontface ="bold")
            p2 <- p2 + ggrepel::geom_text_repel(data = d_order[d_order$SNP %in% 
                annotate_snp & d_order$Location == "Bottom", 
                ], aes(pos_index, pval, label = SNP),size = 1.5,fontface ="bold")
        }
    }
    p1 <- p1 + ylab(yaxislab1)
    p2 <- p2 + ylab(yaxislab2)
    if (chrblocks == TRUE) {
        if (freey == TRUE) {
            print("Sorry, drawing chrblocks with freey=TRUE is currently unsupported and will be ignored.")
        }
        else {
            p1 <- p1 + theme(axis.text.x = element_text(vjust = 1), 
                axis.ticks.x = element_blank()) + ylim(c(yaxismin1, 
                yaxismax1))
            p2 <- p2 + scale_y_reverse(limits = c(yaxismax2, 
                yaxismin2)) + theme(axis.text.x = element_blank(), 
                axis.ticks.x = element_blank())
        }
    }
    else {
        p1 <- p1 + theme(axis.text.x = element_text(vjust = 1), 
            axis.ticks.x = element_blank()) + scale_y_continuous(limits = c(yaxismin1, 
            yaxismax1), expand = expansion(mult = c(0, 0.1)))
        p2 <- p2 + scale_y_reverse(limits = c(yaxismax2, yaxismin2), 
            expand = expansion(mult = c(0.1, 0))) + theme(axis.text.x = element_blank(), 
            axis.ticks.x = element_blank())
    }
    if (background == "white") {
        p1 <- p1 + theme(panel.background = element_rect(fill = "white"))
        p2 <- p2 + theme(panel.background = element_rect(fill = "white"))
    }
    p1 <- p1 + guides(fill = "none")
    p2 <- p2 + guides(fill = "none")



    # Add horizontal lines at y = 5
    #p1 <- p1 + geom_hline(yintercept = 5, color = "#848586", alpha = 0.3)  # Customize the color as needed
    #p2 <- p2 + geom_hline(yintercept = 5, color = "#848586", alpha = 0.3)  # Customize the color as needed

    p1 <- p1 + ggrepel::geom_text_repel(data = d_order[!is.na(d_order$SNP2) & d_order$Location == "Top", ], 
        aes(pos_index, pval, label = SNP2),size = 1.4,fontface ="bold",color="red",,max.overlaps = 100)
    p2 <- p2 + ggrepel::geom_text_repel(data = d_order[!is.na(d_order$SNP2) & d_order$Location == "Bottom", 
        ], aes(pos_index, pval, label = SNP2),size = 1.4,fontface ="bold",color="red",,max.overlaps = 100)


tryCatch({
   p1 <- p1 + ggrepel::geom_text_repel(data = d_order[!is.na(d_order$SNP3) & d_order$Location == "Top", ], 
       aes(pos_index, pval, label = SNP3),size = 1.4,fontface ="bold",color="#251bbb",,max.overlaps = 100)
   }, error = function(e) {
   # Error handling code
   cat("An error occurred:", conditionMessage(e), "\n")
   # Optionally, you can add code here to handle or log the error
   })

tryCatch({
   p2 <- p2 + ggrepel::geom_text_repel(data = d_order[!is.na(d_order$SNP3) & d_order$Location == "Bottom", 
       ], aes(pos_index, pval, label = SNP3),size = 1.4,fontface ="bold",color="#251bbb",,max.overlaps = 100)
   }, error = function(e) {
   # Error handling code
   cat("An error occurred:", conditionMessage(e), "\n")
   # Optionally, you can add code here to handle or log the error
   })


head(d_order)
    print(paste0("Saving plot to ", file, ".", type))
    p <- grid.arrange(arrangeGrob(p1, top = toptitle), arrangeGrob(p2, 
        bottom = bottomtitle), padding = 0, heights = c(hgtratio, 
        1 - hgtratio))
    ggsave(p, filename = paste0(file, ".", type), dpi = res, 
        units = "in", height = hgt, width = wi)
    return(p)
}
