
Variablefeatureremove_mouse_genelst <- c("Trav", "Trbv", "Trdv", "Trgv", "Igkc", "Igkv", "Iglc", "Iglv",
                                         "Ighg", "Ighv", "Igha", "Ighm", "Jchain", "Rps", "Rpl", "mt-")

`%||%` <- function(lhs, rhs) {
  if (!is.null(x = lhs)) {
    return(lhs)
  } else {
    return(rhs)
  }
}

left_join2 <- function(x, y, by = NULL) {
  tmp <- nrow(x)
  output <- x %>%
    left_join(y, by = by)
  if (nrow(output) != tmp) {
    print("something wrong")
    warning("something wrong")
  }
  return(output)
}

rps_trv_ig_mt_mouse_exclusion <- function(input) {
  tmp <- input[!str_detect(rownames(input), "Rps") &
                 !str_detect(rownames(input), "Rpl") &
                 !str_detect(rownames(input), "mt-") &
                 !str_detect(rownames(input), "Trav") &
                 !str_detect(rownames(input), "Trbv") &
                 !str_detect(rownames(input), "Trdv") &
                 !str_detect(rownames(input), "Trgv") &
                 !str_detect(rownames(input), "Igkc") &
                 !str_detect(rownames(input), "Igkv") &
                 !str_detect(rownames(input), "Iglc") &
                 !str_detect(rownames(input), "Iglv") &
                 !str_detect(rownames(input), "Ighg") &
                 !str_detect(rownames(input), "Igha") &
                 !str_detect(rownames(input), "Ighm") &
                 !str_detect(rownames(input), "Jchain") &
                 !str_detect(rownames(input), "Ighv"), ]
  return(tmp)
}

VariableFeatures_genelst <- function(object, genelst) {
  for (this_gene in genelst) {
    VariableFeatures(object) <- VariableFeatures(object)[!str_detect(VariableFeatures(object), this_gene)]
  }
  return(object)
}

ggColorHue <- function(n, l = 65) {
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = l, c = 100)[1:n]
}

RedColorHue <- function(n) {
  hues <- seq(0, 0.9, length = n + 1)
  rgb(r = 1, g = hues, b = hues)[1:n]
}
BlueColorHue <- function(n) {
  hues <- seq(0, 0.9, length = n + 1)
  rgb(r = hues, g = hues, b = 1)[1:n]
}
GreenColorHue <- function(n) {
  hues <- seq(0, 0.9, length = n + 1)
  rgb(r = hues, g = 1, b = hues)[1:n]
}

DotPlot2 <- function (object,
                      assay = NULL,
                      features,
                      cols = c("lightgrey", "blue"),
                      col.min = -2.5,
                      col.max = 2.5,
                      dot.min = 0,
                      dot.scale = 6,
                      idents = NULL,
                      group.by = NULL,
                      split.by = NULL,
                      cluster.idents = FALSE,
                      scale = TRUE,
                      scale.by = "radius",
                      scale.min = NA,
                      scale.max = NA)
{
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  split.colors <- !is.null(x = split.by) && !any(cols %in%
                                                   rownames(x = brewer.pal.info))
  scale.func <- switch(EXPR = scale.by,
                       size = scale_size,
                       radius = scale_radius,
                       stop("'scale.by' must be either 'size' or 'radius'"))
  feature.groups <- NULL
  if (is.list(features) | any(!is.na(names(features)))) {
    feature.groups <- unlist(x = sapply(
      X = 1:length(features),
      FUN = function(x) {
        return(rep(x = names(x = features)[x], each = length(features[[x]])))
      }
    ))
    if (any(is.na(x = feature.groups))) {
      warning("Some feature groups are unnamed.",
              call. = FALSE,
              immediate. = TRUE)
    }
    features <- unlist(x = features)
    names(x = feature.groups) <- features
  }
  cells <- unlist(x = CellsByIdentities(object = object, idents = idents))
  data.features <- FetchData(object = object,
                             vars = features,
                             cells = cells)
  data.features$id <- if (is.null(x = group.by)) {
    Idents(object = object)[cells, drop = TRUE]
  }
  else {
    object[[group.by, drop = TRUE]][cells, drop = TRUE]
  }
  if (!is.factor(x = data.features$id)) {
    data.features$id <- factor(x = data.features$id)
  }
  id.levels <- levels(x = data.features$id)
  data.features$id <- as.vector(x = data.features$id)
  if (!is.null(x = split.by)) {
    splits <- object[[split.by, drop = TRUE]][cells, drop = TRUE]
    if (split.colors) {
      if (length(x = unique(x = splits)) > length(x = cols)) {
        stop("Not enough colors for the number of groups")
      }
      cols <- cols[1:length(x = unique(x = splits))]
      names(x = cols) <- unique(x = splits)
    }
    data.features$id <- paste(data.features$id, splits, sep = "_")
    unique.splits <- unique(x = splits)
    id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)),
                        "_",
                        rep(x = unique(x = splits), times = length(x = id.levels)))
  }
  data.plot <- lapply(
    X = unique(x = data.features$id),
    FUN = function(ident) {
      data.use <- data.features[data.features$id == ident, 1:(ncol(x = data.features) - 1), drop = FALSE]
      avg.exp <- apply(
        X = data.use,
        MARGIN = 2,
        FUN = function(x) {
          return(mean(x = expm1(x = x)))
        }
      )
      pct.exp <- apply(
        X = data.use,
        MARGIN = 2,
        FUN = PercentAbove,
        threshold = 0
      )
      return(list(avg.exp = avg.exp, pct.exp = pct.exp))
    }
  )
  names(x = data.plot) <- unique(x = data.features$id)
  if (cluster.idents) {
    mat <- do.call(what = rbind,
                   args = lapply(X = data.plot, FUN = unlist))
    mat <- scale(x = mat)
    id.levels <- id.levels[hclust(d = dist(x = mat))$order]
  }
  data.plot <- lapply(
    X = names(x = data.plot),
    FUN = function(x) {
      data.use <- as.data.frame(x = data.plot[[x]])
      data.use$features.plot <- rownames(x = data.use)
      data.use$id <- x
      return(data.use)
    }
  )
  data.plot <- do.call(what = "rbind", args = data.plot)
  if (!is.null(x = id.levels)) {
    data.plot$id <- factor(x = data.plot$id, levels = id.levels)
  }
  ngroup <- length(x = levels(x = data.plot$id))
  if (ngroup == 1) {
    scale <- FALSE
    warning(
      "Only one identity present, the expression values will be not scaled",
      call. = FALSE,
      immediate. = TRUE
    )
  }
  else if (ngroup < 5 & scale) {
    warning(
      "Scaling data with a low number of groups may produce misleading results",
      call. = FALSE,
      immediate. = TRUE
    )
  }
  avg.exp.scaled <- sapply(
    X = unique(x = data.plot$features.plot),
    FUN = function(x) {
      data.use <- data.plot[data.plot$features.plot ==
                              x, "avg.exp"]
      if (scale) {
        data.use <- scale(x = data.use)
        data.use <- MinMax(data = data.use,
                           min = col.min,
                           max = col.max)
      }
      else {
        data.use <- log1p(x = data.use)
        data.use <- MinMax(data = data.use,
                           min = col.min,
                           max = col.max)  
      }
      return(data.use)
    }
  )
  avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
  if (split.colors) {
    avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled, breaks = 20))
  }
  data.plot$avg.exp.scaled <- avg.exp.scaled
  data.plot$features.plot <- factor(x = data.plot$features.plot, levels = features)
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100
  if (split.colors) {
    splits.use <- vapply(
      X = as.character(x = data.plot$id),
      FUN = gsub,
      FUN.VALUE = character(length = 1L),
      pattern = paste0("^((", paste(
        sort(x = levels(x = object), decreasing = TRUE), collapse = "|"
      ), ")_)"),
      replacement = "",
      USE.NAMES = FALSE
    )
    data.plot$colors <- mapply(
      FUN = function(color, value) {
        return(colorRampPalette(colors = c("grey", color))(20)[value])
      },
      color = cols[splits.use],
      value = avg.exp.scaled
    )
  }
  color.by <- ifelse(test = split.colors, yes = "colors", no = "avg.exp.scaled")
  if (!is.na(x = scale.min)) {
    data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
  }
  if (!is.na(x = scale.max)) {
    data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
  }
  if (!is.null(x = feature.groups)) {
    data.plot$feature.groups <- factor(x = feature.groups[data.plot$features.plot], levels = unique(x = feature.groups))
  }
  plot <- ggplot(data = data.plot,
                 mapping = aes_string(x = "features.plot", y = "id")) + geom_point(mapping = aes_string(size = "pct.exp", color = color.by)) +
    scale.func(range = c(0, dot.scale),
               limits = c(scale.min, scale.max)) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    guides(size = guide_legend(title = "Percent Expressed")) +
    labs(x = "Features",
         y = ifelse(
           test = is.null(x = split.by),
           yes = "Identity",
           no = "Split Identity"
         )) + theme_cowplot()
  if (!is.null(x = feature.groups)) {
    plot <- plot + facet_grid(
      facets = ~ feature.groups,
      scales = "free_x",
      space = "free_x",
      switch = "y"
    ) + theme(panel.spacing = unit(x = 1, units = "lines"),
              strip.background = element_blank())
  }
  if (split.colors) {
    plot <- plot + scale_color_identity()
  }
  else if (length(x = cols) == 1) {
    plot <- plot + scale_color_distiller(palette = cols)
  }
  else {
    plot <- plot + scale_color_gradient(low = cols[1], high = cols[2])
  }
  if (!split.colors) {
    plot <- plot + guides(color = guide_colorbar(title = "Average Expression"))
  }
  return(plot)
}

##' GSEA plot that mimic the plot generated by broad institute's GSEA software
##'
##'
##' @title gseaplot2
##' @param x gseaResult object
##' @param geneSetID gene set ID
##' @param title plot title
##' @param color color of running enrichment score line
##' @param base_size base font size
##' @param rel_heights relative heights of subplots
##' @param subplots which subplots to be displayed
##' @param pvalue_table whether add pvalue table
##' @param ES_geom geom for plotting running enrichment score,
##' one of 'line' or 'dot'
##' @return plot
##' @export
##' @importFrom ggplot2 theme_classic
##' @importFrom ggplot2 element_line
##' @importFrom ggplot2 element_text
##' @importFrom ggplot2 element_blank
##' @importFrom ggplot2 element_rect
##' @importFrom ggplot2 scale_x_continuous
##' @importFrom ggplot2 scale_y_continuous
##' @importFrom ggplot2 scale_color_manual
##' @importFrom ggplot2 theme_void
##' @importFrom ggplot2 geom_rect
##' @importFrom ggplot2 margin
##' @importFrom ggplot2 annotation_custom
##' @importFrom stats quantile
##' @importFrom RColorBrewer brewer.pal
##' @author Guangchuang Yu
gseaplot2_2 <- function(x,
                        geneSetID,
                        title = "",
                        color = "green",
                        base_size = 11,
                        rel_heights = c(1.5, .5, 1),
                        subplots = 1:3,
                        pvalue_table = FALSE,
                        ES_geom = "line",
                        ylimit = NULL,
                        pval_pos_x = 200,
                        pval_pos_y = 0) {
  ES_geom <- match.arg(ES_geom, c("line", "dot"))
  
  geneList <- position <- NULL ## to satisfy codetool
  
  gsInfo <- function(object, geneSetID) {
    geneList <- object@geneList
    
    if (is.numeric(geneSetID))
      geneSetID <- object@result[geneSetID, "ID"]
    
    geneSet <- object@geneSets[[geneSetID]]
    exponent <- object@params[["exponent"]]
    df <- gseaScores(geneList, geneSet, exponent, fortify = TRUE)
    df$ymin <- 0
    df$ymax <- 0
    pos <- df$position == 1
    h <- diff(range(df$runningScore)) / 20
    df$ymin[pos] <- -h
    df$ymax[pos] <- h
    df$geneList <- geneList
    
    df$Description <- object@result[geneSetID, "Description"]
    return(df)
  }
  
  gseaScores <- getFromNamespace("gseaScores", "DOSE")
  
  if (length(geneSetID) == 1) {
    gsdata <- gsInfo(x, geneSetID)
  } else {
    gsdata <- do.call(rbind, lapply(geneSetID, gsInfo, object = x))
  }
  
  p <- ggplot(gsdata, aes_(x = ~ x)) + xlab(NULL) +
    theme_classic(base_size) +
    scale_x_continuous(expand = c(0, 0))
  
  if (ES_geom == "line") {
    es_layer <- geom_line(aes_(y = ~ runningScore, color = ~ Description), size =
                            1)
  } else {
    es_layer <- geom_point(
      aes_(y = ~ runningScore, color = ~ Description),
      size = 1,
      data = subset(gsdata, position == 1)
    )
  }
  
  p.res <- p + es_layer +
    theme(
      legend.position = c(.8, .8),
      legend.title = element_blank(),
      legend.background = element_rect(fill = "transparent")
    )
  
  if (!is.null(ylimit)) {
    p.res <- p.res + scale_y_continuous(limits = ylimit)
  }
  
  p.res <- p.res +
    geom_hline(yintercept = 0, linetype = "dashed") +
    ylab("Running Enrichment Score") +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
      plot.margin = margin(
        t = .2,
        r = .2,
        b = 0,
        l = .2,
        unit = "cm"
      ),
      panel.border = element_rect(fill = NA, size = 1)
    )
  
  i <- 0
  for (term in unique(gsdata$Description)) {
    idx <- which(gsdata$ymin != 0 & gsdata$Description == term)
    gsdata[idx, "ymin"] <- i
    gsdata[idx, "ymax"] <- i + 1
    i <- i + 1
  }
  
  p2 <- ggplot(gsdata, aes_(x = ~ x)) +
    geom_linerange(aes_(
      ymin =  ~ ymin,
      ymax =  ~ ymax,
      color =  ~ Description
    )) +
    xlab(NULL) + ylab(NULL) + theme_classic(base_size) +
    theme(
      legend.position = "none",
      plot.margin = margin(t = -.1, b = 0, unit = "cm"),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.line.x = element_blank(),
      panel.border = element_rect(fill = NA, size = 1)
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0))
  
  if (length(geneSetID) == 1) {
    v <- seq(1, sum(gsdata$position), length.out = 9)
    inv <- findInterval(rev(cumsum(gsdata$position)), v)
    if (min(inv) == 0)
      inv <- inv + 1
    
    col <- c(rev(brewer.pal(5, "Blues")), brewer.pal(5, "Reds"))
    
    ymin <- min(p2$data$ymin)
    yy <- max(p2$data$ymax - p2$data$ymin) * .3
    xmin <- which(!duplicated(inv))
    xmax <- xmin + as.numeric(table(inv)[as.character(unique(inv))])
    d <- data.frame(
      ymin = ymin,
      ymax = yy,
      xmin = xmin,
      xmax = xmax,
      col = col[unique(inv)]
    )
    p2 <- p2 + geom_rect(
      aes_(
        xmin =  ~ xmin,
        xmax =  ~ xmax,
        ymin =  ~ ymin,
        ymax =  ~ ymax,
        fill =  ~ I(col)
      ),
      data = d,
      # alpha=.9,
      inherit.aes = FALSE
    )
  }
  
  
  df2 <- p$data
  df2$y <- p$data$geneList[df2$x]
  p.pos <- p + geom_segment(data = df2,
                            aes_(
                              x =  ~ x,
                              xend =  ~ x,
                              y =  ~ y,
                              yend = 0
                            ),
                            color = "grey")
  p.pos <- p.pos + ylab("Ranked List Metric") +
    xlab("Rank in ordered gene list") +
    theme(
      plot.margin = margin(
        t = -.1,
        r = .2,
        b = .2,
        l = .2,
        unit = "cm"
      ),
      panel.border = element_rect(fill = NA, size = 1)
    )
  
  if (!is.null(title) && !is.na(title) && title != "")
    p.res <- p.res + ggtitle(title)
  
  if (length(color) == length(geneSetID)) {
    p.res <- p.res + scale_color_manual(values = color)
    if (length(color) == 1) {
      p.res <- p.res + theme(legend.position = "none")
      p2 <- p2 + scale_color_manual(values = "black")
    } else {
      p2 <- p2 + scale_color_manual(values = color)
    }
  }
  
  if (pvalue_table) {
    pd <- x[geneSetID, c("Description", "pvalue", "qvalue")]
    # pd <- pd[order(pd[,1], decreasing=FALSE),]
    
    if (length(geneSetID) == 1) {
      rownames(pd) <- geneSetID
    } else {
      rownames(pd) <- pd$Description
    }
    
    pd <- pd[, -1]
    # pd <- round(pd, 4)
    for (i in seq_len(ncol(pd))) {
      pd[, i] <- format(pd[, i], digits = 4)
    }
    
    tableGrob2 <- function(d, p = NULL) {
      # has_package("gridExtra")
      d <- d[order(rownames(d)), ]
      tp <- gridExtra::tableGrob(d)
      if (is.null(p)) {
        return(tp)
      }
      
      # Fix bug: The 'group' order of lines and dots/path is different
      p_data <- ggplot_build(p)$data[[1]]
      # pcol <- unique(ggplot_build(p)$data[[1]][["colour"]])
      p_data <- p_data[order(p_data[["group"]]), ]
      pcol <- unique(p_data[["colour"]])
      ## This is fine too
      ## pcol <- unique(p_data[["colour"]])[unique(p_data[["group"]])]
      j <- which(tp$layout$name == "rowhead-fg")
      
      for (i in seq_along(pcol)) {
        tp$grobs[j][[i + 1]][["gp"]] <- gpar(col = pcol[i])
      }
      return(tp)
    }
    
    if (length(geneSetID) == 1) {
      ES <- format(x[geneSetID, "enrichmentScore"], digits = 4)
      NES <- format(x[geneSetID, "NES"], digits = 4)
      tp <- str_glue("p-value: {pd$pvalue[1]} \n q-value: {pd$qvalue[1]} \n ES: {ES} \n NES: {NES}")
      p.res <- p.res + theme(legend.position = "none") +
        annotate("text",
                 label = tp,
                 x = pval_pos_x,
                 y = pval_pos_y)
    } else {
      tp <- tableGrob2(pd, p.res)
      p.res <- p.res + theme(legend.position = "none") +
        annotation_custom(
          tp,
          xmin = quantile(p.res$data$x, .5),
          xmax = quantile(p.res$data$x, .95),
          ymin = quantile(p.res$data$runningScore, .75),
          ymax = quantile(p.res$data$runningScore, .9)
        )
    }
    
  }
  
  plotlist <- list(p.res, p2, p.pos)[subplots]
  n <- length(plotlist)
  plotlist[[n]] <- plotlist[[n]] +
    theme(
      axis.line.x = element_line(),
      axis.ticks.x = element_line(),
      axis.text.x = element_text()
    )
  
  if (length(subplots) == 1)
    return(plotlist[[1]] + theme(plot.margin = margin(
      t = .2,
      r = .2,
      b = .2,
      l = .2,
      unit = "cm"
    )))
  
  
  if (length(rel_heights) > length(subplots))
    rel_heights <- rel_heights[subplots]
  
  aplot::gglist(gglist = plotlist,
                ncol = 1,
                heights = rel_heights)
}
