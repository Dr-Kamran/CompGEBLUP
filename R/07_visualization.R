#' Visualization Functions for CompGEBLUP
#'
#' Publication-quality plots for genomic prediction results
#'
#' @name visualization
#' @keywords internal
NULL

#' Plot GEBV Distribution
#'
#' @param model GBLUPModel object or data frame with GEBV column
#' @param type Type of plot: "histogram", "density", "boxplot"
#' @param by_env Logical, plot separately by environment
#' @param colors Colors for environments/groups
#' @param notch Logical, draw notched boxplots (default TRUE)
#' @param show_mean Logical, show mean as a point on boxplots
#' @param main Plot title
#' @param ... Additional arguments passed to plotting functions
#' @return Invisible NULL
#' @export
plot_gebv <- function(model, 
                      type = c("histogram", "density", "boxplot"),
                      by_env = FALSE,
                      colors = NULL,
                      notch = TRUE,
                      show_mean = TRUE,
                      main = NULL,
                      ...) {
  
  type <- match.arg(type)
  
  # Extract GEBV data
  if (inherits(model, "GBLUPModel")) {
    gebv <- model@gebv
  } else if (is.data.frame(model) && "GEBV" %in% names(model)) {
    gebv <- model
  } else {
    stop("model must be a GBLUPModel object or data frame with GEBV column")
  }
  
  if (is.null(main)) {
    main <- paste("GEBV Distribution")
  }
  
  # Default colors - publication quality palette
  if (is.null(colors)) {
    colors <- c("#6BAED6", "#74C476", "#FD8D3C", "#FDDC00", 
                "#9E9AC8", "#FC9272", "#A1D99B", "#9ECAE1")
  }
  
  if (type == "histogram") {
    if (by_env && "ENV" %in% names(gebv)) {
      envs <- unique(gebv$ENV)
      n_env <- length(envs)
      par(mfrow = c(ceiling(n_env/2), min(2, n_env)), mar = c(4, 4, 3, 1))
      
      for (i in seq_along(envs)) {
        env_data <- gebv$GEBV[gebv$ENV == envs[i]]
        hist(env_data,
             main = paste("Environment:", envs[i]),
             xlab = "GEBV",
             col = colors[i],
             border = "white",
             las = 1,
             ...)
        abline(v = mean(env_data, na.rm = TRUE), col = "red", lwd = 2, lty = 2)
      }
      par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)
    } else {
      hist(gebv$GEBV,
           main = main,
           xlab = "GEBV",
           col = colors[1],
           border = "white",
           las = 1,
           ...)
      abline(v = mean(gebv$GEBV, na.rm = TRUE), col = "red", lwd = 2, lty = 2)
    }
    
  } else if (type == "density") {
    if (by_env && "ENV" %in% names(gebv)) {
      envs <- unique(gebv$ENV)
      all_gebv <- gebv$GEBV[!is.na(gebv$GEBV)]
      first_env_data <- gebv$GEBV[gebv$ENV == envs[1] & !is.na(gebv$GEBV)]
      
      plot(density(first_env_data),
           main = main, xlab = "GEBV", ylab = "Density",
           col = colors[1], lwd = 3, xlim = range(all_gebv), las = 1, ...)
      
      for (i in 2:length(envs)) {
        env_data <- gebv$GEBV[gebv$ENV == envs[i] & !is.na(gebv$GEBV)]
        if (length(env_data) >= 2) {
          lines(density(env_data), col = colors[i], lwd = 3)
        }
      }
      legend("topright", legend = envs, col = colors[seq_along(envs)], lwd = 3, bty = "n")
    } else {
      gebv_valid <- gebv$GEBV[!is.na(gebv$GEBV)]
      plot(density(gebv_valid), main = main, xlab = "GEBV", ylab = "Density",
           col = colors[1], lwd = 3, las = 1, ...)
      polygon(density(gebv_valid), col = rgb(0.4, 0.7, 0.9, 0.3), border = colors[1])
    }
    
  } else if (type == "boxplot") {
    if (by_env && "ENV" %in% names(gebv)) {
      envs <- unique(gebv$ENV)
      n_env <- length(envs)
      
      boxplot(GEBV ~ ENV, data = gebv, main = main, xlab = "", ylab = "GEBV",
              col = colors[seq_len(n_env)], notch = notch, las = 1,
              outline = TRUE, outpch = 1, outcex = 0.7, ...)
      
      if (show_mean) {
        means <- tapply(gebv$GEBV, gebv$ENV, mean, na.rm = TRUE)
        points(seq_along(means), means, pch = 16, col = "red", cex = 1.2)
      }
    } else {
      boxplot(gebv$GEBV, main = main, ylab = "GEBV", col = colors[1],
              notch = notch, las = 1, ...)
      if (show_mean) {
        points(1, mean(gebv$GEBV, na.rm = TRUE), pch = 16, col = "red", cex = 1.2)
      }
    }
  }
  
  invisible(NULL)
}


#' Plot Phenotype Boxplots (Publication Quality)
#'
#' Creates notched boxplots with mean points like in wheat phenology studies.
#'
#' @param pheno Data frame with phenotype columns or GBLUPData object
#' @param traits Character vector of trait column names to plot
#' @param colors Named vector of colors for each trait
#' @param notch Logical, draw notched boxplots (default TRUE)
#' @param show_mean Logical, show mean as a red point
#' @param main Plot title (auto-generated if NULL)
#' @param ylab Y-axis label
#' @param ... Additional arguments passed to boxplot
#'
#' @return Invisible NULL
#' @export
plot_phenotypes <- function(pheno, traits = NULL, colors = NULL,
                            notch = TRUE, show_mean = TRUE,
                            main = NULL, ylab = "Days", ...) {
  
  if (inherits(pheno, "GBLUPData")) pheno <- pheno@phenotypes
  if (!is.data.frame(pheno)) stop("pheno must be a data frame or GBLUPData object")
  
  if (is.null(traits)) {
    meta_cols <- c("GID", "ENV", "ID", "id", "Genotype", "Environment", "Rep", "Block")
    numeric_cols <- names(pheno)[sapply(pheno, is.numeric)]
    traits <- setdiff(numeric_cols, meta_cols)
  }
  
  if (is.null(colors)) {
    colors <- c("#6BAED6", "#74C476", "#FD8D3C", "#FDDC00",
                "#9E9AC8", "#FC9272", "#A1D99B", "#9ECAE1")[seq_along(traits)]
  }
  
  if (is.null(main)) main <- paste("Boxplots of", paste(traits, collapse = ", "))
  
  plot_list <- lapply(traits, function(tr) pheno[[tr]])
  names(plot_list) <- traits
  
  bp <- boxplot(plot_list, main = main, ylab = ylab, col = colors,
                notch = notch, las = 1, outline = TRUE, outpch = 1,
                outcex = 0.7, boxwex = 0.6, ...)
  
  if (show_mean) {
    means <- sapply(plot_list, mean, na.rm = TRUE)
    points(seq_along(means), means, pch = 16, col = "red", cex = 1.3)
  }
  
  invisible(bp)
}


#' Plot Observed vs Predicted
#'
#' @param model GBLUPModel object from fit_gblup()
#' @param by_env Logical, create separate plots by environment
#' @param add_line Add regression line
#' @param add_identity Add 1:1 identity line
#' @param point_color Point color
#' @param main Plot title
#' @param ... Additional arguments passed to plot
#' @return Invisible NULL
#' @export
plot_obs_vs_pred <- function(model, by_env = FALSE, add_line = TRUE,
                             add_identity = TRUE, point_color = "#3182BD",
                             main = NULL, ...) {
  
  if (!inherits(model, "GBLUPModel")) stop("model must be a GBLUPModel object")
  gebv <- model@gebv
  if (is.null(main)) main <- "Observed vs Predicted"
  
  if (by_env && "ENV" %in% names(gebv)) {
    envs <- unique(gebv$ENV)
    colors <- c("#3182BD", "#E6550D", "#31A354", "#756BB1", 
                "#636363", "#6BAED6", "#FD8D3C", "#74C476")
    par(mfrow = c(ceiling(length(envs)/2), min(2, length(envs))), mar = c(4, 4, 3, 1))
    
    for (i in seq_along(envs)) {
      env_data <- gebv[gebv$ENV == envs[i], ]
      env_data <- env_data[complete.cases(env_data[, c("GEBV", "observed")]), ]
      if (nrow(env_data) < 2) next
      
      plot(env_data$GEBV, env_data$observed, xlab = "Predicted (GEBV)", ylab = "Observed",
           main = paste("Environment:", envs[i]), pch = 16,
           col = adjustcolor(colors[i], alpha.f = 0.6), cex = 1.2, las = 1, ...)
      
      if (add_identity) abline(0, 1, col = "gray50", lty = 2, lwd = 2)
      if (add_line) abline(lm(observed ~ GEBV, data = env_data), col = colors[i], lwd = 2)
      
      r <- cor(env_data$GEBV, env_data$observed, use = "complete.obs")
      legend("topleft", legend = paste("r =", round(r, 3)), bty = "n", cex = 1.1)
    }
    par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)
  } else {
    gebv_complete <- gebv[complete.cases(gebv[, c("GEBV", "observed")]), ]
    
    plot(gebv_complete$GEBV, gebv_complete$observed, xlab = "Predicted (GEBV)",
         ylab = "Observed", main = main, pch = 16,
         col = adjustcolor(point_color, alpha.f = 0.6), cex = 1.2, las = 1, ...)
    
    if (add_identity) abline(0, 1, col = "gray50", lty = 2, lwd = 2)
    if (add_line) abline(lm(observed ~ GEBV, data = gebv_complete), col = point_color, lwd = 2)
    
    r <- cor(gebv_complete$GEBV, gebv_complete$observed, use = "complete.obs")
    legend("topleft", legend = paste("Accuracy =", round(r, 3)), bty = "n", cex = 1.2)
  }
  
  invisible(NULL)
}


#' Plot Variance Components
#' @param model GBLUPModel object
#' @param type Type: "barplot" or "pie"
#' @param colors Color palette
#' @param main Plot title
#' @param ... Additional arguments
#' @return Invisible NULL
#' @export
plot_variance_components <- function(model, type = c("barplot", "pie"),
                                     colors = NULL, main = NULL, ...) {
  type <- match.arg(type)
  if (!inherits(model, "GBLUPModel")) stop("model must be a GBLUPModel object")
  
  vc <- model@varcomp
  if (nrow(vc) == 0) stop("No variance components found")
  if (is.null(main)) main <- "Variance Components"
  
  effect_names <- if ("Component" %in% names(vc)) vc$Component else paste("Effect", seq_len(nrow(vc)))
  
  if (is.null(colors)) {
    colors <- c("#3182BD", "#E6550D", "#31A354", "#756BB1", 
                "#636363", "#6BAED6", "#FD8D3C", "#74C476")[seq_len(nrow(vc))]
  }
  
  if (type == "barplot") {
    par(mar = c(7, 4, 4, 2) + 0.1)
    bp <- barplot(vc$Variance, names.arg = effect_names, main = main,
                  ylab = "Variance", col = colors, las = 2, border = "white", ...)
    
    total_var <- sum(vc$Variance)
    if (total_var > 0) {
      percentages <- round(vc$Variance / total_var * 100, 1)
      text(x = bp, y = vc$Variance + max(vc$Variance) * 0.05,
           labels = paste0(percentages, "%"), pos = 3, cex = 0.9, font = 2)
    }
    par(mar = c(5, 4, 4, 2) + 0.1)
  } else {
    percentages <- round(vc$Variance / sum(vc$Variance) * 100, 1)
    labels <- paste0(effect_names, "\n(", percentages, "%)")
    pie(vc$Variance, labels = labels, main = main, col = colors, border = "white", ...)
  }
  
  invisible(NULL)
}


#' Plot Heritability
#' @param h2 Named numeric vector or GBLUPModel
#' @param se Optional standard errors
#' @param colors Bar colors
#' @param main Plot title
#' @param ylim Y-axis limits
#' @param ... Additional arguments
#' @return Invisible NULL
#' @export
plot_heritability <- function(h2, se = NULL, colors = NULL, main = "Heritability",
                              ylim = c(0, 1), ...) {
  if (inherits(h2, "GBLUPModel")) h2 <- heritability(h2)
  if (!is.numeric(h2)) stop("h2 must be numeric")
  
  # Remove n_env if present (it's metadata, not a heritability value)
  h2 <- h2[!names(h2) %in% c("n_env")]
  
  # Remove NA values
  h2 <- h2[!is.na(h2)]
  
  if (length(h2) == 0) stop("No valid heritability values to plot")
  
  if (is.null(colors)) colors <- rep("#3182BD", length(h2))
  
  par(mar = c(7, 4, 4, 2) + 0.1)
  bp <- barplot(h2, main = main, ylab = expression(h^2), ylim = ylim,
                col = colors, las = 2, border = "white", ...)
  
  if (!is.null(se) && length(se) == length(h2)) {
    arrows(bp, h2 - se, bp, h2 + se, angle = 90, code = 3, length = 0.1, lwd = 2)
  }
  
  abline(h = c(0.3, 0.5, 0.7), lty = 3, col = "gray60")
  text(bp, h2 + 0.03, labels = round(h2, 2), pos = 3, cex = 0.9, font = 2)
  par(mar = c(5, 4, 4, 2) + 0.1)
  
  invisible(NULL)
}


#' Plot Relationship Matrix Heatmap
#' @param K Relationship matrix
#' @param cluster Cluster individuals
#' @param main Plot title
#' @param col Color palette
#' @param max_ind Maximum individuals to plot
#' @param ... Additional arguments
#' @return Invisible NULL
#' @export
plot_relationship_matrix <- function(K, cluster = TRUE, main = "Genomic Relationship Matrix",
                                     col = NULL, max_ind = 100, ...) {
  if (!is.matrix(K)) stop("K must be a matrix")
  
  if (is.null(col)) {
    col <- colorRampPalette(c("#2166AC", "#4393C3", "#92C5DE", "#F7F7F7", 
                              "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"))(100)
  }
  
  if (nrow(K) > max_ind) {
    K <- K[1:max_ind, 1:max_ind]
    main <- paste(main, "(first", max_ind, "individuals)")
  }
  
  if (cluster) {
    hc <- hclust(as.dist(1 - K), method = "ward.D2")
    K <- K[hc$order, hc$order]
  }
  
  par(mar = c(1, 1, 4, 5))
  image(1:nrow(K), 1:ncol(K), K, col = col, xlab = "", ylab = "", main = main, axes = FALSE, ...)
  par(mar = c(5, 4, 4, 2) + 0.1)
  
  invisible(NULL)
}


#' Plot Selection Response
#' @param gebv Data frame with GEBV column
#' @param top_percent Percentage of top individuals
#' @param main Plot title
#' @param ... Additional arguments
#' @return Invisible NULL
#' @export
plot_selection_response <- function(gebv, top_percent = 0.2, main = "Selection Response", ...) {
  if (!is.data.frame(gebv) || !"GEBV" %in% names(gebv)) {
    stop("gebv must be a data frame with 'GEBV' column")
  }
  
  if ("GID" %in% names(gebv)) {
    unique_gebv <- aggregate(GEBV ~ GID, data = gebv, mean, na.rm = TRUE)
  } else {
    unique_gebv <- data.frame(GID = seq_len(nrow(gebv)), GEBV = gebv$GEBV)
  }
  
  unique_gebv <- unique_gebv[order(-unique_gebv$GEBV), ]
  unique_gebv <- unique_gebv[!is.na(unique_gebv$GEBV), ]
  
  n_select <- max(1, ceiling(nrow(unique_gebv) * top_percent))
  threshold <- unique_gebv$GEBV[n_select]
  bar_colors <- ifelse(unique_gebv$GEBV >= threshold, "#E41A1C", "#969696")
  
  plot(seq_len(nrow(unique_gebv)), unique_gebv$GEBV, type = "h",
       xlab = "Individuals (ranked)", ylab = "GEBV", main = main,
       col = bar_colors, lwd = 2, las = 1, ...)
  
  abline(h = threshold, col = "#377EB8", lty = 2, lwd = 2)
  
  mean_selected <- mean(unique_gebv$GEBV[1:n_select])
  mean_population <- mean(unique_gebv$GEBV)
  
  abline(h = mean_selected, col = "#E41A1C", lty = 2, lwd = 2)
  abline(h = mean_population, col = "black", lty = 2, lwd = 2)
  
  legend("topright",
         legend = c(paste("Top", round(top_percent * 100), "%"),
                    paste("Mean selected:", round(mean_selected, 2)),
                    paste("Population mean:", round(mean_population, 2)),
                    paste("Sel. diff:", round(mean_selected - mean_population, 2))),
         col = c("#E41A1C", "#E41A1C", "black", NA),
         lty = c(1, 2, 2, 0), lwd = 2, bty = "n", cex = 0.9)
  
  invisible(NULL)
}


#' Plot PCA (2D or 3D)
#' @param gdata GBLUPData object or genotype/relationship matrix
#' @param K Optional relationship matrix
#' @param groups Grouping factor for coloring
#' @param colors Named vector of colors
#' @param pcs Which PCs to plot
#' @param plot_3d Create 3D plot
#' @param point_size Point size
#' @param show_labels Show point labels
#' @param title Plot title
#' @param ... Additional arguments
#' @return List with PCA results
#' @export
plot_pca <- function(gdata, K = NULL, groups = NULL, colors = NULL,
                     pcs = c(1, 2), plot_3d = FALSE, point_size = 1.5,
                     show_labels = FALSE, title = NULL, ...) {
  
  if (inherits(gdata, "GBLUPData")) {
    M <- gdata@genotypes
  } else if (is.matrix(gdata)) {
    M <- gdata
  } else {
    stop("gdata must be a GBLUPData object or matrix")
  }
  
  is_K <- !is.null(K) || (nrow(M) == ncol(M) && isSymmetric(M))
  
  if (!is.null(K) || is_K) {
    if (is.null(K)) K <- M
    eig <- eigen(K, symmetric = TRUE)
    scores <- eig$vectors %*% diag(sqrt(pmax(eig$values, 0)))
    var_explained <- (eig$values / sum(eig$values)) * 100
    ids <- rownames(K)
  } else {
    M_scaled <- scale(M, center = TRUE, scale = FALSE)
    svd_result <- svd(M_scaled, nu = max(pcs), nv = 0)
    scores <- svd_result$u %*% diag(svd_result$d[1:ncol(svd_result$u)])
    var_explained <- (svd_result$d^2 / sum(svd_result$d^2)) * 100
    ids <- rownames(M)
  }
  
  if (is.null(ids)) ids <- paste0("Ind", seq_len(nrow(scores)))
  
  if (!is.null(groups)) {
    if (is.null(colors)) {
      unique_groups <- unique(groups)
      colors <- setNames(
        c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
          "#FFFF33", "#A65628", "#F781BF")[seq_along(unique_groups)],
        unique_groups)
    }
    point_colors <- colors[groups]
  } else {
    point_colors <- rep("red", nrow(scores))
  }
  
  pc_labels <- paste0("PC", pcs, " (", round(var_explained[pcs], 1), "%)")
  if (is.null(title)) title <- "Principal Component Analysis"
  
  if (plot_3d && length(pcs) >= 3) {
    if (requireNamespace("scatterplot3d", quietly = TRUE)) {
      s3d <- scatterplot3d::scatterplot3d(
        x = scores[, pcs[1]], y = scores[, pcs[2]], z = scores[, pcs[3]],
        color = point_colors, pch = 16, cex.symbols = point_size,
        xlab = pc_labels[1], ylab = pc_labels[2], zlab = pc_labels[3],
        main = title, box = TRUE, grid = TRUE, angle = 40)
      
      if (show_labels) {
        s3d.coords <- s3d$xyz.convert(scores[, pcs[1]], scores[, pcs[2]], scores[, pcs[3]])
        text(s3d.coords$x, s3d.coords$y, labels = ids, pos = 4, cex = 0.6)
      }
      if (!is.null(groups)) legend("topright", legend = names(colors), col = colors, pch = 16, bty = "n")
      
      return(invisible(list(scores = scores, var_explained = var_explained, ids = ids)))
    } else {
      message("Package 'scatterplot3d' required for 3D. Using 2D.")
    }
  }
  
  plot(scores[, pcs[1]], scores[, pcs[2]], xlab = pc_labels[1], ylab = pc_labels[2],
       main = title, pch = 16, col = point_colors, cex = point_size, las = 1, ...)
  
  if (show_labels) text(scores[, pcs[1]], scores[, pcs[2]], labels = ids, pos = 4, cex = 0.6)
  if (!is.null(groups)) legend("topright", legend = names(colors), col = colors, pch = 16, bty = "n")
  
  invisible(list(scores = scores, var_explained = var_explained, ids = ids))
}


#' Plot Population Structure Dendrogram
#' @param gdata GBLUPData object or relationship matrix
#' @param K Optional relationship matrix
#' @param groups Grouping factor for coloring
#' @param colors Named vector of colors
#' @param type "circular" or "rectangular"
#' @param method Clustering method
#' @param show_labels Show tip labels
#' @param title Plot title
#' @return Invisible hclust object
#' @export
plot_dendrogram <- function(gdata, K = NULL, groups = NULL, colors = NULL,
                            type = c("circular", "rectangular"),
                            method = "ward.D2", show_labels = TRUE,
                            title = "Population Structure") {
  type <- match.arg(type)
  
  if (is.null(K)) {
    if (inherits(gdata, "GBLUPData")) {
      M <- gdata@genotypes
      M_scaled <- scale(M, center = TRUE, scale = FALSE)
      K <- tcrossprod(M_scaled) / ncol(M)
    } else if (is.matrix(gdata)) {
      K <- gdata
    } else {
      stop("Provide GBLUPData object or relationship matrix")
    }
  }
  
  K_scaled <- (K - min(K)) / (max(K) - min(K))
  D <- as.dist(1 - K_scaled)
  hc <- hclust(D, method = method)
  
  if (!is.null(groups)) {
    if (is.null(colors)) {
      unique_groups <- unique(groups)
      colors <- setNames(
        c("#6B8E9F", "#E8A838", "#C4A3BF", "#F5C6C6", "#B8D4B8",
          "#FFE4B5", "#DDA0DD", "#87CEEB")[seq_along(unique_groups)],
        unique_groups)
    }
    label_colors <- colors[groups[hc$order]]
  } else {
    label_colors <- "black"
  }
  
  if (type == "circular" && requireNamespace("ape", quietly = TRUE)) {
    phylo <- ape::as.phylo(hc)
    par(mar = c(1, 1, 2, 1))
    ape::plot.phylo(phylo, type = "fan", tip.color = label_colors,
                    cex = 0.5, show.tip.label = show_labels, main = title)
    if (!is.null(groups)) {
      legend("topright", legend = names(colors), fill = colors, border = NA, bty = "n", cex = 0.7)
    }
    par(mar = c(5, 4, 4, 2) + 0.1)
  } else {
    if (type == "circular") message("Package 'ape' required for circular dendrogram. Using rectangular.")
    par(mar = c(2, 4, 2, 1))
    plot(hc, main = title, xlab = "", sub = "", labels = if (show_labels) NULL else FALSE)
    if (!is.null(groups)) {
      legend("topright", legend = names(colors), fill = colors, border = NA, bty = "n", cex = 0.7)
    }
    par(mar = c(5, 4, 4, 2) + 0.1)
  }
  
  invisible(hc)
}
