vec2mat <- function(vec, n_times) {
  return(matrix(rep(vec, n_times), ncol = n_times, byrow = FALSE))
}

match_quantiles <- function(counts_sub, old_mu, old_phi, new_mu, new_phi) {
  new_counts_sub <- matrix(NA, nrow = nrow(counts_sub), ncol = ncol(counts_sub))
  for (a in 1:nrow(counts_sub)) {
    for (b in 1:ncol(counts_sub)) {
      if (counts_sub[a, b] <= 1) {
        new_counts_sub[a, b] <- counts_sub[a, b]
      }
      else {
        tmp_p <- pnbinom(counts_sub[a, b] - 1, mu = old_mu[
          a,
          b
        ], size = 1 / old_phi[a])
        if (abs(tmp_p - 1) < 1e-04) {
          new_counts_sub[a, b] <- counts_sub[a, b]
        }
        else {
          new_counts_sub[a, b] <- 1 + qnbinom(tmp_p,
            mu = new_mu[a, b], size = 1 / new_phi[a]
          )
        }
      }
    }
  }
  return(new_counts_sub)
}

ComBat_seq_custom <- function(counts, batch, group = NULL, covar_mod = NULL, full_mod = TRUE,
                              shrink = FALSE, shrink.disp = FALSE, gene.subset.n = NULL) {
  batch <- as.factor(batch)
  if (any(table(batch) <= 1)) {
    stop("ComBat-seq doesn't support 1 sample per batch yet")
  }
  keep_lst <- lapply(levels(batch), function(b) {
    which(apply(counts[, batch == b], 1, function(x) {
      !all(x == 0)
    }))
  })
  keep <- Reduce(intersect, keep_lst)
  rm <- setdiff(1:nrow(counts), keep)
  countsOri <- counts
  counts <- counts[keep, ]
  dge_obj <- DGEList(counts = counts)
  n_batch <- nlevels(batch)
  batches_ind <- lapply(1:n_batch, function(i) {
    which(batch == levels(batch)[i])
  })
  n_batches <- sapply(batches_ind, length)
  n_sample <- sum(n_batches)
  cat("Found", n_batch, "batches\n")
  batchmod <- model.matrix(~ -1 + batch)
  group <- as.factor(group)
  if (full_mod & nlevels(group) > 1) {
    cat("Using full model in ComBat-seq.\n")
    mod <- model.matrix(~group)
  } else {
    cat("Using null model in ComBat-seq.\n")
    mod <- model.matrix(~1, data = as.data.frame(t(counts)))
  }
  if (!is.null(covar_mod)) {
    if (is.data.frame(covar_mod)) {
      covar_mod <- do.call(cbind, lapply(
        1:ncol(covar_mod),
        function(i) {
          model.matrix(~ covar_mod[, i])
        }
      ))
    }
    covar_mod <- covar_mod[, !apply(covar_mod, 2, function(x) {
      all(x == 1)
    })]
  }
  mod <- cbind(mod, covar_mod)
  design <- cbind(batchmod, mod)
  check <- apply(design, 2, function(x) all(x == 1))
  design <- as.matrix(design[, !check])
  cat("Adjusting for", ncol(design) - ncol(batchmod), "covariate(s) or covariate level(s)\n")
  if (qr(design)$rank < ncol(design)) {
    if (ncol(design) == (n_batch + 1)) {
      stop("The covariate is confounded with batch! Remove the covariate and rerun ComBat-Seq")
    }
    if (ncol(design) > (n_batch + 1)) {
      if ((qr(design[, -c(1:n_batch)])$rank < ncol(design[
        ,
        -c(1:n_batch)
      ]))) {
        stop("The covariates are confounded! Please remove one or more of the covariates so the design is not confounded")
      } else {
        stop("At least one covariate is confounded with batch! Please remove confounded covariates and rerun ComBat-Seq")
      }
    }
  }
  NAs <- any(is.na(counts))
  if (NAs) {
    cat(c("Found", sum(is.na(counts)), "Missing Data Values\n"),
      sep = " "
    )
  }
  cat("Estimating dispersions\n")
  disp_common <- sapply(1:n_batch, function(i) {
    if ((n_batches[i] <= ncol(design) - ncol(batchmod) +
      1) | qr(mod[batches_ind[[i]], ])$rank < ncol(mod)) {
      return(estimateGLMCommonDisp(counts[, batches_ind[[i]]],
        design = NULL, subset = nrow(counts)
      ))
    } else {
      return(estimateGLMCommonDisp(counts[, batches_ind[[i]]],
        design = mod[batches_ind[[i]], ], subset = nrow(counts)
      ))
    }
  })
  genewise_disp_lst <- lapply(1:n_batch, function(j) {
    if ((n_batches[j] <= ncol(design) - ncol(batchmod) +
      1) | qr(mod[batches_ind[[j]], ])$rank < ncol(mod)) {
      return(rep(disp_common[j], nrow(counts)))
    } else {
      return(estimateGLMTagwiseDisp(counts[, batches_ind[[j]]],
        design = mod[batches_ind[[j]], ], dispersion = disp_common[j],
        prior.df = 0
      ))
    }
  })
  names(genewise_disp_lst) <- paste0("batch", levels(batch))
  phi_matrix <- matrix(NA, nrow = nrow(counts), ncol = ncol(counts))
  for (k in 1:n_batch) {
    phi_matrix[, batches_ind[[k]]] <- vec2mat(
      genewise_disp_lst[[k]],
      n_batches[k]
    )
  }
  cat("Fitting the GLM model\n")
  glm_f <- glmFit(dge_obj,
    design = design, dispersion = phi_matrix,
    prior.count = 1e-04
  )
  alpha_g <- glm_f$coefficients[, 1:n_batch] %*% as.matrix(n_batches / n_sample)
  new_offset <- t(vec2mat(getOffset(dge_obj), nrow(counts))) +
    vec2mat(alpha_g, ncol(counts))
  glm_f2 <- glmFit.default(dge_obj$counts,
    design = design,
    dispersion = phi_matrix, offset = new_offset, prior.count = 1e-04
  )
  gamma_hat <- glm_f2$coefficients[, 1:n_batch]
  mu_hat <- glm_f2$fitted.values
  phi_hat <- do.call(cbind, genewise_disp_lst)
  if (shrink) {
    cat("Apply shrinkage - computing posterior estimates for parameters\n")
    mcint_fun <- monte_carlo_int_NB
    monte_carlo_res <- lapply(1:n_batch, function(ii) {
      if (ii == 1) {
        mcres <- mcint_fun(
          dat = counts[, batches_ind[[ii]]],
          mu = mu_hat[, batches_ind[[ii]]],
          gamma = gamma_hat[, ii],
          phi = phi_hat[, ii],
          gene.subset.n = gene.subset.n
        )
      } else {
        invisible(capture.output(mcres = mcint_fun(
          dat = counts[, batches_ind[[ii]]], mu = mu_hat[, batches_ind[[ii]]],
          gamma = gamma_hat[, ii], phi = phi_hat[, ii],
          gene.subset.n = gene.subset.n
        )))
      }
      return(mcres)
    })
    names(monte_carlo_res) <- paste0("batch", levels(batch))
    gamma_star_mat <- lapply(monte_carlo_res, function(res) {
      res$gamma_star
    })
    gamma_star_mat <- do.call(cbind, gamma_star_mat)
    phi_star_mat <- lapply(monte_carlo_res, function(res) {
      res$phi_star
    })
    phi_star_mat <- do.call(cbind, phi_star_mat)
    if (!shrink.disp) {
      cat("Apply shrinkage to mean only\n")
      phi_star_mat <- phi_hat
    }
  } else {
    cat("Shrinkage off - using GLM estimates for parameters\n")
    gamma_star_mat <- gamma_hat
    phi_star_mat <- phi_hat
  }
  mu_star <- matrix(NA, nrow = nrow(counts), ncol = ncol(counts))
  for (jj in 1:n_batch) {
    mu_star[, batches_ind[[jj]]] <- exp(log(mu_hat[, batches_ind[[jj]]]) -
      vec2mat(gamma_star_mat[, jj], n_batches[jj]))
  }
  phi_star <- rowMeans(phi_star_mat)
  cat("Adjusting the data\n")
  adjust_counts <- matrix(NA, nrow = nrow(counts), ncol = ncol(counts))
  for (kk in 1:n_batch) {
    cat("Progress: ",kk,"/",n_batch,"\r",sep = "")
    counts_sub <- counts[, batches_ind[[kk]]]
    old_mu <- mu_hat[, batches_ind[[kk]]]
    old_phi <- phi_hat[, kk]
    new_mu <- mu_star[, batches_ind[[kk]]]
    new_phi <- phi_star
    adjust_counts[, batches_ind[[kk]]] <- match_quantiles(
      counts_sub = counts_sub,
      old_mu = old_mu, old_phi = old_phi, new_mu = new_mu,
      new_phi = new_phi
    )
  }
  colnames(adjust_counts) <- colnames(countsOri)
  rownames(adjust_counts) <- rownames(countsOri)[keep]
  adjust_counts_whole <- rbind(adjust_counts, countsOri[rm, ])
  return(adjust_counts_whole)
}

cleanY <- function(y, mod, svs) {
  X <- cbind(mod, svs)
  Hat <- solve(t(X) %*% X) %*% t(X)
  beta <- (Hat %*% t(y))
  rm(Hat)
  gc()
  P <- ncol(mod)
  return(y - t(as.matrix(X[, -c(1:P)]) %*% beta[-c(1:P), ]))
}


custom_point <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    stat_density2d(geom = "tile", aes(fill = ..density..^0.25, alpha = 1), contour = FALSE, show.legend = F) +
    stat_density2d(geom = "tile", aes(fill = ..density..^0.25, alpha = ifelse(..density..^0.25 < 0.4, 0, 1)), contour = FALSE, show.legend = F) +
    geom_point(shape = 16, color = "#E18727FF", size = 1, alpha = 0.4) +
    geom_rug(alpha = 0.4) +
    geom_smooth(method = "lm", formula = "y ~ x", color = "red", alpha = 0.5) +
    geom_abline(intercept = 0, slope = 1, size = 1, alpha = 0.5, color = "black", linetype = 2) +
    scale_fill_gradientn(colours = colorRampPalette(c("white", "#BC3C29FF"))(256)) +
    # scale_x_continuous(trans = "log10", labels = trans_format("log10", math_format(10^.x)))+
    # scale_y_continuous(trans = "log10", labels = trans_format("log10", math_format(10^.x)))+
    theme_pubr(border = T) +
    theme(axis.text = element_text(size = 10))
}

custom_density <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_histogram(color = "black", fill = "#E18727FF", bins = 30) +
    # scale_x_continuous(trans = "log10", labels = trans_format("log10", math_format(10^.x)))+
    theme_pubr(border = T) +
    theme(axis.text = element_text(size = 10))
}

custom_cor <- function(data, mapping, ...) {
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  cor <- round(cor(x = x, y = y, method = "pearson"), 2)
  color <- color_cor[as.character(cor)]
  ggplot() +
    labs(x = NULL, y = NULL) +
    geom_text(aes(
      x = 0.5, y = 0.5,
      label = paste0("Corr:\n", cor)
    ),
    size = 5, fontface = 2, color = "white"
    ) +
    theme_pubr(border = T) +
    theme(panel.background = element_rect(fill = color))
}

aplotGrob <- function(x) {
  mp <- x$plotlist[[1]]
  if ( length(x$plotlist) == 1) {
    return(ggplotGrob(mp))
  }
  
  for (i in x$layout[, x$main_col]) {
    if (is.na(i)) next
    if (i == 1) next
    x$plotlist[[i]] <- suppressMessages(x$plotlist[[i]] + xlim2(mp))
  }
  for (i in x$layout[x$main_row,]) {
    if(is.na(i)) next
    if (i == 1) next
    x$plotlist[[i]] <- suppressMessages(x$plotlist[[i]] + ylim2(mp))
  }
  
  idx <- as.vector(x$layout)
  idx[is.na(idx)] <- x$n + 1 
  x$plotlist[[x$n+1]] <- ggplot() + theme_void() # plot_spacer()
  plotlist <- x$plotlist[idx]
  
  pp <- plotlist[[1]] + theme_no_margin()
  for (i in 2:length(plotlist)) {
    pp <- pp + (plotlist[[i]] + theme_no_margin())
  }
  
  res <- pp + plot_layout(byrow=F, ncol=ncol(x$layout), 
                          widths = x$width,
                          heights= x$height,
                          guides = 'collect')
  patchworkGrob(res)
}

theme_no_margin <- function(...) {
  ggplot2::theme(plot.margin = ggplot2::margin(), ...)
}
