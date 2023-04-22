#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------- #
# rx_identifier.R - Compute an Rx ratio from a single individual, using the    #
#                   output of samtools idxstats.                               #
# ---------------------------------------------------------------------------- #
# Author: Maël Lefeuvre                                                        #
# ---------------------------------------------------------------------------- #
# - Based on the methods developed in (Mitnik A, et al - 2016)                 #
#   https://doi.org/10.1371/journal.pone.0163019                               #
#                                                                              #
# - Original script: https://doi.org/10.1371/journal.pone.0163019.s003         #
# ---------------------------------------------------------------------------- #
# Usage: rx_identifier.R <output-prefix> [input.flagstats]

rx_to_classification <- function(lower_tail, upper_tail) {
  if (lower_tail > 0.8) {
    "XX"
  } else if (upper_tail < 0.6) {
    "XY"
  } else if (lower_tail > 0.6 && upper_tail > 0.8) {
    "XX_consistent/XY_inconsistent"
  } else if (lower_tail < 0.6 && upper_tail < 0.8) {
    "XY_consistent/XX_inconsistent"
  } else {
    "Unassignable"
  }
}

#' Create a linear model between chromosome length and mapped reads.
#' Return the R squared coefficient + the model's p.value.
get_mapping_significance <- function(idxstats) {
  model     <- summary(lm(idxstats$sequence.length ~ idxstats$mapped))
  list(r.squared = model$r.squared, p.value = model$coefficients[2, 4])
}


find_one_chromosome_idx <- function(idxstats, regex) {
  chr_index <- grep(regex, idxstats$chr)
  if (length(chr_index) > 1) {
    msg <- paste(
      "Found multiple chromosome candidates when using regex pattern",
      regex
    )
    stop(msg)
  }
  chr_index
}

parse_idxstats_file <- function(path) {
  # ---- Open idxstats file.
  idxstats <- read.table(
    path,
    header = FALSE,
    col.names = c("chr", "sequence.length", "mapped", "unmapped")
  )

  # ---- Keep only autosomes + X + Y, and sort lexicographically
  regex    <- "^(chr){0,1}([0-9]+)|X|Y"
  idxstats <- idxstats[grepl(regex, idxstats$chr, perl = TRUE), ]
  idxstats[order(idxstats$chr), ]
}

idxstats_to_rx <- function(idxstats) {

  total_length <- sum(as.numeric(idxstats[, 2]))
  total_mapped <- sum(as.numeric(idxstats[, 3]))


  # ---- Compute per-chromosome mapped vs. length ratio.
  idxstats$Rt <- (
    (idxstats[, 3] / total_mapped) / (idxstats[, 2] / total_length)
  )

  # ---- Find the index of chromosome X and Y
  chr_x_idx <- find_one_chromosome_idx(idxstats, regex = "^(chr){0,1}(23|X)$")
  chr_y_idx <- find_one_chromosome_idx(idxstats, regex = "^(chr){0,1}(23|Y)$")

  # ---- Compute Rt_autosomes / Rt_X
  norm_rx <- (
    idxstats$Rt[chr_x_idx] / idxstats$Rt[-c(chr_x_idx, chr_y_idx)]
  )

  if (length(norm_rx) > 22) {
    msg <- paste(
      "Vector of normalized Rxs is greater than the expected",
      "number of autosomes. Expected", 22,
      ". Got", length(norm_rx)
    )
    stop(msg)
  }

  average_rx <- mean(norm_rx)

  # ---- Compute the confidence interval
  ci_95      <- 1.96 * (sd(norm_rx) / sqrt(length(norm_rx)))
  lower_tail <- average_rx - ci_95
  upper_tail <- average_rx + ci_95

  # ---- Assign XY / XX classification, based on CI thresholds
  sex_assignment <- rx_to_classification(lower_tail, upper_tail)

  # ---- Ensure there is a correlation between mapped reads and correlations
  mapping_correlation <- get_mapping_significance(idxstats)
  if (mapping_correlation$p.value > 0.05) {

    msg <- paste(
      "Lack of correlation between chromosome length and amount of mapped",
      "reads."
    )
    warning(msg)
  }

  list(
    Rx                        = average_rx,
    confidence_interval       = ci_95,
    assignment                = sex_assignment,
    mapping_corr_p_val        = mapping_correlation$p.value,
    mapping_corr_R_squared    = mapping_correlation$r.squared
  )
}

ry_to_classification <- function(
  ry, ci, female_threshold = 0.016, male_threshold = 0.075
) {
  if (ry < female_threshold && ry > male_threshold) {
    "Unassignable"
  } else if (ry == 0.0) {
    "XX_consistent"
  } else if (ry + ci < female_threshold) {
    "XX"
  } else if (ry - ci > male_threshold) {
    "XY"
  } else if (ry - ci > female_threshold && ry + ci > male_threshold) {
    "XY_consistent/XX_inconsistent"
  } else if (ry - ci < female_threshold && ry + ci < male_threshold) {
    "XX_consistent/XY_inconsistent"
  } else {
    "Unassignable"
  }

}

# Inspired Pontus Skoglund
idxstats_to_ry <- function(idxstats) {

  # Search for the index of chromosomes X and Y
  chr_x_idx <- find_one_chromosome_idx(idxstats, regex = "^(chr){0,1}(23|X)$")
  chr_y_idx <- find_one_chromosome_idx(idxstats, regex = "^(chr){0,1}(23|Y)$")

  # Compute R_y with 95% confidence interval
  n                   <- idxstats$mapped[chr_x_idx] + idxstats$mapped[chr_y_idx]
  ry                  <- 1.0 * idxstats$mapped[chr_y_idx] / n
  standard_error      <- sqrt((ry * (1.0 - ry)) / n)
  confidence_interval <- 1.96 * standard_error

  list(
    n_mapped   = sum(idxstats$mapped),
    n_XY_seq   = n,
    n_Y        = idxstats$mapped[chr_y_idx],
    Ry         = ry,
    ci         = confidence_interval,
    Assignment = ry_to_classification(ry, confidence_interval)
  )
}


if (!interactive()) {

  args          <- commandArgs(trailingOnly = TRUE)
  output_prefix <- args[1]
  input_files   <- args[2:length(args)]

  idxstats_dfs  <- sapply(
    X         = input_files,
    simplify  = FALSE,
    USE.NAMES = TRUE,
    FUN       = parse_idxstats_file
  )

  rx <- lapply(input_files, FUN = function(path) {
    sample   <- basename(path)
    idxstats <- idxstats_dfs[[path]]
    cbind(sample = as.data.frame(sample), idxstats_to_rx(idxstats))
  })

  rx <- do.call("rbind", rx)

  ry <- lapply(input_files, FUN = function(path) {
    sample   <- basename(path)
    idxstats <- idxstats_dfs[[path]]
    cbind(sample = as.data.frame(sample), idxstats_to_ry(idxstats))
  })
  ry <- do.call("rbind", ry)


  write.table(
    x      = rx,
    file   = paste0(output_prefix, "-Rx-ratios.tsv"),
    quote  = FALSE,
    sep    = "\t",
    row.names = FALSE
  )

  write.table(
    x      = ry,
    file   = paste0(output_prefix, "-Ry-ratios.tsv"),
    quote  = FALSE,
    sep    = "\t",
    row.names = FALSE
  )
}
