#' @title plot_manhattan
#' @description Generate a Manhattan plot from a PLINK assoc.logistic data.frame
#' @param df A data.frame containing at least the following columns: 
#'           - CHR: chromosome identifier  
#'           - BP: base-pair position  
#'           - P: p-value  
#' @param chr_col Name of the chromosome column (string). Default: "CHR"
#' @param bp_col Name of the base-pair position column (string). Default: "BP"
#' @param p_col  Name of the p-value column (string). Default: "P"
#' @param sig_line Significance threshold for p-values. Default: 5e-8
#' @param title   Plot title. Default: "Manhattan Plot"
#' @return A ggplot2 object representing the Manhattan plot
#' @export
plot_manhattan <- function(df,
                           chr_col  = "CHR",
                           bp_col   = "BP",
                           p_col    = "P",
                           sig_line = 5e-8,
                           title    = "Manhattan Plot") {
  
  # 1. Filter out rows with missing p-values
  df2 <- df %>%
    filter(!is.na(.data[[p_col]]))
  
  # 2. Compute chromosome lengths and cumulative start positions
  chr_info <- df2 %>%
    group_by(chr = .data[[chr_col]]) %>%
    summarise(chr_len = max(.data[[bp_col]]), .groups = "drop") %>%
    arrange(chr) %>%
    mutate(cum_start = cumsum(lag(chr_len, default = 0)))
  
  # 3. Join cumulative start and compute cumulative SNP positions
  df2 <- df2 %>%
    left_join(chr_info %>% select(chr, cum_start),
              by = setNames("chr", chr_col)) %>%
    arrange(.data[[chr_col]], .data[[bp_col]]) %>%
    mutate(pos_cum = .data[[bp_col]] + cum_start)
  
  # 4. Prepare axis labels (chromosome centers)
  axis_df <- df2 %>%
    group_by(chr = .data[[chr_col]]) %>%
    summarise(center = (min(pos_cum) + max(pos_cum)) / 2, .groups = "drop")
  
  y_data_max <- max(-log10(df2[[p_col]]), na.rm = TRUE)
  y_sig      <- -log10(sig_line)
  y_upper    <- max(y_data_max, y_sig, 10)
  # 5. Build the Manhattan plot
  p <- ggplot(df2, aes(x = pos_cum, y = -log10(.data[[p_col]]),
                       color = as.factor(.data[[chr_col]]))) +
    geom_point(alpha = 0.6, size = 0.8) +
    scale_x_continuous(breaks = axis_df$center,
                       labels = axis_df$chr) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(ylim = c(0, y_upper)) +
    labs(
      x     = "Chromosome",
      y     = expression(-log[10](P)),
      title = title
    ) +
    geom_hline(yintercept = -log10(sig_line),
               linetype = "dashed", color = "red") +
    theme_bw() +
    theme(
      legend.position        = "none",
      panel.grid.major.x     = element_blank(),
      panel.grid.minor.x     = element_blank()
    )
  
  return(p)
}


#' @title Prepare Data for GWAS Manhattan Plot
#' @description Process a GWAS result data.frame to compute cumulative SNP positions, -log10(p-values),
#'              and interactive tooltips, as well as x- and y-axis ranges for an interactive Manhattan plot.
#'
#' @param df A data.frame containing at least the following columns:
#'   - CHR: chromosome identifier
#'   - BP: base-pair position
#'   - P: p-value
#' @param chr_col Name of the chromosome column (string). Default: "CHR"
#' @param bp_col Name of the base-pair position column (string). Default: "BP"
#' @param p_col  Name of the p-value column (string). Default: "P"
#'
#' @return A list with components:
#'   - df: The input data.frame augmented with `pos_cum`, `logP`, and `tooltip` columns.
#'   - axis_df: A data.frame of chromosome centers for x-axis breaks (`chr`, `center`).
#'   - xrange: Numeric vector of length 2 specifying the cumulative position range.
#'   - yrange: Numeric vector of length 2 specifying the y-axis plotting range.
#'
#' @examples
#' # Synthetic example
#' df_example <- data.frame(
#'   CHR = rep(1:3, each = 5),
#'   BP  = rep(seq(1e6, 5e6, by = 1e6), 3),
#'   P   = runif(15, 1e-9, 1)
#' )
#' res <- prep_manhattan(df_example)
#' str(res)
#'
#' # Real GWAS dataset example
#' assoc <- read.table("output_gwas_results.assoc.logistic",
#'                     header = TRUE, stringsAsFactors = FALSE)
#'
#' # Prepare tooltip
#' assoc$tooltip <- paste0(
#'   "CHR: ", assoc$CHR,
#'   "\\nBP: ", assoc$BP,
#'   "\\nP: ", signif(assoc$P, 3),
#'   if ("SNP" %in% names(assoc)) paste0("\\nSNP: ", assoc$SNP) else ""
#' )
#'
#' # Process cumulative position
#' manhattan_data <- prep_manhattan(assoc)
#'
#' # Plot using a basic Manhattan plot function
#' manhattan_plot <- plot_manhattan_bg(assoc)
#' print(manhattan_plot)
#'
#' # Save to file
#' ggsave("manhattan_plot.png", manhattan_plot, width = 10, height = 5, dpi = 300)
#'
#' @export
prep_manhattan <- function(df,
                           chr_col = "CHR",
                           bp_col  = "BP",
                           p_col   = "P") {
  df2 <- df %>% filter(!is.na(.data[[p_col]]))
  chr_info <- df2 %>%
    group_by(chr = .data[[chr_col]]) %>%
    summarise(chr_len = max(.data[[bp_col]]), .groups = "drop") %>%
    arrange(chr) %>%
    mutate(cum_start = cumsum(as.numeric(lag(chr_len, default = 0))))
  df2 <- df2 %>%
    left_join(chr_info %>% select(chr, cum_start),
              by = setNames("chr", chr_col)) %>%
    arrange(.data[[chr_col]], .data[[bp_col]]) %>%
    mutate(
      pos_cum = as.numeric(.data[[bp_col]]) + as.numeric(cum_start),
      logP    = -log10(.data[[p_col]]),
      tooltip = paste0(
        "CHR: ", .data[[chr_col]],
        "<br>BP: ", .data[[bp_col]],
        "<br>P: ", signif(.data[[p_col]], 3),
        if ("SNP" %in% names(.)) paste0("<br>SNP: ", .data[["SNP"]]) else ""
      )
    )
  axis_df <- df2 %>%
    group_by(chr = .data[[chr_col]]) %>%
    summarise(center = (min(pos_cum) + max(pos_cum)) / 2, .groups = "drop")
  list(
    df       = df2,
    axis_df  = axis_df,
    xrange   = range(df2$pos_cum, na.rm = TRUE),
    yrange   = c(0, max(df2$logP, na.rm = TRUE) * 1.05)
  )
}




