###########################################################
# Correlation Tests Power Analysis                        #
# Robust implementation with proper dependencies          #
###########################################################

# Load required packages
library(future)
library(future.apply)
library(MASS)
library(mvtnorm)
library(copula)
library(extraDistr)
library(dplyr)
library(ggplot2)
library(kableExtra)
library(progress)
library(scales)
library(tibble)
library(purrr)
library(tidyr)
library(bench)

# Set seed for reproducibility
set.seed(2023)

# Configure parallel processing
plan(multisession, workers = max(1, availableCores() - 1))

## Section 2: Core Function Definitions ------------------------------------

#' Calculate correlation p-values for all three methods with error handling
#'
#' @param x Numeric vector
#' @param y Numeric vector
#' @return List with Pearson, Spearman, and Kendall p-values
#' @examples
#' calculate_cor_pvalues(rnorm(100), rnorm(100))
calculate_cor_pvalues <- function(x, y) {
  # Input validation
  stopifnot(is.numeric(x), is.numeric(y),
            length(x) == length(y),
            length(x) > 3)  # Minimum sample size for correlation
  
  # Safely compute correlations with error handling
  safe_cor_test <- function(method) {
    tryCatch({
      test <- cor.test(x, y, method = method)
      test$p.value
    }, error = function(e) NA_real_)
  }
  
  list(
    pearson = safe_cor_test("pearson"),
    spearman = safe_cor_test("spearman"),
    kendall = safe_cor_test("kendall")
  )
}

#' Perform power simulation for a given distribution with progress reporting
#'
#' @param n Sample size
#' @param rho True correlation parameter
#' @param n_sim Number of simulations
#' @param dist_fun Function that generates bivariate data
#' @param alpha Significance level (default 0.05)
#' @return Data frame with power estimates for each test
#' @examples
#' simulate_power(50, 0.3, 100, generate_bvn)
simulate_power <- function(n, rho, n_sim = 1000, dist_fun, alpha = 0.05) {
  # Input validation
  stopifnot(n > 3, 
            rho >= -1 && rho <= 1,
            n_sim > 0,
            is.function(dist_fun))
  
  # Use a global variable for progress tracking with future
  progress_file <- tempfile()
  cat(0, file = progress_file)
  
  # Run simulations
  results <- future.apply::future_replicate(n_sim, {
    # Update progress (in a way compatible with future)
    progress <- as.integer(readLines(progress_file))
    progress <- progress + 1
    cat(progress, file = progress_file)
    
    data <- tryCatch(dist_fun(n, rho), error = function(e) NULL)
    if (is.null(data)) return(c(pearson = NA, spearman = NA, kendall = NA))
    
    pvals <- calculate_cor_pvalues(data$x, data$y)
    c(pearson = pvals$pearson < alpha,
      spearman = pvals$spearman < alpha,
      kendall = pvals$kendall < alpha)
  }, simplify = "matrix", future.seed = TRUE)
  
  # Calculate power (proportion of rejections)
  power <- rowMeans(results, na.rm = TRUE)
  n_failed <- sum(is.na(results[1, ]))
  
  data.frame(
    test = c("Pearson", "Spearman", "Kendall"),
    power = power,
    n = n,
    rho = rho,
    distribution = deparse(substitute(dist_fun)),
    n_sim = n_sim,
    n_failed = n_failed,
    stringsAsFactors = FALSE
  )
}

## Section 3: Distribution Generators --------------------------------------

#' Generate bivariate normal data with parameter validation
#'
#' @param n Sample size
#' @param rho Correlation parameter
#' @return Data frame with x and y columns
#' @examples
#' generate_bvn(100, 0.5)
generate_bvn <- function(n, rho) {
  stopifnot(n > 0, rho >= -1 && rho <= 1)
  sigma <- matrix(c(1, rho, rho, 1), ncol = 2)
  data <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = sigma)
  data.frame(x = data[, 1], y = data[, 2])
}

#' Generate bivariate t-distribution data (heavy-tailed)
#'
#' @param n Sample size
#' @param rho Correlation parameter
#' @param df Degrees of freedom (default 3)
#' @return Data frame with x and y columns
#' @examples
#' generate_bvt(100, 0.5)
generate_bvt <- function(n, rho, df = 3) {
  stopifnot(n > 0, rho >= -1 && rho <= 1, df > 0)
  sigma <- matrix(c(1, rho, rho, 1), ncol = 2)
  data <- mvtnorm::rmvt(n, sigma = sigma, df = df)
  data.frame(x = data[, 1], y = data[, 2])
}

#' Generate data with quadratic relationship
#'
#' @param n Sample size
#' @param rho Base correlation parameter (affects noise)
#' @return Data frame with x and y columns
#' @examples
#' generate_quadratic(100, 0.2)
generate_quadratic <- function(n, rho) {
  stopifnot(n > 0, rho >= -1 && rho <= 1)
  x <- rnorm(n)
  y <- x^2 + rho * rnorm(n)
  data.frame(x = x, y = y)
}

#' Generate data with Frank copula (nonlinear dependence)
#'
#' @param n Sample size
#' @param rho Correlation parameter
#' @return Data frame with x and y columns
#' @examples
#' generate_frank_copula(100, 0.5)
generate_frank_copula <- function(n, rho) {
  stopifnot(n > 0, rho >= -1 && rho <= 1)
  tau <- rho  # Approximation
  theta <- copula::iTau(copula::frankCopula(), tau)
  cop <- copula::frankCopula(theta)
  data <- copula::rCopula(n, cop)
  data.frame(
    x = qnorm(data[, 1]),
    y = qnorm(data[, 2])
  )
}

#' Generate data with heavy-tailed margins
#'
#' @param n Sample size
#' @param rho Correlation parameter
#' @param df Degrees of freedom for t-distribution (default 3)
#' @return Data frame with x and y columns
#' @examples
#' generate_heavy_tails(100, 0.5)
generate_heavy_tails <- function(n, rho, df = 3) {
  stopifnot(n > 0, rho >= -1 && rho <= 1, df > 0)
  sigma <- matrix(c(1, rho, rho, 1), ncol = 2)
  normal_data <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = sigma)
  data.frame(
    x = extraDistr::qtpn(pnorm(normal_data[, 1]), df = df),
    y = extraDistr::qtpn(pnorm(normal_data[, 2]), df = df)
  )
}

## Section 4: Power Analysis Framework -------------------------------------

#' Run comprehensive power analysis across scenarios
#'
#' @param n_vec Vector of sample sizes to test
#' @param rho_vec Vector of correlation values to test
#' @param n_sim Number of simulations per scenario
#' @param distributions List of distribution functions to test
#' @return Data frame with power results
#' @examples
#' run_full_power_analysis(n_vec = c(20, 50), rho_vec = c(0.1, 0.3))
run_full_power_analysis <- function(n_vec = c(20, 50, 100),
                                   rho_vec = c(0.1, 0.3, 0.5),
                                   n_sim = 1000,
                                   distributions = list(
                                     bvn = generate_bvn,
                                     bvt = generate_bvt,
                                     quadratic = generate_quadratic,
                                     frank = generate_frank_copula,
                                     heavy = generate_heavy_tails
                                   )) {
  # Input validation
  stopifnot(all(n_vec > 3),
            all(rho_vec >= -1 & rho_vec <= 1),
            n_sim > 0,
            is.list(distributions),
            all(sapply(distributions, is.function)))
  
  # Create all combinations of parameters
  scenarios <- expand.grid(
    n = n_vec,
    rho = rho_vec,
    distribution = names(distributions),
    stringsAsFactors = FALSE
  )
  
  # Initialize result storage
  results <- vector("list", nrow(scenarios))
  
  # Run simulations for all scenarios
  for (i in seq_len(nrow(scenarios))) {
    scenario <- scenarios[i, ]
    dist_fun <- distributions[[scenario$distribution]]
    
    message(sprintf(
      "Running scenario %d/%d: n=%d, rho=%.1f, dist=%s",
      i, nrow(scenarios),
      scenario$n, scenario$rho, scenario$distribution
    ))
    
    # Special handling for distributions with additional parameters
    if (scenario$distribution == "bvt") {
      power <- simulate_power(
        n = scenario$n,
        rho = scenario$rho,
        n_sim = n_sim,
        dist_fun = function(n, rho) generate_bvt(n, rho, df = 3)
      )
    } else if (scenario$distribution == "heavy") {
      power <- simulate_power(
        n = scenario$n,
        rho = scenario$rho,
        n_sim = n_sim,
        dist_fun = function(n, rho) generate_heavy_tails(n, rho, df = 3)
      )
    } else {
      power <- simulate_power(
        n = scenario$n,
        rho = scenario$rho,
        n_sim = n_sim,
        dist_fun = dist_fun
      )
    }
    
    # Add distribution name properly
    power$distribution <- scenario$distribution
    results[[i]] <- power
  }
  
  # Combine all results
  dplyr::bind_rows(results)
}

## Section 5: Visualization and Reporting ----------------------------------

#' Create publication-quality power comparison plot
#'
#' @param power_data Data frame with power results
#' @param facet_var Variable to facet by ("n" or "rho")
#' @return ggplot object
#' @examples
#' plot_power_comparison(power_results, "n")
plot_power_comparison <- function(power_data, facet_var = "n") {
  # Clean up distribution names
  dist_labels <- c(
    bvn = "Bivariate Normal",
    bvt = "Bivariate t (df=3)",
    quadratic = "Quadratic",
    frank = "Frank Copula",
    heavy = "Heavy-tailed Margins"
  )
  
  # Handle missing distribution mappings (use original value)
  for (dist in unique(power_data$distribution)) {
    if (!(dist %in% names(dist_labels))) {
      dist_labels[dist] <- dist
    }
  }
  
  # Ensure power is numeric and complete cases
  power_data <- power_data %>%
    filter(!is.na(power)) %>%
    mutate(distribution = factor(
      distribution,
      levels = names(dist_labels),
      labels = dist_labels[names(dist_labels) %in% names(dist_labels)]
    ))
  
  # Create plot
  ggplot(power_data, aes(x = factor(rho), y = power, color = test, group = interaction(test, distribution))) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    facet_wrap(as.formula(paste("~", facet_var)), ncol = 1) +
    labs(
      x = "True Correlation (ρ)",
      y = "Power (Probability of Rejection)",
      color = "Correlation Test",
      title = "Statistical Power Comparison of Correlation Tests",
      subtitle = "Performance across different distributions and sample sizes"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "bottom",
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_line(color = "gray95"),
      strip.background = element_rect(fill = "gray90", color = NA),
      strip.text = element_text(face = "bold")
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.1),
      labels = scales::percent_format(accuracy = 1)
    ) +
    scale_color_manual(
      values = c("#1f77b4", "#ff7f0e", "#2ca02c"),  # Colorblind-friendly palette
      labels = c("Pearson", "Spearman", "Kendall")
    ) +
    guides(color = guide_legend(nrow = 1, byrow = TRUE))
}

#' Create interactive power table with formatting
#'
#' @param power_data Data frame with power results
#' @return Formatted kable table
#' @examples
#' create_power_table(power_results)
create_power_table <- function(power_data) {
  # Summarize any potential duplicates to avoid the warning
  summarized_data <- power_data %>%
    group_by(distribution, n, rho, test) %>%
    summarize(
      power = mean(power, na.rm = TRUE),
      n_failed = sum(n_failed),
      n_sim = first(n_sim),
      .groups = "drop"
    )
  
  # Create the formatted table
  table_data <- summarized_data %>%
    mutate(
      power = round(power, 3),
      n_failed = ifelse(n_failed > 0, 
                       paste0(n_sim - n_failed, "/", n_sim),
                       as.character(n_sim))
    ) %>%
    pivot_wider(
      names_from = test,
      values_from = c(power, n_failed),
      values_fn = list(power = mean, n_failed = first)
    ) %>%
    arrange(distribution, n, rho)
  
  # Calculate color range for proper spec_color
  color_values <- c(
    table_data$power_Pearson, 
    table_data$power_Spearman, 
    table_data$power_Kendall
  )
  
  # Create the kable table
  kable_table <- table_data %>%
    kable(
      caption = "Power Comparison of Correlation Tests Across Different Scenarios",
      col.names = c("Distribution", "n", "ρ", "n_sim",
                   "Pearson Power", "Spearman Power", "Kendall Power",
                   "Pearson Sims", "Spearman Sims", "Kendall Sims")
    ) %>%
    kable_styling(
      bootstrap_options = c("striped", "hover", "condensed"),
      full_width = FALSE,
      font_size = 14
    ) %>%
    column_spec(1, bold = TRUE) %>%
    column_spec(5:7, color = "white", 
               background = spec_color(color_values, option = "B", direction = -1)) %>%
    add_header_above(c(" " = 4, "Power" = 3, "Successful Simulations" = 3)) %>%
    footnote(
      general = "Power estimates based on proportion of significant tests (α = 0.05)",
      general_title = "Note:"
    )
  
  return(kable_table)
}

## Section 6: Statistical Analysis Utilities ------------------------------

#' Compare tests statistically by computing relative efficiency
#'
#' @param power_data Data frame with power results
#' @return Data frame with relative efficiency statistics
#' @examples
#' compare_test_efficiency(power_results)
compare_test_efficiency <- function(power_data) {
  # First handle potential duplicate entries
  power_data <- power_data %>%
    group_by(distribution, n, rho, test) %>%
    summarize(
      power = mean(power, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Now create the comparison
  power_data %>%
    # Use complete to ensure all test combinations exist
    tidyr::complete(distribution, n, rho, test) %>%
    group_by(distribution, n, rho) %>%
    mutate(
      pearson_power = power[test == "Pearson"]
    ) %>%
    # Filter cases where we can compute relative efficiency
    filter(test != "Pearson", !is.na(power), !is.na(pearson_power), pearson_power > 0) %>%
    mutate(relative_efficiency = power / pearson_power) %>%
    group_by(distribution, test) %>%
    summarise(
      mean_relative_efficiency = mean(relative_efficiency, na.rm = TRUE),
      median_relative_efficiency = median(relative_efficiency, na.rm = TRUE),
      min_relative_efficiency = min(relative_efficiency, na.rm = TRUE),
      max_relative_efficiency = max(relative_efficiency, na.rm = TRUE),
      sd_relative_efficiency = sd(relative_efficiency, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(distribution, desc(mean_relative_efficiency))
}

#' Estimate required sample size for target power using logistic regression
#'
#' @param power_data Data frame with power results
#' @param target_power Desired power level (default 0.8)
#' @return Data frame with estimated required sample sizes
#' @examples
#' estimate_required_n(power_results)
estimate_required_n <- function(power_data, target_power = 0.8) {
  # Handle potential duplicate entries first
  power_data <- power_data %>%
    group_by(distribution, n, rho, test) %>%
    summarize(
      power = mean(power, na.rm = TRUE),
      .groups = "drop"
    )
  
  power_data %>%
    filter(!is.na(power)) %>%
    group_by(distribution, test, rho) %>%
    group_modify(~ {
      if (nrow(.x) < 2 || sd(.x$power, na.rm = TRUE) == 0) return(tibble())
      
      model <- tryCatch({
        glm(
          power ~ n,
          data = .x,
          family = quasibinomial()
        )
      }, error = function(e) NULL)
      
      if(is.null(model)) return(tibble())
      
      # Predict over range of n values
      pred_range <- seq(min(.x$n), max(.x$n) * 2, length.out = 100)
      pred_power <- predict(model, 
                          newdata = data.frame(n = pred_range),
                          type = "response")
      
      # Find smallest n achieving target power
      sufficient_idx <- which(pred_power >= target_power)
      if(length(sufficient_idx) == 0) return(tibble())
      sufficient_n <- pred_range[sufficient_idx[1]]
      
      tibble(
        required_n = sufficient_n,
        predicted_power = pred_power[which(pred_range == sufficient_n)[1]],
        model_r2 = 1 - model$deviance/model$null.deviance
      )
    }) %>%
    filter(!is.na(required_n)) %>%
    arrange(distribution, test, rho)
}

## Section 7: Benchmarking and Performance Analysis ------------------------

#' Benchmark correlation test computation times
#'
#' @param n Sample sizes to test
#' @param rho Correlation values to test
#' @param n_iter Number of iterations per benchmark
#' @return Data frame with timing results
#' @examples
#' benchmark_correlation_tests()
benchmark_correlation_tests <- function(n = c(20, 50, 100, 200),
                                       rho = c(0.1, 0.3, 0.5),
                                       n_iter = 100) {
  scenarios <- expand.grid(n = n, rho = rho)
  
  map_dfr(seq_len(nrow(scenarios)), function(i) {
    n_val <- scenarios$n[i]
    rho_val <- scenarios$rho[i]
    
    data <- generate_bvn(n_val, rho_val)
    
    bench <- bench::mark(
      pearson = cor.test(data$x, data$y, method = "pearson"),
      spearman = cor.test(data$x, data$y, method = "spearman"),
      kendall = cor.test(data$x, data$y, method = "kendall"),
      iterations = n_iter,
      check = FALSE
    )
    
    bench %>%
      as_tibble() %>%
      mutate(
        n = n_val,
        rho = rho_val,
        method = c("pearson", "spearman", "kendall")
      ) %>%
      select(method, n, rho, median_time = median, mem_alloc)
  })
}

## Section 8: Main Analysis Execution --------------------------------------

# Create directories if they don't exist
dir.create("data/results", recursive = TRUE, showWarnings = FALSE)
dir.create("figs", recursive = TRUE, showWarnings = FALSE)

# Run the full power analysis (time-consuming step)
if (!file.exists("data/results/power_results.rds")) {
  power_results <- run_full_power_analysis(
    n_vec = c(20, 50, 100),
    rho_vec = c(0.1, 0.3, 0.5),
    n_sim = 2000
  )
  saveRDS(power_results, "data/results/power_results.rds")
} else {
  power_results <- readRDS("data/results/power_results.rds")
}

# Generate visualizations
power_plot_n <- plot_power_comparison(power_results, "n")
ggsave("figs/power_by_n.png", power_plot_n, width = 10, height = 8, dpi = 300)

power_plot_rho <- plot_power_comparison(power_results, "rho")
ggsave("figs/power_by_rho.png", power_plot_rho, width = 10, height = 8, dpi = 300)

# Create formatted tables
power_table <- create_power_table(power_results)
saveRDS(power_table, "data/results/power_table.rds")  # Save as R object first

# Use tryCatch to handle potential errors with save_kable gracefully
tryCatch({
  save_kable(power_table, "figs/power_table.html")
}, error = function(e) {
  message("Warning: Could not save HTML table. Saving plain text version instead.")
  write(kable_styling(power_table, full_width = FALSE), "figs/power_table.txt")
})

# Perform statistical comparisons with error handling
tryCatch({
  efficiency_results <- compare_test_efficiency(power_results)
  saveRDS(efficiency_results, "data/results/efficiency_results.rds")
}, error = function(e) {
  message("Warning: Error in efficiency comparison: ", e$message)
})

tryCatch({
  required_n_results <- estimate_required_n(power_results)
  saveRDS(required_n_results, "data/results/required_n_results.rds")
}, error = function(e) {
  message("Warning: Error in sample size estimation: ", e$message)
})

# Run benchmarks
tryCatch({
  timing_results <- benchmark_correlation_tests(n_iter = 200)
  saveRDS(timing_results, "data/results/timing_results.rds")
}, error = function(e) {
  message("Warning: Error in benchmarking: ", e$message)
})

# Print a summary of the results
cat("Power analysis completed successfully.\n")
cat("Results saved in data/results/ directory.\n")
cat("Plots saved in figs/ directory.\n")
