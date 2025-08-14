# GAM Spatial Analysis for Pacific Cod Project
# Implementation of Generalized Additive Models for spatial variation in fish growth
# Author: AFSC-NOAA Pacific Cod Project
# Date: 2024

# ============================================================================
# PACKAGE INSTALLATION AND LOADING
# ============================================================================

# Load required packages (install separately if needed)
required_packages <- c("mgcv", "ggplot2", "dplyr", "maps", "mapdata", 
                      "viridis", "gridExtra", "sp", "raster", "fields")

# Function to safely load packages
load_packages <- function(packages) {
  loaded <- c()
  missing <- c()
  
  for (pkg in packages) {
    if (require(pkg, character.only = TRUE, quietly = TRUE)) {
      loaded <- c(loaded, pkg)
    } else {
      missing <- c(missing, pkg)
    }
  }
  
  if (length(loaded) > 0) {
    cat("Loaded packages:", paste(loaded, collapse = ", "), "\n")
  }
  
  if (length(missing) > 0) {
    cat("Missing packages (install with install.packages()):", paste(missing, collapse = ", "), "\n")
    cat("Note: Some functions may not work without these packages.\n")
  }
  
  return(list(loaded = loaded, missing = missing))
}

cat("Loading required packages...\n")
package_status <- load_packages(required_packages)

# ============================================================================
# CORE GAM FUNCTIONS
# ============================================================================

#' Fit GAM model for length-age-spatial relationships
#' 
#' @param data Data frame with columns: length, age, latitude, longitude
#' @param k_age Basis dimension for age smooth (default: 10)
#' @param k_spatial Basis dimension for spatial smooth (default: 20)
#' @return Fitted GAM model object
fit_spatial_gam <- function(data, k_age = 10, k_spatial = 20) {
  cat("Fitting GAM model: length ~ s(age) + s(latitude, longitude)...\n")
  
  # Fit GAM with age and spatial smooths
  gam_model <- mgcv::gam(
    length ~ s(age, k = k_age) + s(latitude, longitude, k = k_spatial),
    data = data,
    method = "REML"
  )
  
  # Print model summary
  cat("\nGAM Model Summary:\n")
  print(summary(gam_model))
  
  return(gam_model)
}

#' Calculate first derivatives of spatial smooth term
#' 
#' @param gam_model Fitted GAM model
#' @param data Original data used to fit model
#' @param n_grid Number of grid points for derivative calculation (default: 50)
#' @return List containing derivative information
calculate_spatial_derivatives <- function(gam_model, data, n_grid = 50) {
  cat("Calculating first derivatives of spatial smooth term...\n")
  
  # Create prediction grid
  lat_range <- range(data$latitude, na.rm = TRUE)
  lon_range <- range(data$longitude, na.rm = TRUE)
  
  lat_seq <- seq(lat_range[1], lat_range[2], length.out = n_grid)
  lon_seq <- seq(lon_range[1], lon_range[2], length.out = n_grid)
  
  pred_grid <- expand.grid(
    latitude = lat_seq,
    longitude = lon_seq,
    age = median(data$age, na.rm = TRUE)  # Fix age at median
  )
  
  # Calculate derivatives using finite differences
  delta_lat <- lat_seq[2] - lat_seq[1]
  delta_lon <- lon_seq[2] - lon_seq[1]
  
  # Predict on grid
  predictions <- predict(gam_model, newdata = pred_grid, type = "response")
  pred_matrix <- matrix(predictions, nrow = n_grid, ncol = n_grid)
  
  # Calculate gradients
  grad_lat <- array(NA, dim = c(n_grid, n_grid))
  grad_lon <- array(NA, dim = c(n_grid, n_grid))
  
  # Central differences for interior points
  for (i in 2:(n_grid-1)) {
    for (j in 2:(n_grid-1)) {
      grad_lat[i, j] <- (pred_matrix[i+1, j] - pred_matrix[i-1, j]) / (2 * delta_lat)
      grad_lon[i, j] <- (pred_matrix[i, j+1] - pred_matrix[i, j-1]) / (2 * delta_lon)
    }
  }
  
  # Calculate gradient magnitude
  grad_magnitude <- sqrt(grad_lat^2 + grad_lon^2)
  
  return(list(
    latitude = lat_seq,
    longitude = lon_seq,
    pred_matrix = pred_matrix,
    grad_lat = grad_lat,
    grad_lon = grad_lon,
    grad_magnitude = grad_magnitude,
    pred_grid = pred_grid
  ))
}

#' Detect spatial breakpoints based on derivative analysis
#' 
#' @param derivatives Output from calculate_spatial_derivatives
#' @param threshold Threshold for gradient magnitude (default: quantile 0.8)
#' @return Data frame of detected breakpoints
detect_spatial_breakpoints <- function(derivatives, threshold = NULL) {
  cat("Detecting spatial breakpoints...\n")
  
  # Set threshold if not provided
  if (is.null(threshold)) {
    threshold <- quantile(derivatives$grad_magnitude, 0.8, na.rm = TRUE)
  }
  
  cat(paste("Using gradient magnitude threshold:", round(threshold, 4), "\n"))
  
  # Find points exceeding threshold
  high_grad_points <- which(derivatives$grad_magnitude > threshold, arr.ind = TRUE)
  
  if (nrow(high_grad_points) == 0) {
    cat("No breakpoints detected with current threshold.\n")
    return(data.frame())
  }
  
  # Convert indices to coordinates
  breakpoints <- data.frame(
    latitude = derivatives$latitude[high_grad_points[, 1]],
    longitude = derivatives$longitude[high_grad_points[, 2]],
    grad_magnitude = derivatives$grad_magnitude[high_grad_points]
  )
  
  cat(paste("Detected", nrow(breakpoints), "potential breakpoints.\n"))
  return(breakpoints)
}

# ============================================================================
# SIMULATION FUNCTIONS
# ============================================================================

#' Simulate spatially structured fish growth data
#' 
#' @param n_fish Number of fish to simulate (default: 1000)
#' @param lat_range Latitude range (default: c(30, 60))
#' @param lon_range Longitude range (default: c(-180, -120))
#' @param n_zones Number of spatial zones (default: 3)
#' @param age_range Age range (default: c(1, 15))
#' @return List containing simulated data and true zone boundaries
simulate_fish_growth_data <- function(n_fish = 1000, lat_range = c(30, 60), 
                                     lon_range = c(-180, -120), n_zones = 3,
                                     age_range = c(1, 15)) {
  cat("Simulating spatially structured fish growth data...\n")
  
  # Create spatial zones
  lat_breaks <- seq(lat_range[1], lat_range[2], length.out = n_zones + 1)
  zone_boundaries <- lat_breaks[2:n_zones]  # Interior boundaries
  
  # Simulate fish locations
  latitude <- runif(n_fish, lat_range[1], lat_range[2])
  longitude <- runif(n_fish, lon_range[1], lon_range[2])
  age <- runif(n_fish, age_range[1], age_range[2])
  
  # Assign zones based on latitude
  zone <- cut(latitude, breaks = lat_breaks, labels = FALSE, include.lowest = TRUE)
  
  # Zone-specific growth parameters
  zone_effects <- c(-0.5, 0, 0.5)[1:n_zones]  # Different growth rates by zone
  length_base <- 20 + 3 * age  # Base length-age relationship
  
  # Add zone effects and random noise
  length <- length_base + zone_effects[zone] * age + rnorm(n_fish, 0, 2)
  
  # Create data frame
  sim_data <- data.frame(
    fish_id = 1:n_fish,
    latitude = latitude,
    longitude = longitude,
    age = age,
    length = length,
    true_zone = zone
  )
  
  cat(paste("Simulated", n_fish, "fish in", n_zones, "zones.\n"))
  cat("Zone boundaries (latitude):", paste(round(zone_boundaries, 2), collapse = ", "), "\n")
  
  return(list(
    data = sim_data,
    zone_boundaries = zone_boundaries,
    zone_effects = zone_effects
  ))
}

#' Validate detected zones against known boundaries
#' 
#' @param detected_breakpoints Output from detect_spatial_breakpoints
#' @param true_boundaries True zone boundaries from simulation
#' @param tolerance Tolerance for matching (degrees latitude)
#' @return Validation results
validate_zone_detection <- function(detected_breakpoints, true_boundaries, tolerance = 1.0) {
  cat("Validating detected zones against known boundaries...\n")
  
  if (nrow(detected_breakpoints) == 0) {
    return(list(
      n_detected = 0,
      n_true = length(true_boundaries),
      n_matched = 0,
      match_rate = 0,
      false_positives = 0
    ))
  }
  
  # Find detected breakpoints near true boundaries
  detected_lats <- detected_breakpoints$latitude
  n_matched <- 0
  matched_boundaries <- c()
  
  for (true_lat in true_boundaries) {
    distances <- abs(detected_lats - true_lat)
    min_distance <- min(distances)
    
    if (min_distance <= tolerance) {
      n_matched <- n_matched + 1
      matched_boundaries <- c(matched_boundaries, true_lat)
      cat(paste("Matched boundary at", round(true_lat, 2), 
                "with detected point at", round(detected_lats[which.min(distances)], 2), "\n"))
    }
  }
  
  # Calculate validation metrics
  n_detected <- length(unique(round(detected_lats, 1)))
  n_true <- length(true_boundaries)
  false_positives <- max(0, n_detected - n_matched)
  match_rate <- n_matched / n_true
  
  cat(paste("Validation Results:\n"))
  cat(paste("  True boundaries:", n_true, "\n"))
  cat(paste("  Detected zones:", n_detected, "\n"))
  cat(paste("  Matched boundaries:", n_matched, "\n"))
  cat(paste("  Match rate:", round(match_rate * 100, 1), "%\n"))
  cat(paste("  False positives:", false_positives, "\n"))
  
  return(list(
    n_detected = n_detected,
    n_true = n_true,
    n_matched = n_matched,
    match_rate = match_rate,
    false_positives = false_positives,
    matched_boundaries = matched_boundaries
  ))
}

# ============================================================================
# VISUALIZATION FUNCTIONS
# ============================================================================

#' Create comprehensive visualization of spatial analysis
#' 
#' @param data Fish data
#' @param gam_model Fitted GAM model
#' @param derivatives Derivative calculations
#' @param breakpoints Detected breakpoints
#' @return List of ggplot objects
create_spatial_plots <- function(data, gam_model, derivatives, breakpoints) {
  cat("Creating spatial visualization plots...\n")
  
  # Plot 1: Raw data distribution
  p1 <- ggplot(data, aes(x = longitude, y = latitude, color = length)) +
    geom_point(alpha = 0.6) +
    scale_color_viridis_c(name = "Length") +
    labs(title = "Fish Length Distribution",
         x = "Longitude", y = "Latitude") +
    theme_minimal()
  
  # Plot 2: GAM predictions
  pred_df <- data.frame(
    longitude = rep(derivatives$longitude, each = length(derivatives$latitude)),
    latitude = rep(derivatives$latitude, times = length(derivatives$longitude)),
    prediction = as.vector(derivatives$pred_matrix)
  )
  
  p2 <- ggplot(pred_df, aes(x = longitude, y = latitude, fill = prediction)) +
    geom_raster() +
    scale_fill_viridis_c(name = "Predicted\nLength") +
    labs(title = "GAM Spatial Predictions",
         x = "Longitude", y = "Latitude") +
    theme_minimal()
  
  # Plot 3: Gradient magnitude
  grad_df <- data.frame(
    longitude = rep(derivatives$longitude, each = length(derivatives$latitude)),
    latitude = rep(derivatives$latitude, times = length(derivatives$longitude)),
    gradient = as.vector(derivatives$grad_magnitude)
  )
  
  p3 <- ggplot(grad_df, aes(x = longitude, y = latitude, fill = gradient)) +
    geom_raster() +
    scale_fill_viridis_c(name = "Gradient\nMagnitude", option = "plasma") +
    labs(title = "Spatial Gradient Magnitude",
         x = "Longitude", y = "Latitude") +
    theme_minimal()
  
  # Plot 4: Detected breakpoints
  p4 <- p3
  if (nrow(breakpoints) > 0) {
    p4 <- p4 + 
      geom_point(data = breakpoints, aes(x = longitude, y = latitude), 
                 color = "red", size = 2, inherit.aes = FALSE) +
      labs(title = "Detected Breakpoints")
  }
  
  return(list(data_plot = p1, predictions = p2, gradients = p3, breakpoints = p4))
}

#' Create summary visualization combining all plots
#' 
#' @param plots List of plots from create_spatial_plots
#' @return Combined plot
create_summary_plot <- function(plots) {
  gridExtra::grid.arrange(
    plots$data_plot, plots$predictions,
    plots$gradients, plots$breakpoints,
    ncol = 2, nrow = 2
  )
}

# ============================================================================
# REAL DATA APPLICATION FUNCTIONS
# ============================================================================

#' Create example Northeast Pacific sablefish dataset
#' 
#' This function creates a realistic example dataset based on known
#' sablefish distribution patterns in the Northeast Pacific
create_example_sablefish_data <- function(n_fish = 2000) {
  cat("Creating example Northeast Pacific sablefish dataset...\n")
  
  # Define realistic spatial boundaries for NE Pacific
  # Southern California Bight: ~32-35°N
  # North Pacific Current transition: ~40-42°N
  # Alaska Gyre boundary: ~50-52°N
  
  set.seed(123)  # For reproducibility
  
  # Create spatially clustered sampling
  lat_centers <- c(33.5, 41, 51)  # Major oceanographic regions
  lon_centers <- c(-120, -125, -140)
  
  # Simulate fish locations with clustering
  region <- sample(1:3, n_fish, replace = TRUE, prob = c(0.3, 0.4, 0.3))
  
  latitude <- numeric(n_fish)
  longitude <- numeric(n_fish)
  
  for (i in 1:n_fish) {
    r <- region[i]
    latitude[i] <- rnorm(1, lat_centers[r], 2)
    longitude[i] <- rnorm(1, lon_centers[r], 3)
  }
  
  # Constrain to realistic ranges
  latitude <- pmax(30, pmin(60, latitude))
  longitude <- pmax(-160, pmin(-110, longitude))
  
  # Simulate ages (sablefish can live 90+ years, but survey typically catches 2-20)
  age <- rpois(n_fish, lambda = 8) + 1
  age <- pmax(1, pmin(25, age))
  
  # Simulate lengths with regional effects
  # Southern region: warmer water, faster early growth
  # Central region: moderate growth
  # Northern region: slower growth, larger maximum size
  
  base_length <- 15 + 2.5 * age  # Base von Bertalanffy-like growth
  
  regional_effects <- c(0.8, 0, -0.3)[region]  # Regional growth modifiers
  temp_effect <- -0.1 * (latitude - 45)  # Temperature gradient effect
  
  length <- base_length + regional_effects * sqrt(age) + 
           temp_effect * age + rnorm(n_fish, 0, 3)
  
  # Ensure positive lengths
  length <- pmax(5, length)
  
  # Create data frame
  sablefish_data <- data.frame(
    fish_id = 1:n_fish,
    latitude = latitude,
    longitude = longitude,
    age = age,
    length = length,
    region = region,
    depth = rnorm(n_fish, 300, 100),  # Typical sablefish depths
    year = sample(2010:2023, n_fish, replace = TRUE)
  )
  
  cat(paste("Created dataset with", n_fish, "sablefish observations.\n"))
  cat("Latitude range:", round(range(latitude), 1), "\n")
  cat("Longitude range:", round(range(longitude), 1), "\n")
  cat("Age range:", range(age), "\n")
  cat("Length range:", round(range(length), 1), "\n")
  
  return(sablefish_data)
}

#' Apply complete GAM analysis workflow to real data
#' 
#' @param data Fish data (must have columns: length, age, latitude, longitude)
#' @param plot_results Whether to create plots (default: TRUE)
#' @return Complete analysis results
apply_gam_workflow <- function(data, plot_results = TRUE) {
  cat("\n=== APPLYING COMPLETE GAM SPATIAL ANALYSIS WORKFLOW ===\n")
  
  # Step 1: Fit GAM model
  gam_model <- fit_spatial_gam(data)
  
  # Step 2: Calculate derivatives
  derivatives <- calculate_spatial_derivatives(gam_model, data)
  
  # Step 3: Detect breakpoints
  breakpoints <- detect_spatial_breakpoints(derivatives)
  
  # Step 4: Create visualizations
  if (plot_results && nrow(data) > 0) {
    plots <- create_spatial_plots(data, gam_model, derivatives, breakpoints)
    
    # Display summary plot
    cat("Creating summary visualization...\n")
    create_summary_plot(plots)
    
    # Save individual plots
    ggsave("spatial_data_distribution.png", plots$data_plot, width = 10, height = 8)
    ggsave("gam_predictions.png", plots$predictions, width = 10, height = 8)
    ggsave("spatial_gradients.png", plots$gradients, width = 10, height = 8)
    ggsave("detected_breakpoints.png", plots$breakpoints, width = 10, height = 8)
    
    cat("Individual plots saved as PNG files.\n")
  }
  
  # Step 5: Create results summary
  results <- list(
    model = gam_model,
    derivatives = derivatives,
    breakpoints = breakpoints,
    data = data
  )
  
  if (plot_results) {
    results$plots <- plots
  }
  
  cat("\n=== ANALYSIS COMPLETE ===\n")
  return(results)
}

# ============================================================================
# DEMONSTRATION AND TESTING
# ============================================================================

#' Run complete demonstration of GAM spatial analysis
run_demo <- function() {
  cat("\n=== GAM SPATIAL ANALYSIS DEMONSTRATION ===\n")
  
  # Part 1: Simulation Testing
  cat("\n--- PART 1: SIMULATION TESTING ---\n")
  
  # Generate simulated data
  sim_results <- simulate_fish_growth_data(n_fish = 1500, n_zones = 3)
  sim_data <- sim_results$data
  
  # Apply GAM workflow to simulated data
  sim_analysis <- apply_gam_workflow(sim_data, plot_results = FALSE)
  
  # Validate results
  validation <- validate_zone_detection(
    sim_analysis$breakpoints, 
    sim_results$zone_boundaries,
    tolerance = 1.5
  )
  
  # Part 2: Real Data Application
  cat("\n--- PART 2: REAL DATA APPLICATION ---\n")
  
  # Create example sablefish data
  sablefish_data <- create_example_sablefish_data(n_fish = 2000)
  
  # Apply GAM workflow
  real_analysis <- apply_gam_workflow(sablefish_data, plot_results = TRUE)
  
  # Print detected breakpoints with oceanographic context
  if (nrow(real_analysis$breakpoints) > 0) {
    cat("\nDetected breakpoints in Northeast Pacific:\n")
    for (i in 1:nrow(real_analysis$breakpoints)) {
      lat <- real_analysis$breakpoints$latitude[i]
      lon <- real_analysis$breakpoints$longitude[i]
      
      # Classify by known oceanographic features
      if (lat >= 32 && lat <= 35) {
        region <- "Southern California Bight region"
      } else if (lat >= 40 && lat <= 42) {
        region <- "North Pacific Current transition"
      } else if (lat >= 50 && lat <= 52) {
        region <- "Alaska Gyre boundary region"
      } else {
        region <- "Other region"
      }
      
      cat(sprintf("  Breakpoint %d: %.2f°N, %.2f°W (%s)\n", 
                  i, lat, abs(lon), region))
    }
  }
  
  # Return combined results
  return(list(
    simulation = list(
      data = sim_data,
      analysis = sim_analysis,
      validation = validation,
      true_boundaries = sim_results$zone_boundaries
    ),
    real_data = list(
      data = sablefish_data,
      analysis = real_analysis
    )
  ))
}

# ============================================================================
# MAIN EXECUTION
# ============================================================================

cat("GAM Spatial Analysis Package Loaded Successfully!\n")
cat("Key functions available:\n")
cat("  - fit_spatial_gam(): Fit GAM model\n")
cat("  - calculate_spatial_derivatives(): Calculate spatial gradients\n")
cat("  - detect_spatial_breakpoints(): Find spatial breakpoints\n")
cat("  - simulate_fish_growth_data(): Create test data\n")
cat("  - apply_gam_workflow(): Complete analysis pipeline\n")
cat("  - run_demo(): Full demonstration\n")
cat("\nTo run the complete demonstration, use: run_demo()\n")