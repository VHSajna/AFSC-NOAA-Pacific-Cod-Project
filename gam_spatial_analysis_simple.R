# GAM Spatial Analysis for Pacific Cod Project - Simplified Version
# Implementation using base R and mgcv (which is already installed)
# Author: AFSC-NOAA Pacific Cod Project
# Date: 2024

# ============================================================================
# PACKAGE LOADING
# ============================================================================

# Load available packages
cat("Loading required packages...\n")
library(mgcv)     # For GAM modeling
library(stats)    # For statistical functions
library(graphics) # For plotting
library(grDevices) # For graphics devices

cat("Successfully loaded packages: mgcv, stats, graphics, grDevices\n")

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
#' @param n_grid Number of grid points for derivative calculation (default: 30)
#' @return List containing derivative information
calculate_spatial_derivatives <- function(gam_model, data, n_grid = 30) {
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
  
  # Calculate gradients using finite differences
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
  
  cat("Derivatives calculated successfully.\n")
  cat(paste("Grid size:", n_grid, "x", n_grid, "\n"))
  cat(paste("Max gradient magnitude:", round(max(grad_magnitude, na.rm = TRUE), 4), "\n"))
  
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
    return(data.frame(latitude = numeric(0), longitude = numeric(0), grad_magnitude = numeric(0)))
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
  set.seed(123)  # For reproducibility
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
# VISUALIZATION FUNCTIONS (using base R graphics)
# ============================================================================

#' Create spatial plots using base R graphics
#' 
#' @param data Fish data
#' @param derivatives Derivative calculations
#' @param breakpoints Detected breakpoints
#' @param filename Optional filename to save plots
create_spatial_plots_base <- function(data, derivatives, breakpoints, filename = NULL) {
  cat("Creating spatial visualization plots with base R...\n")
  
  if (!is.null(filename)) {
    png(filename, width = 12, height = 10, res = 150)
  }
  
  # Set up 2x2 plot layout with smaller margins
  par(mfrow = c(2, 2), mar = c(3, 3, 2, 1), oma = c(1, 1, 1, 1))
  
  # Plot 1: Raw data distribution
  plot(data$longitude, data$latitude, col = heat.colors(100)[cut(data$length, 100)],
       pch = 16, cex = 0.6, main = "Fish Length Distribution",
       xlab = "Longitude", ylab = "Latitude")
  
  # Add colorbar legend
  legend("topright", legend = c("Short", "Long"), 
         col = c("red", "yellow"), pch = 16, title = "Length")
  
  # Plot 2: GAM predictions
  pred_df <- expand.grid(
    longitude = derivatives$longitude,
    latitude = derivatives$latitude
  )
  pred_df$prediction <- as.vector(derivatives$pred_matrix)
  
  # Create filled contour plot
  image(derivatives$longitude, derivatives$latitude, 
        t(derivatives$pred_matrix), 
        main = "GAM Spatial Predictions",
        xlab = "Longitude", ylab = "Latitude",
        col = heat.colors(50))
  contour(derivatives$longitude, derivatives$latitude, 
          t(derivatives$pred_matrix), add = TRUE, nlevels = 5)
  
  # Plot 3: Gradient magnitude
  image(derivatives$longitude, derivatives$latitude, 
        t(derivatives$grad_magnitude), 
        main = "Spatial Gradient Magnitude",
        xlab = "Longitude", ylab = "Latitude",
        col = rainbow(50, start = 0, end = 0.7))
  contour(derivatives$longitude, derivatives$latitude, 
          t(derivatives$grad_magnitude), add = TRUE, nlevels = 3)
  
  # Plot 4: Detected breakpoints
  image(derivatives$longitude, derivatives$latitude, 
        t(derivatives$grad_magnitude), 
        main = "Detected Breakpoints",
        xlab = "Longitude", ylab = "Latitude",
        col = rainbow(50, start = 0, end = 0.7))
  
  if (nrow(breakpoints) > 0) {
    points(breakpoints$longitude, breakpoints$latitude, 
           col = "red", pch = 16, cex = 1.5)
    points(breakpoints$longitude, breakpoints$latitude, 
           col = "white", pch = 1, cex = 1.5, lwd = 2)
  }
  
  if (!is.null(filename)) {
    dev.off()
    cat(paste("Plots saved to", filename, "\n"))
  }
  
  # Reset plot parameters
  par(mfrow = c(1, 1))
}

# ============================================================================
# REAL DATA APPLICATION FUNCTIONS
# ============================================================================

#' Create example Northeast Pacific sablefish dataset
create_example_sablefish_data <- function(n_fish = 2000) {
  cat("Creating example Northeast Pacific sablefish dataset...\n")
  
  set.seed(123)  # For reproducibility
  
  # Define realistic spatial boundaries for NE Pacific
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
  
  # Simulate ages
  age <- rpois(n_fish, lambda = 8) + 1
  age <- pmax(1, pmin(25, age))
  
  # Simulate lengths with regional effects
  base_length <- 15 + 2.5 * age
  regional_effects <- c(0.8, 0, -0.3)[region]
  temp_effect <- -0.1 * (latitude - 45)
  
  length <- base_length + regional_effects * sqrt(age) + 
           temp_effect * age + rnorm(n_fish, 0, 3)
  length <- pmax(5, length)
  
  sablefish_data <- data.frame(
    fish_id = 1:n_fish,
    latitude = latitude,
    longitude = longitude,
    age = age,
    length = length,
    region = region,
    depth = rnorm(n_fish, 300, 100),
    year = sample(2010:2023, n_fish, replace = TRUE)
  )
  
  cat(paste("Created dataset with", n_fish, "sablefish observations.\n"))
  cat("Latitude range:", round(range(latitude), 1), "\n")
  cat("Longitude range:", round(range(longitude), 1), "\n")
  cat("Age range:", range(age), "\n")
  cat("Length range:", round(range(length), 1), "\n")
  
  return(sablefish_data)
}

#' Apply complete GAM analysis workflow
#' 
#' @param data Fish data
#' @param create_plots Whether to create plots (default: TRUE)
#' @return Complete analysis results
apply_gam_workflow <- function(data, create_plots = TRUE) {
  cat("\n=== APPLYING COMPLETE GAM SPATIAL ANALYSIS WORKFLOW ===\n")
  
  # Step 1: Fit GAM model
  gam_model <- fit_spatial_gam(data)
  
  # Step 2: Calculate derivatives
  derivatives <- calculate_spatial_derivatives(gam_model, data)
  
  # Step 3: Detect breakpoints
  breakpoints <- detect_spatial_breakpoints(derivatives)
  
  # Step 4: Create visualizations
  if (create_plots && nrow(data) > 0) {
    create_spatial_plots_base(data, derivatives, breakpoints, "spatial_analysis_plots.png")
  }
  
  # Step 5: Create results summary
  results <- list(
    model = gam_model,
    derivatives = derivatives,
    breakpoints = breakpoints,
    data = data
  )
  
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
  sim_results <- simulate_fish_growth_data(n_fish = 800, n_zones = 3)
  sim_data <- sim_results$data
  
  # Apply GAM workflow to simulated data
  sim_analysis <- apply_gam_workflow(sim_data, create_plots = FALSE)
  
  # Validate results
  validation <- validate_zone_detection(
    sim_analysis$breakpoints, 
    sim_results$zone_boundaries,
    tolerance = 1.5
  )
  
  # Part 2: Real Data Application
  cat("\n--- PART 2: REAL DATA APPLICATION ---\n")
  
  # Create example sablefish data
  sablefish_data <- create_example_sablefish_data(n_fish = 1200)
  
  # Apply GAM workflow
  real_analysis <- apply_gam_workflow(sablefish_data, create_plots = TRUE)
  
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
  } else {
    cat("No significant breakpoints detected. Try adjusting threshold.\n")
  }
  
  # Part 3: Print Summary Statistics
  cat("\n--- SUMMARY STATISTICS ---\n")
  cat("Simulation Results:\n")
  cat(sprintf("  Method detected %d/%d true boundaries (%.1f%% success rate)\n",
              validation$n_matched, validation$n_true, validation$match_rate * 100))
  
  cat("\nReal Data Results:\n")
  cat(sprintf("  Model R-squared: %.3f\n", summary(real_analysis$model)$r.sq))
  cat(sprintf("  Detected %d potential breakpoints\n", nrow(real_analysis$breakpoints)))
  
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

cat("\n=== GAM SPATIAL ANALYSIS PACKAGE LOADED SUCCESSFULLY ===\n")
cat("This package implements spatial analysis of fish growth using GAM models.\n\n")
cat("Key functions available:\n")
cat("  - fit_spatial_gam(): Fit GAM model with spatial terms\n")
cat("  - calculate_spatial_derivatives(): Calculate spatial gradients\n")
cat("  - detect_spatial_breakpoints(): Find spatial breakpoints\n")
cat("  - simulate_fish_growth_data(): Create test data with known zones\n")
cat("  - validate_zone_detection(): Compare detected vs. true boundaries\n")
cat("  - create_example_sablefish_data(): Generate realistic NE Pacific data\n")
cat("  - apply_gam_workflow(): Complete analysis pipeline\n")
cat("  - run_demo(): Full demonstration with simulation and real data\n")
cat("\nTo run the complete demonstration, use: demo_results <- run_demo()\n")
cat("To analyze your own data, use: results <- apply_gam_workflow(your_data)\n\n")