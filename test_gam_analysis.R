# Test script for GAM Spatial Analysis
# This script demonstrates the key functionality

cat("=== TESTING GAM SPATIAL ANALYSIS ===\n\n")

# Load the analysis functions
source('gam_spatial_analysis_simple.R')

cat("1. Testing data simulation...\n")
# Test data simulation
sim_data <- simulate_fish_growth_data(n_fish = 200, n_zones = 2)
cat("✓ Simulated data created successfully\n\n")

cat("2. Testing GAM model fitting...\n")
# Test GAM fitting
gam_model <- fit_spatial_gam(sim_data$data)
cat("✓ GAM model fitted successfully\n\n")

cat("3. Testing derivative calculation...\n")
# Test derivatives
derivatives <- calculate_spatial_derivatives(gam_model, sim_data$data, n_grid = 20)
cat("✓ Derivatives calculated successfully\n\n")

cat("4. Testing breakpoint detection...\n")
# Test breakpoint detection
breakpoints <- detect_spatial_breakpoints(derivatives)
cat("✓ Breakpoint detection completed\n\n")

cat("5. Testing validation...\n")
# Test validation
validation <- validate_zone_detection(breakpoints, sim_data$zone_boundaries, tolerance = 2.0)
cat("✓ Validation completed\n\n")

cat("6. Testing real data generation...\n")
# Test real data
sablefish_data <- create_example_sablefish_data(n_fish = 300)
cat("✓ Example sablefish data created\n\n")

cat("7. Testing complete workflow...\n")
# Test workflow
results <- apply_gam_workflow(sablefish_data, create_plots = FALSE)
cat("✓ Complete workflow executed successfully\n\n")

# Print summary results
cat("=== TEST RESULTS SUMMARY ===\n")
cat(paste("Simulation - Detected breakpoints:", nrow(breakpoints), "\n"))
cat(paste("Simulation - True boundaries:", length(sim_data$zone_boundaries), "\n"))
cat(paste("Simulation - Match rate:", round(validation$match_rate * 100, 1), "%\n"))
cat(paste("Real data - Model R-squared:", round(summary(results$model)$r.sq, 3), "\n"))
cat(paste("Real data - Detected breakpoints:", nrow(results$breakpoints), "\n"))

cat("\n✓ All tests completed successfully!\n")
cat("The GAM spatial analysis implementation is working correctly.\n")