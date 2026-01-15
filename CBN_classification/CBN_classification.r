library(clustNet)

# Read the binary matrix, skipping the first line and removing first two columns
binary_matrix <- as.matrix(read.table("../dataset/sample_matrices/sample_matrix_Pan.txt", header = FALSE, skip = 1)[, -c(1,2)])

args <- commandArgs(trailingOnly = TRUE)
k <- as.integer(args[1])

sprintf("%.0f clusters", k)

result <- bestAICsearch(binary_matrix, minK = k, maxK = k, chiVec = 0, startseed = 42, nIterations = 20, AICrange = 100, plot_heatmap = TRUE)

# View the results
sprintf("AIC for %.0f clusters: %f", k, result$aics[1,1])

write.table(t(result$output[[1]]$relativeweights), "sample_Ps_CBN.dat", sep=" ", row.names=FALSE, col.names=FALSE)
