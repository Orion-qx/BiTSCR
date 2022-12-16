#' BiTSC_v1
#'
#' First version of BiTSC; developed Dec 10, 2022
#'
#' @param rho_0 the porportion to be sampled
#' @param niter_0 number of iterations
#'
#' @return the same x
#'
#' @import readr dplyr matrixLaplacian wordspace matrixStats rARPACK Rfast Matrix
#'
#' @export
#'

bitsc <- function(rho_0, niter_0, urlfile1, urlfile2, urlfile3) {

  # urlfile1 = "https://raw.githubusercontent.com/edensunyidan/BiTSC/master/data/node_covariate_one.csv"
  # urlfile2 = "https://raw.githubusercontent.com/edensunyidan/BiTSC/master/data/node_covariate_two.csv"
  # urlfile3 = "https://raw.githubusercontent.com/edensunyidan/BiTSC/master/data/edge_one_two.csv"

  X1.orig = read_csv(url(urlfile1))
  X2.orig = read_csv(url(urlfile2))
  A.orig = read_csv(url(urlfile3))

  X1 = as.matrix(scale(log2(X1.orig[,-1] + 1)))
  X2 = as.matrix(scale(log2(X2.orig[,-1] + 1)))

  gene1 = X1.orig$gene_id
  gene2 = X2.orig$gene_id

  num_m = nrow(X1)
  num_n = nrow(X2)
  p1 = ncol(X1)
  p2 = ncol(X2)

  A <- matrix(0, num_m, num_n, dimnames = list(gene1, gene2))
  A = sparseMatrix(i = match(A.orig$S1, rownames(A)),
                   j = match(A.orig$S2, colnames(A)),
                   x = A.orig$Value)

  # euclid1_xij = as.matrix(Dist(X1, method = "euclidean", diag = TRUE, upper = TRUE))
  euclid1_xij = as.matrix(Dist(X1, method = "euclidean"))
  mat_K1 = exp(-(euclid1_xij)^2/p1)

  # euclid2_xij = as.matrix(Dist(X2, method = "euclidean", diag = TRUE, upper = TRUE))
  euclid2_xij = as.matrix(Dist(X2, method = "euclidean"))
  mat_K2 = exp(-(euclid2_xij)^2/p2)

  # set Tau value:
  Tau1 = 1; Tau2 = 1

  # start_time <- Sys.time()
  B = mat_K1 %*% A %*% mat_K2 + Tau1*A %*% mat_K2 + Tau2*mat_K1 %*% A + Tau1*Tau2*A
  # browser()
  #
  B <- as.matrix(B)
  # print(typeof(B))
  B_T = t(B)

  B0 = matrix(0, nrow = num_m, ncol = num_m)
  B2 = B
  B1 = B_T
  B3 = matrix(0, nrow = num_n, ncol = num_n)
  W1 = cbind(B0,B2)
  W2 = cbind(B1, B3)
  W = rbind(W1,W2)

  diagnal = c(rowSums(B),colSums(B))

  # start_time <- Sys.time()
  D_sqrt_inv = Diagonal(num_m + num_n, 1/sqrt(diagnal))
  L1 = diag(num_m + num_n) - D_sqrt_inv %*% W %*% D_sqrt_inv

  sample_func <- function(niter, rho, K0 = 15) {
    avg.M = 0
    for (h in 1:niter) {

      m_sim = floor(rho*num_m)
      n_sim = floor(rho*num_n)

      gene1_sim = sample(gene1, size = m_sim)
      gene2_sim = sample(gene2, size = n_sim)

      match_gene1_idx <- match(gene1_sim, gene1)
      match_gene2_idx <- match(gene2_sim, gene2)

      X1_sim = X1[match_gene1_idx,]
      X2_sim = X2[match_gene2_idx,]

      w_col_select <- match_gene2_idx + num_m # first m cols are 0
      w_row_select <- match_gene1_idx + num_n # first n rows are 0

      L1_sim <- L1[c(match_gene2_idx, w_row_select), c(match_gene1_idx,w_col_select)]
      eigen_L1_sim <- eigs_sym(L1_sim, k = K0)
      eigen_vec_sim = Re(eigen_L1_sim$vectors)

      U = eigen_vec_sim
      cluster_sim = kmeans(U, K0)

      # Reassign the unsampled part
      gene1_unsampled <- gene1[-match_gene1_idx] # name of unsampled gene
      gene2_unsampled <- gene2[-match_gene2_idx]
      gene1_uns_idx <- (1:num_m)[-match_gene1_idx] # index of unsampled gene
      gene2_uns_idx <- (1:num_n)[-match_gene2_idx]

      MCV1 <- matrix(nrow = p1, ncol = K0)
      MCV2 <- matrix(nrow = p2, ncol = K0)
      # long_gene_idx <- c(match_gene2_idx, match_gene1_idx)
      for (i in 1:K0) {
        # split clusters into 2 genes:

        cur_clusters <- which(cluster_sim$cluster == i) # the index of nodes in cluster that belongs to this cluster i
        is_gene2 <- cur_clusters <= n_sim # if the index is smaller than
        cur_gene2_cluster <- cur_clusters[is_gene2]
        cur_gene1_cluster <- cur_clusters[!is_gene2]

        if (length(cur_gene2_cluster) >1 ) {
          MCV2[,i] <- colMeans(X2[match_gene2_idx[cur_gene2_cluster],])
        } else if (length(cur_gene2_cluster) == 1) {
          MCV2[,i] <- X2[match_gene2_idx[cur_gene2_cluster],]
        }

        if (length(cur_gene1_cluster) > 1 ) {
          MCV1[,i] <- colMeans(X1[match_gene1_idx[cur_gene1_cluster - n_sim],])
        } else if (length(cur_gene1_cluster) == 1 ) {
          MCV1[,i] <- X1[match_gene1_idx[cur_gene1_cluster - n_sim],]
        }
      }

      X1_uns_sim = X1[gene1_uns_idx,]
      X2_uns_sim = X2[gene2_uns_idx,]

      # unassigned_cluster1 <- dist.matrix(X1_uns_sim, t(MCV1), method = "euclidean")
      unassigned_cluster1 <- dista(X1_uns_sim, t(MCV1), type = "euclidean")
      cluster_uns_g1 <- vector(length = nrow(unassigned_cluster1))
      # for (j in 1:nrow(unassigned_cluster1)) {
      #   cluster_uns_g1[j] <- which.min(unassigned_cluster1[j,])
      # }
      cluster_uns_g1 <- rowMins(unassigned_cluster1, value = FALSE) # changed na.rm = T

      # unassigned_cluster2 <- dist.matrix(X2_uns_sim, t(MCV2), method = "euclidean")
      unassigned_cluster2 <- dista(X2_uns_sim, t(MCV2), type = "euclidean")
      cluster_uns_g2 <- vector(length = nrow(unassigned_cluster2))


      # stt <- Sys.time()
      # browser()
      cluster_uns_g2 <- rowMins(unassigned_cluster2, value = FALSE)
      # print(Sys.time() - stt)
      # stt <- Sys.time()
      # for (j in 1:nrow(unassigned_cluster2)) {
      #   cluster_uns_g2[j] <- which.min(unassigned_cluster2[j,])
      # }
      # print(Sys.time() - stt)


      unsampled_gene1_cluster <- cbind(cluster_uns_g1, gene1_uns_idx)
      sampled_gene1_cluster <- cbind(cluster_sim$cluster[(n_sim+1) : (n_sim+m_sim)], match_gene1_idx)
      combined_gene1_cluster <- rbind(unsampled_gene1_cluster, sampled_gene1_cluster)
      sorted_gene1_cluster <- combined_gene1_cluster[order(combined_gene1_cluster[,2]),]

      unsampled_gene2_cluster <- cbind(cluster_uns_g2, gene2_uns_idx)
      sampled_gene2_cluster <- cbind(cluster_sim$cluster[1 : n_sim], match_gene2_idx)
      combined_gene2_cluster <-rbind(unsampled_gene2_cluster, sampled_gene2_cluster)
      sorted_gene2_cluster <- combined_gene2_cluster[order(combined_gene2_cluster[,2]),]
      sorted_gene2_cluster[,2] <- sorted_gene2_cluster[,2] + num_m

      combined_g1_g2_cluster <- rbind(sorted_gene1_cluster, sorted_gene2_cluster)

      # add Xuechun & Zhen's M #
      X_clust <- combined_g1_g2_cluster[,1]
      clust = Outer(as.numeric(X_clust),as.numeric( X_clust), oper  = "/")
      clust[clust != 1] <- 0
      M = clust
      M.diff = (M - avg.M)
      avg.M = avg.M + M.diff / h
    }
    return(avg.M)
  }

  set.seed(20221128)
  avg.M <- sample_func(niter_0, rho_0)

  Hclust1 = hclust(as.dist(1-avg.M), method = "complete")
  # plot dendrogram
  hclust_plt = plot(Hclust1)
  tree = cutree(Hclust1, h = 1 - alpha)
  # plot ecdf function
  Fn = ecdf(avg.M)
  ecdf_plt = plot(Fn, verticals = FALSE,  col.hor = "red", xlab = "concensus value", ylab = "empirical CDF", main = "", cex = 0.5)

  # txt file of co-clusters:
  #cut_avg = data.frame(tree)
  #cut_avg = cbind(GeneId = rownames(cut_avg), cut_avg)
  #rownames(cut_avg) = 1:nrow(cut_avg)
  #colnames(cut_avg) = c("GeneId","Cluster")
  #output = NULL
  #for(i in 1:len_cluster){
  #  output[[i]] = cut_avg2$GeneId[cut_avg$Cluster == i]
  #} 
  #opt_file = write.table(output, file = "C://Users//15588//OneDrive//Desktop//cutree04.txt",sep = "\t", row.name = TRUE)
  
  my_return = list("ecdf_plt" = ecdf_plt, "avg.M" = avg.M, "Fn" = Fn, "hclust_plt" = hclust_plt)
  return(avg.M)
}
# library(readr)
# library(dplyr)
# library(Matrix)
# library(matrixLaplacian)
# library(wordspace)
# library(matrixStats)
# library(profvis)
# library(rARPACK)
# library(bench)
# library(Rfast)
# library(IOHanalyzer)
# tst_m <- bitsc_v1(0.01, 3)
# devtools::install_github("Orion-qx/BiTSCR")
