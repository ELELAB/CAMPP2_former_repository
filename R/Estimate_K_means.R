#' @title Estimate Kmeans
#' @description A function for estimating K-means.
#' Number of kmeans will be based on the bayesian information criterion(BIC);
#'
#' (this is more for runEstimateKmeans)
#' In case of data sets with many genes,
#' multiple samples of 3000 variables will be generated
#' and tested to overcome issues with computational time and the
#' consensus of best n kmeans will be returned.
#'
#' Clustering may only be performed one dataset at a time.
#'
#' (The number of clusters
#' tested will be based on number of samples, fewer samples will result in fewer
#' kmeans tested.)
#'
#' @param data a dataframe of feature (e.g. gene) counts
#' @export
#' @import mclust
#' @seealso
#' @return clustering r
#' @examples \dontrun{
#' ...
#' }

EstimateKmeans <- function(data) {
    BIC <- mclustBIC(data)
    plot(BIC)
    mod1 <- Mclust(data, x = BIC) #obtained model
    # fviz_nbclust(my_data, kmeans, method = "silhouette") #alternative approach
    n.clusters <- as.numeric(mod1$G)  #number of clusters according to the best model
    cat(paste0("number of clusters according to the best model, G: ",n.clusters,"\n"))
    summary.clust<-summary(mod1, parameters = FALSE) #summary for the best model
    print(summary.clust)
    # mclust::print.summary.Mclust(mod1)
    # return(list(clusters, mod1, summary.clust))

    return(n.clusters)

}
