#' @title Estimate Kmeans
#' @description A function for estimating K-means.
#' Number of kmeans will be based on the bayesian information criterion(BIC)
#' provided by mclust package. Summary describing the best model is printed on
#' the screen during the calculation.
#' @param data a dataframe of feature (e.g. gene) counts
#' @export
#' @import mclust
#' @seealso
#' @return number of clusters
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

    return(n.clusters)

}
