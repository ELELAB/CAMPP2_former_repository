#' @title Estimate Kmeans
#' @description A function for estimating the k-means clusters. The clusters
#' are visualised on a PCA graph.
#' In case of data sets with >1000 features,
#' multiple sets (number of sets = number of features/1000; rounded up to the
#' integer) will be generated
#' and tested for the number of clusters estimation using BIC.
#' Consensus of best n k-means will be returned. Number of randomly selected
#' features in 1 set
#' is limited to 2000 (for large sets having more than 2000 features).
#' In case k-means estimation
#' using BIC fails, number of clusters will be based on the number of samples
#' (2-6 clusters for dataset with less than 100 samples; 2-11 for datasets with
#' 101-500 samples and 2-16 clusters for data sets with more than 500 samples).
#' The clusters visualized on a PCA graph are saved into the png.
#' @param data a data.frame of feature (e.g. gene) counts
#' @param PCA.labels a text ("all", "none) specifying the elements to be
#' labelled. Default value is "none".
#' @param cols a vector of colors (one color for each group)
#' @param prefix a character defining a prefix of output file.
#' @export
#' @import factoextra
#' @import FactoMineR
#' @seealso
#' @return a data.frame with cluster info assigned to each sample;
#' Figures saved into png.
#' Results from PCA
#' @examples \dontrun{
#' runKmeans(campp2_brca_1_normalized$E, PCA.labels = FALSE, cols=NULL, prefix="test_run")
#' }

runKmeans <- function(data, PCA.labels = FALSE, cols=NULL, prefix=NULL){

     if(is.null(prefix)){
        stop(print("Please, provide a prefix for the result files."))
     }

    # Number of sample sets to generate
    n.subsets <- 1:ceiling(nrow(data)/1000)
    if(length(n.subsets) > 10) {
        n.subsets <- 1:10
    }

    # Number of variables (genes) in each sample set
    setsize <- nrow(data)
    if (setsize > 2000) {
        setsize <- 2000
    }

    # Number of kmeans to try if k-means estimation using BIC is not working
    if(ncol(data) <= 100) {
        nks <- 2:6
    } else if (ncol(data) > 100 && ncol(data) <= 500) {
        nks <- 2:11
    } else {
        nks <- 2:16
    }


    list.of.dfs <- list()

    for (idx in 1:length(n.subsets)) {
        df <- t(data[sample(nrow(data), setsize), ])  #make random gene selections
        list.of.dfs[[idx]] <- df  #create list of random selections
    }
    cat(paste0("\n", length(n.subsets), " clustering runs will be done in total...\n"))
    n.clusters.list <- lapply(list.of.dfs, function(x) EstimateKmeans(x)) #estimate how many clusters in the data subsets


    nclus <- unique(unlist(n.clusters.list))

    if (unique(is.na(nclus)) == TRUE) {
        cat(paste0("Number of clusters could not be determined using BIC. There may be little or poor clustering of samples. Alternatively, based on size of dataset, ", length(n.subsets), " sample sets will be generated of size ", setsize, " and ", length(nks), " clusters will be tested. \nRunning..."))

        nclus <- nks
    }
    nclus <- sort(nclus)
    paste0("Numbers of clusters being tested: ",nclus)

    res.pca <- PCA(t(data),  graph = FALSE, ncp=10, scale = FALSE) # principal component analysis
    res.list <- list()

    for (idx in 1:length(nclus)) {
        set.seed(10)
        Kclus <- kmeans(t(data), nclus[[idx]])
        Clusters <- as.factor(paste0("C",data.frame(Kclus$cluster)$Kclus.cluster))
        res.list[[idx]] <- Clusters
        if(length(cols) < nclus[[idx]]){
            print("Number of colours is defined by the number of groups and is smaller than the number of clusters. Colour scheme will be defined automatically.")
            cols <- NULL
        }
        fviz_pca_ind(res.pca,
                     label = PCA.labels, # show/hide individual labels; labels are taken from feature counts matrix automatically
                     habillage = as.factor(Clusters), # color by groups (clusters in this case)
                     palette = cols,
                     addEllipses = TRUE, # concentration ellipses
                     repel=TRUE,
                     ggtheme = theme_classic(),
                     title = paste0("k-means ", nclus[[idx]], " clusters"),
                     labelsize = 2

        )
        ggsave(paste0(prefix,"_PCA_Kmeans_C",nclus[[idx]],".png"))

    }
    names(res.list) <- paste0("Clus", nclus)
    res.clusters <- as.data.frame(res.list)
    return(list("res.clusters"=res.clusters,"res.pca"=res.pca))

}
