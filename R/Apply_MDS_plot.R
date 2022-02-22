#' @title Correlation analysis
#' @description If two datasets are provided (for correlation and/or network analysis), this option should be specified as a comma separated list (without quotes or parenthesis!) of length two,
#' first entry referring to data file 1 and second entry referring to the data file 2.
#' @param my.d1 first dataframe with values to be correlated. These must have the same dimensions and rownames must the same.
#' @param my.d2 second dataframe with values to be correlated. These must have the same dimensions and rownames must the same.
#' @param my.filename a name of output plot
#' @export
#' @import ggplot2
#' @import plyr
#' @seealso
#' @return correlation plots
#' @examples \dontrun{
#' ...
#' }

applyMDSPlot<-function(data=data, plotmds=plotmds, databatch=databatch, data.batch=data.batch, group=group, ids=ids, prefix=prefix, technology=technology, colors=colors){

    MDScolors <- gsub(pattern = "FF", replacement = "", x = colors)


    if (plotmds == TRUE && databatch == TRUE){
        mdsplot <- MDSPlot(data.batch, group, ids, MDScolors)
        ggsave(paste0(prefix, "_MDSplot_batchcorr.pdf"), plot = mdsplot, dpi = 300, width = 8, height = 8)

    } else if (plotmds == TRUE && databatch == FALSE){
        if (technology[1] == "seq") {
            mdsplot <- MDSPlot(data.frame(data$E), group, ids, MDScolors)
        } else {
            mdsplot <- MDSPlot(data, group, ids, MDScolors)
        }
        ggsave(paste0(prefix, "_MDSplot.pdf"), plot = mdsplot, dpi = 300, width = 8, height = 8)

        rm(mdsplot)

    } else {
        cat("\n- No preliminary plot requested.\n")
    }
}
