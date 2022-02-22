#' @title Apply MDS plot
#' @description Function for creating MDS plots
#' @param data a gene counts matrix
#' @param plotmds
#' @param databatch TRUE/FALSE parameter for batch correction
#' @param data.batch batch corrected gene counts
#' @param group a group information for each sample
#' @param ids sample identifiers from metadata table
#' @param prefix a name of the results folder
#' @param technology a technology
#' @param colors colors
#' @export
#' @seealso
#' @return MDS plots
#' @examples \dontrun{
#' ...
#' }
#'

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
