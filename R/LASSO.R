    #' @title LASSO function
    #' @description Run LASSO function to obtain DA features
    #' @param seed a randomly generated vector of seeds
    #' @param data a data frame of expression/abundance counts
    #' @param group a vector of integers specifying group
    #' @param alpha a hyperparameter alpha. This value must be set to 0.0 < x < 1.0 for Elastic Net (0.5 is default) or to 1.0 for LASSO regression.
    #' @param validation perform validation TRUE/FALSE
    #' @param multinorm perfom multinorm TRUE/FALSE; IF multinorm=TRUE, then analysis will be multinomial, IF multinorm=FALSE, then analysis will be binomial
    #' @export
    #' @import glmnet
    #' @import parallel
    #' @import doMC
    #' @seealso
    #' @return DA/DE features from LASSO (add format)
    #' @examples \dontrun{
    #' ...
    #' }

    ## Is setting the seed necessary?
    LASSOFeature2 <- function(seed, data, group, alpha, validation=FALSE, multinorm=TRUE) {

        if (validation == TRUE) {  ###REALLY NOT SURE IF THIS IS NECESSARY

            ll <- list()
            llev <- levels(as.factor(group))

            ##creates a list of indeces for each group
            for (idx in 1:length(llev)) {
                pos <- which(group == as.character(llev[idx]))
                ll[[idx]] <- pos
            }

            ##create a random sample selection (1/4 of the samples from each group in the input data). This is only for the calc of errors
            samp <- unlist(lapply(ll, function(x) sample(x, ceiling((length(x)/4)))))

            data.val <- data[,samp]
            group.val <- as.character(group[samp])  ##originally as.integer

            data <- data[,-samp]
            group <- as.character(group[-samp])  ##originally as.integer
        }

        ## Is OK for this step that data are sub-sampled in case of validation is TRUE??? (probably yes - we need a dataset which is large enough)
        ## - cross validation is working BUT results generate less genes
        if(multinorm == TRUE) {
            set.seed(seed)   ###do we need this?
            nth <- detectCores(logical = TRUE)
            registerDoMC(cores=nth)
            fit <- cv.glmnet(x = t(data), y = group, family="multinomial", type.multinomial = "grouped", nfolds = 10, alpha = alpha, parallel=TRUE)
            print("fit")
            print(head(fit))
            print(fit$lambda.1se)
            print(fit$lambda.min)
            ###add printing of minimum cvm and

            coef <- coef(fit, s=fit$lambda.1se)
            # print("coef")
            # print(head(coef))
            ma <- as(coef[[1]], "matrix")
            print("ma")
            print(head(ma))
        } else {
            set.seed(seed)
            fit <- cv.glmnet(x = t(data), y = group, family = "binomial", type.measure = "class", nfolds = 10, alpha = alpha)
            coef <- coef(fit, s=fit$lambda.1se)
            ma <- as(coef, "matrix")
        }

        ###Double cross-validation
        if (validation == TRUE) {
            ## Here we calculate miss classification error.
            meanerror <- round(as.numeric(mean(predict(fit, t(data.val), s=fit$lambda.1se, type="class") != group.val))*100, digits = 2)

            print("PREDICT")
            print(predict(fit, t(data.val), s=fit$lambda.1se, type="class"))
            print(predict(fit, t(data.val), s=fit$lambda.1se, type="class") != group.val)
            print(mean(predict(fit, t(data.val), s=fit$lambda.1se, type="class") != group.val))
            print("cvm")
            print(fit$cvm)

            cat(paste0("\nOne LASSO/EN run out completed. Classification error rate on test data was = ", meanerror, "%\n"))
            ## here the fit from the data is applied on data.val
            ## I think cross validation error should be mentioned here, not in the case below
            ## Study this part more - fit variable from the previous steps

        } else {
            meanerror <- round(as.numeric(mean(predict(fit, t(data), s=fit$lambda.1se, type="class") != group))*100, digits = 2)
            cat(paste0("\nOne LASSO/EN run out completed. Mean cross-validation error was = ", meanerror, "%\n"))
            ## here the fit from the data is applied on data again which results in 0.
            ## I think this should not be done in this case. Instead of this we should report that mean cross validation error was not calculated due to the insufficient nubmer of samples
        }


        ##remove !=0; removing genes with coef=0 could be done in the main function
        ma <- names(ma[ma[,1] != 0, ])

        ## return genes and coefficients
        return(list(ma, meanerror))

    }


