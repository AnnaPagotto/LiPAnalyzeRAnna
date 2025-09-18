## Anna: in the contrast, run wilcox test instead of t-statistics that you do with linear regression. In case after the RUV the residuals are not normally distributed.

globalVariables(names=c("XPep", "XProt", "Y"))


runModel_wilcox <- function(quantityList, annotS=NULL, formulaRUV="Y~XPep+XProt",
                     formulaContrast="Y~Condition", lowRUV=c(-1e9, 0, 0),
                     upRUV=c(Inf, Inf, Inf), addRUVbounds=FALSE,
                     returnRUVmodels=FALSE, returnContrastmodels=FALSE){

    ## if necessary transforming formulas into character
    ## remove spaces in formula
    if(!is.null(formulaRUV)){
        if(inherits(formulaRUV,"formula")){
            formulaRUV <- Reduce(paste, deparse(formulaRUV))
        }
        formulaRUV <- gsub(" ", "", formulaRUV)
    }
    if(!is.null(formulaContrast)){
        if(inherits(formulaContrast,"formula")){
            formulaContrast <- Reduce(paste, deparse(formulaContrast))
        }
        formulaContrast <- gsub(" ", "", formulaContrast)
    }

    ## check input format of formulas
    if(!(is.null(formulaRUV)|is.character(formulaRUV))&
       (is.null(formulaContrast)|is.character(formulaContrast))){
        stop("Please provide 'formulaRUV' and 'formulaConstrast' in the correct
class. Use 'character' or 'formula'.")
    }
    ## stop if no formulas are provided
    if(is.null(formulaRUV) & is.null(formulaContrast)){
        stop("No RUV or contrast formula provided. Please provide at least one
of them for running models.")
    }
    resAll <- NULL

    ## checking if quantiyList is in data.frame/matrix format
    if(!all(unlist(lapply(quantityList, \(x)
                          inherits(x, c("matrix","data.frame")))))){
        stop("Elements of 'quantityList' have to be data.frames or
matrices.")
    }
    quantityList <- lapply(quantityList, as.data.frame)


    ## checking and potentially changing names of quantityList matrices
    if(sum(!names(quantityList) %in% c("LiPPep", "TrPPep", "TrPProt", "LiPProt",
                                       "Y", "XPep", "XProt"))>0){
        stop("Names of matrices in 'quantityList' not permitted. Please change
accordingly.")
    }

    if(paste(names(quantityList), collapse="") == c("LiPPepTrPPepTrPProt")){
        names(quantityList) <- c("Y", "XPep", "XProt")
    }

    if(paste(names(quantityList), collapse="") == c("LiPPepTrPPep")){
        names(quantityList) <- c("Y", "XPep")
    }

    else if(paste(names(quantityList), collapse="") == c("LiPPepLiPProt")|
            paste(names(quantityList), collapse="") == c("LiPPepTrPProt")){
        names(quantityList) <- c("Y", "XProt")
    }

    ## checking that max 1 Y, XPep and XProt are provided in quantityList
    if(sum(duplicated(names(quantityList)))>0){
        if(sum(grepl("Y", names(quantityList)))>1){
            stop("Multiple 'LiPPep'/'Y' matrices provided in the quantityList.
Please only povide one.")
        }
        if(sum(grepl("XPep", names(quantityList)))>1){
            stop("Multiple 'TrPPep'/'XPep' matrices provided in the
quantityList. Please only povide one.")
        }
        if(sum(grepl("XPep", names(quantityList)))>1){
            stop("Multiple 'TrPProt'/'LiPProt'/'XProt' matrices provided in the
quantityList. Please only povide one.")
        }
    }

    ###quantityList <- lapply(quantityList, as.data.frame)

    ## assuring row.names and colnames fit over all input data
    feat <- Reduce(intersect, lapply(quantityList, row.names))
    samples <- Reduce(intersect, lapply(quantityList, colnames))

    if(!is.null(annotS)){
        samples <- intersect(row.names(annotS), samples)
    }
    if(length(feat) == 0|length(samples) == 0){
        stop(length(feat), " peptides/proteins in ", length(samples), " samples
detected.\nPlease check your input as well as column and row names of all
provided data." )
    }
    message(length(feat), " quantities and ", length(samples),
            " samples used in models.")
    quantityList <- lapply(quantityList, function(x){
        x[feat, samples]
    })
    if(!is.null(annotS)){
        annotS <- annotS[samples, ]
    }

    ## running RUV models
    if(!is.null(formulaRUV)){
        message("Running RUV models.")
        modelMat <- createModelMatrix(quantityList, formulaRUV, annotS, samples)
        modelRUV <- runRUV(formulaRUV, modelMat, lowRUV, upRUV,
                             addRUVbounds)
        resRUV <- extractRUV(modelRUV, samples)
        if(is.null(formulaContrast)){
            resAll <- resRUV
            message("Returning RUV results including residuals and estimated
coefficients.")
        }
        else{
            dfRUV <- ncol(resRUV[[2]])-1
            quantityList[["Y"]] <- resRUV[[1]]
        }
    }
    else{
        dfRUV <- NULL
    }

    ## running contrast models
    if(!is.null(formulaContrast)){
        message("Running contrast models.")
        modelMat <- createModelMatrix(quantityList, formulaContrast, annotS,
                                      samples)
        modelContrast <- runContrast(modelMat)
        resContrast <- extractContrast(modelContrast, formulaContrast, dfRUV)
        if(!is.null(formulaRUV)){
            resContrast$modelCoeff <- cbind(resRUV$modelCoeff,
                                            resContrast$modelCoeff)
        }
        resAll <- resContrast
    }

    ## adding feature names to output
    resAll <- lapply(resAll, function(x){
        rownames(x) <- feat
        return(as.data.frame(x))
    })

    ## add message if there TrPPep or TrPProt coefficients are very high
    if(any(grepl("XPep|XProt", colnames(resAll[[1]])))){
        coeffPepProt <- unlist(c(resAll[[1]][, grepl("XPep",
                                                     colnames(resAll[[1]]))],
                                 resAll[[1]][, grepl("XProt",
                                                     colnames(resAll[[1]]))]))
        if(max(stats::na.omit(coeffPepProt))>5){
            message("At least one peptide/protein coefficient estiamted is
higher than 5, please manually check the results of the respective peptides for
plausibility.")
        }
    }

    ## add residuals to output, if not already included
    if(!is.null(formulaRUV)&!(is.null(formulaContrast))){
        resAll[[length(resAll)+1]] <- quantityList[["Y"]]
        names(resAll)[length(resAll)] <- "modelResid"
    }
    ## adding models to results if set in function input
    if(returnRUVmodels){
        names(modelRUV) <- row.names(resAll$modelCoeff)
        resAll[[length(resAll)+1]] <- modelRUV
        names(resAll)[length(resAll)] <- "modelRUV"
    }

    if(returnContrastmodels){
        names(modelContrast) <- row.names(resAll$modelCoeff)
        resAll[[length(resAll)+1]] <- modelContrast
        names(resAll)[length(resAll)] <- "modelContrast"
    }

    return(resAll)

}

#' @title Creating model matrices to run RUV or contrast models
#'
#' @description Creates one model matrices per peptide/protein which are
#' used in the RUV or contrast model function afterwards
#'
#' @usage createModelMatrix(quantityList, formula, annotS, samples)
#'
#' @param quantityList A list of preprocessed matrices, containing quantities of
#' interest(e.g. peptide, modified peptide, precursor) and protein abundances.
#' Rows represent features and columns samples and should match between the
#' different matrices contained in the list.
#' Output from \code{preprocessQuantityMatrix} have to use the variable naming
#' 'Y', 'XPep' and 'XProt'.
#' @param formula A formula providing structure of model matrices created in
#' this function
#' @param annotS A data.frame containing sample annotation. Must contain all
#' columns included in the RUV and contrast models. Rows are samples and must
#' match to columns of the matrices in \code{quantityList}. Must include
#' columns of any further variables used in \code{formulaRUV}.
#' @param samples A character vector providing sample names for the model.
#'
#' @return A list with model matrices for running the RUV or contrast models.

createModelMatrix <- function(quantityList, formula, annotS, samples){

    list2env(quantityList, envir=environment()) ## write matrices to environment

    ## create model matrices for each quantity
    modelMat <- lapply(as.list(seq(1, nrow(Y))), function(i){
        Y <- as.numeric(Y[i, ])
        formula <- stats::as.formula(formula)
        formulaVars <- formula.tools::get.vars(formula)[
            formula.tools::get.vars(formula)!="Y"]
        e <- environment()
        sapply(formulaVars, function(x){
            if(x == "XPep"){
                assign(x, as.numeric(XPep[i,]), envir=e)
            }
            else if(x == "XProt"){
                assign(x, as.numeric(XProt[i,]), envir=e)
            }
            else{
                assign(x, annotS[,x], envir=e)
            }
            return(NULL)
        })

        X <- stats::model.matrix(formula)
        attr(X, "samples") <- samples[as.numeric(row.names(X))]
        return(list(Y=Y[as.numeric(row.names(X))], X=X))
    })

    return(modelMat)
}


#' @title Running RUV models for every quantity provided
#'
#' @description Function to run RUV on a list of model matrices, one element in
#' the list should represent one peptide or protein.
#'
#' @usage runRUV(formula, modelMat, lowRUV, upRUV, addRUVbounds)
#'
#' @param formula A formula used to create \code{modelMat} for RUV models.
#' @param modelMat A list of model matrices to perform RUV using bounded
#' variable least square regression.
#' @param lowRUV A numeric vector defining lower boundaries of the coefficients
#' of the RUV models. Elements refer to definition of \code{formulaRUV}.
#' Default is 'c(-Inf, 0, 0)'.
#' @param upRUV A numeric vector defining upper boundaries of the coefficients
#' of the RUV models. Elements refer to definition of \code{formulaRUV}.
#' Default is 'c(Inf, Inf, Inf)'.
#' @param addRUVbounds A boolean value, if set to 'TRUE' as many bounds as
#' additionally needed based on \code{formulaRUV} in each RUV model are added to
#' \code{lowRUV} and \code{upRUV}. Added boundaries are automatically set to
#' \code{lowRUV = -Inf} and \code{upRUV = Inf}. Important to set to 'TRUE', if
#' you have categories with multiple levels in the RUV model and did not adjust
#' the RUV boundaries based in the number of levels. This might be the case if
#' you have more than two batches you aim to account for in the RUV model.
#' Default is 'FALSE'.
#'
#' @return Returns a list of all RUV models

runRUV <- function(formula, modelMat, lowRUV, upRUV, addRUVbounds){

    ## change -Inf to high negative number instead, since bvls() does not take
    ## -Inf as an input
    lowRUV[lowRUV == -Inf] <- -1e9

    modelRUV <- lapply(modelMat, function(data){

        Y <- data$Y
        X <- data$X

        ## adding as many as necessary -Inf and Inf to the boundaries for bvls()
        ## if 'addRUVbounds' is set to TRUE
        if(addRUVbounds){
            nX <- ncol(X)
            lowRep <- nX-length(lowRUV)
            lowRUV <- c(lowRUV, rep(-1e9, lowRep))
            upRep <- nX-length(upRUV)
            upRUV <- c(upRUV, rep(Inf, upRep))
        }

        ## running RUV
        RUV <- bvls::bvls(A=X,
                           b=Y,
                           bl= lowRUV,
                           bu=upRUV)

        ## add variable and sample annotation to RUV output
        varsRUV <- paste0(dimnames(X)[[2]], "_RUV")
        varsRUV[1] <- "Intercept_RUV"
        attr(RUV, "variables") <- varsRUV
        attr(RUV, "samples") <- attributes(X)$samples
        return(RUV)
    })
    return(modelRUV)
}

#' @title Extracting information from RUV models
#'
#' @description Function to extract the coefficients and residuals from each
#' RUV model
#'
#' @usage extractRUV(mRUV, samples)
#'
#' @param mRUV A list of RUV models estimated using bounded variable least
#' square regression
#' @param samples A character vector providing sample names for the model.
#'
#' @return A list were the first element is a data.frame with residuals from
#; the RUV models and the second element is a data.frame with coefficients
#' from the RUV models.

extractRUV <- function(mRUV, samples){

    ## extract residuals and coefficients from RUV models
    modelResid <- do.call(plyr::rbind.fill, lapply(mRUV, function(x){
        y <- as.data.frame(t(x$residuals))
        colnames(y) <- attributes(x)$samples
        return(y)
    }))

    modelCoeff <- do.call(plyr::rbind.fill, lapply(mRUV, function(x){
        y <- as.data.frame(t(x$x))
        colnames(y) <- attributes(x)$variables
        return(y)
    }))

    ## adding NA columns in residuals if peptide in sample could not be modeled
    startDf <- as.data.frame(matrix(NA, nrow=1, ncol=length(samples)))
    colnames(startDf) <- samples
    modelResid <- plyr::rbind.fill(startDf, modelResid)[-1,]

    return(list(modelResid=modelResid, modelCoeff=modelCoeff))
}


#' @title Running contrast models for every quantity provided
#'
#' @description Function to run contrast on a list of model matrices, one
#' element in the list should represent one peptide or protein.
#'
#' @usage runContrast(modelMat)
#'
#' @param modelMat A list of model matrices to perform contrast modeling on
#' facilitating ordinary least square regression.
#'
#' @return Returns a list of all contrast models

runContrast <- function(modelMat){
    modelRes <- lapply(modelMat, function(data){
        Y <- data$Y
        X <- data$X[, -1, drop=FALSE]
        stats::lm(Y ~ X)
    })
    return(modelRes)
}


#' @title Extracting information from contrast models
#'
#' @description Function to extract the coefficients and residuals from each
#' RUV model. If RUV was run before the contrast model, the p-values estimation
#' taks into account the degrees of freedom already used by the RUV step.
#'
#' @usage extractContrast(mContrast, formulaContrast, dfRUV, coeffPval="Pr(>|t|)",
#' coeffTval="t value")
#'
#' @param mContrast A list of contrast models
#' @param formulaContrast A formula used to create \code{modelMat} for
#' contrast models.
#' @param dfRUV A numberic value providing the number of degrees of freedom used
#' in the RUV models.
#' @param coeffPval A character variable giving the column name were p-values
#' are provided if \code{summary.lm(mRUV)}.
#' Default is "Pr(>|t|)".
#' @param coeffTval A character variable giving the column name were t-values
#' are provided if \code{summary.lm(mRUV)}.
#' Default is "t value".
#'
#' @return A list were the first element is a data.frame with residuals from
#; the contrast models and the second element is a data.frame with coefficients
#' from the contrast models.

extractContrast <- function(mContrast, formulaContrast, dfRUV,
                            coeffPval="Pr(>|t|)", coeffTval="t value"){

    ## extract coefficients from contrast models
    modelCoef <- do.call(plyr::rbind.fill, lapply(mContrast, function(x){
        as.data.frame(t(stats::coef(x)))
    }))

  
    ## if RUV was not run before, extract p-values directly from contrast models
    if(is.null(dfRUV)){
        modelPv <- do.call(plyr::rbind.fill, lapply(mContrast, function(x){
            as.data.frame(t(summary(x)$coefficients[, coeffPval]))
        }))
    }
    ## if RUV was run before, calculate p-values taking degrees of freedom
    ## used in RUV models into account
    else{
        message("Estimating p-values while removing degrees of freedom
previously used in the RUV models.")
        modelPv <- calcualtePvalAfterRUV(mContrast, coeffTval, dfRUV)
    }

    ## adjusting column names of coefficient and p-value data.frame
    if(ncol(modelCoef) == 2){
        formulaContrast <- stats::as.formula(formulaContrast)
        formulaVars <- formula.tools::get.vars(formulaContrast)[
            formula.tools::get.vars(formulaContrast)!="Y"]
        cols <- paste0(c("Intercept", formulaVars), "_Contrast")
    }
    else{
        cols <- colnames(modelCoef)[-1]
        cols <- substring(cols, 2)
        cols <- paste0(c("Intercept", cols), "_Contrast")
    }
    colnames(modelCoef) <- cols
    colnames(modelPv) <- cols
 

    return(list(modelCoeff=modelCoef, modelPv=modelPv)) 
}


#' @title Estimate p-values from contrast model taking degrees of freedom used
#' in RUV step into account
#'
#' @usage calcualtePvalAfterRUV(LM, coeffTval, dfRUV)
#'
#' @param LM A single linear regression model originating from contrast step
#' @param coeffTval A character variable giving the column name were t-values
#' are provided if \code{summary.lm(mRUV)}.
#' @param dfRUV A numberic value providing the number of degrees of freedom used
#' in the RUV models.
#'
#' @return Returns a data.frame with p-values.

calcualtePvalAfterRUV <- function(LM, coeffTval, dfRUV){
    modelPv <- do.call(plyr::rbind.fill.matrix, lapply(LM, function(x){
        x <- stats::summary.lm(x)
        df <- x$df[2] - dfRUV

        if(df<1){
            warning("Not enough degrees of freedom to estimate p-values, model
is not reliable! Returning NAs for affected quantity.")
            pv <- stats::setNames(rep(NA, nrow(x$coefficients)),
                                  row.names(x$coefficients))
            return(t(pv))
        }

        ## Estimating p-values from t-values of contrast models taking degrees
        ## of freedom used in RUV into account
        else{
            pv <- sapply(x$coefficients[, coeffTval], function(y){
                2*stats::pt(abs(y), df, lower.tail=FALSE)
            })
            return(t(pv))
        }
    }))

    return(modelPv)
}
