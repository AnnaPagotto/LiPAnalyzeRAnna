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
