#' Select top features
#'
#' @export
#' @param data.splsda A sPLS-DA object from mixomics.
#' @param feature_names A character vector containing feature names.
#' @param component A character string containing sPLS-DA's component name.
#' @param Y_name A character string indicating which is the variable associated with component specified. Must be a dicothomous variable.
#' @param n_top An integer number indicating how many top associated features to keep. (Default = 15)
#' @param offset A value indicating the amount of justification from axis for possible plots. (Default = 0.005)
#' @param order A boolean value indicating if the loadings should be separately ordered by positives and negatives (Default = TRUE). If FALSE loadings are ordered using their absolute value.
#' @return A data frame containings \code{n_top} feature associated with \code{component} with information about interesting \code{Y_name} variable.
#' @seealso \code{\link{HotLoadings.names}} for extraction of \code{feature_names}

HotLoadings.top_features <- function(data.splsda, feature_names, component = c("comp 1"),Y_name,n_top = 15,offset = 0.005, order = TRUE){
  # Generate a df with loadings and feature names for the selected component
  loadings_df <- data.frame(loadings = data.splsda$loadings$X[,component],feature = feature_names)
  loadings_df <- loadings_df[order(abs(loadings_df$loadings),decreasing = TRUE),]
  # Generate variable associated with component
  # Check Y loadings
  Y_loadings <- data.splsda$loadings$Y[,component]
  stronger_index <- c(which(Y_loadings == max(Y_loadings)),which(Y_loadings == min(Y_loadings)))
  Y_loadings <- Y_loadings[stronger_index]
  loadings_df$comparison <- ifelse(loadings_df$loadings < 0,names(Y_loadings)[which(Y_loadings < 0)],names(Y_loadings)[which(Y_loadings > 0)])
  # Selection top n features
  top_features <- loadings_df[1:n_top,]
  # Ordering
  if(order)
    top_features <- top_features[order(top_features$loadings),]
  # Convert to factor to retain sorted order in possible plot.
  top_features$feature <- factor(top_features$feature, levels = top_features$feature)
  # Generate a variable for text justification
  top_features$just <- ifelse(top_features$loadings>0,"right","left")
  colnames(top_features) <- c("loadings","feature",Y_name,"just")
  top_features$offset <- ifelse(top_features$just=="left",offset,-offset)
  return(top_features)
}
