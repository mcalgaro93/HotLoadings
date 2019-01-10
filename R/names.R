#' Extract feature names
#'
#' @param PSOBJ A phyloseq object.
#' @param format "last" if only the deepest taxonomic characterization is needed. "long" if all taxonomic levels are needed. "short" if last 2 known taxonomic levels are needed.
#' @return A data frame containings feature names in wanted format.
#' @examples
#' library(phyloseq)
#' data("GlobalPatterns")
#' feature_names <- HotLoadings.names(PSOBJ = GlobalPatterns, format = "short")
#' head(feature_names)

HotLoadings.names <- function(PSOBJ,format = c("last","long","short")){
  # Pasting togheter names from Kingdom to the deepest characterisation for each taxa
  extended_names <- apply(PSOBJ@tax_table@.Data,1,function(row){
    index <- which(!is.na(row))
    # 1:Kingdom,2:Phylum,3:Class,4:Order,5:Family,6:Genus,7:Species
    levels <- c("k:","p:","c:","o:","f:","g:","s:")
    switch (format,
      last = return(paste0(levels[max(index)],row[max(index)],collapse = ",")),
      long = return(paste0(levels[index],row[index],collapse = ",")),
      short = return(paste0(levels[c(max(index)-1,max(index))],row[c(max(index)-1,max(index))],collapse = ",")))
  })
  # Pasting SV number too
  if(!is.null(names)){
    extended_names <- paste(names(extended_names),extended_names,sep = "|")
  }
  return(extended_names)
}

