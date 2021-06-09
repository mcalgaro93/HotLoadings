#' Extract feature names
#' @export
#' @param PSOBJ A phyloseq object.
#' @param format "last" if only the deepest taxonomic characterization is needed.
#' "long" if all taxonomic levels are needed. "short" if last 2 known taxonomic
#' levels are needed. "none" for the SV/OTU number only.
#' @return A data frame containings feature names in wanted format.

HotLoadings.names <- function(PSOBJ,format = c("last","long","short","none")){
  # Pasting togheter names from Kingdom to the deepest characterisation for each taxa
  extended_names <- apply(PSOBJ@tax_table@.Data,1,function(row){
    index <- which(!is.na(row))
    # 1:Kingdom,2:Phylum,3:Class,4:Order,5:Family,6:Genus,7:Species
    levels <- c("k:","p:","c:","o:","f:","g:","s:")
    switch (format,
      last = return(paste0(levels[max(index)],row[max(index)],collapse = ",")),
      long = return(paste0(levels[index],row[index],collapse = ",")),
      short = return(paste0(levels[c(max(index)-1,max(index))],row[c(max(index)-1,max(index))],collapse = ",")))
      none = return("")
  })
  # Pasting SV number too
  if(format != "none"){
    extended_names <- paste(names(extended_names),extended_names,sep = "|")
  } else {
    extended_names <- names(extended_names)
  }
  return(extended_names)
}

