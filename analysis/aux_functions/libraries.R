# Load a vector of libraries, installing missing ones
libraries <- function(pkgs) {
    if(!require(BiocManager)) install.packages("BiocManager")
    #-- Get uninstall libraries
    new_libraries <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
    #-- Install using BiocManager
    if(length(new_libraries)) BiocManager::install(new_libraries)
    #-- Load
    invisible(lapply(pkgs, require, character.only = TRUE))
}