#+------------------------------------------------------------------------+
#+                                                                        +
#+          u   U   NN  N  IIIIII  PPPP   RRRR    OOO   TTTTTT            +
#+          U   U   NNN N    II    P   P  R   R  O   O    TT              +
#+          U   U   N N N    II    PPPP   RRRR   O   O    TT              +
#+          U   U   N  NN    II    P      R  R   O   O    TT              +
#+           UuU    N   N  IIIIII  P      R   R   OOO     TT              +
#+                                                                        +
#+                Algorithm for automatic annotation                      +
#+         Modified from Programmatic Access - ID Mapping                 +
#+                            Version 1.0                                 +
#+                      By: Alex Fernando Arita                           +
#+                            2024-04-15                                  +
#+                                                                        +
#+ This script receives a directory/folder where ".txt" files are stored  +
#+ and generates a new folder called "results" where annotations are      +
#+ stored for each file.                                                  +
#+                                                                        +
#+ The files must have the format:                                        +
#+  - INFORMATION: Must have a # at the start of each line,               +
#+  - GENES: Must be named in gene symbols and separated by "\n".         +
#+           Modifications to the code could be made to work with other   +
#+           annotations.                                                 +
#+                                                                        + 
#+------------------------------------------------------------------------+

library(httr)
library(readr)
library(dplyr)

#+------------------------------------------------------------------------+
#+                                                                        +
#+                              FUNCTIONS                                 +
#+                                                                        +
#+------------------------------------------------------------------------+                                                                        
                                                                          
# Function to check if job is ready
isJobReady <- function(jobId) {
  pollingInterval <- 30
  maxTries <- 10
  for (i in 1:maxTries) {
    url <- paste("https://rest.uniprot.org/idmapping/status/", 
                 jobId, 
                 sep = "")
    r <- GET(url = url, accept_json())
    status <- content(r, as = "parsed")
    if (!is.null(status[["results"]]) || !is.null(status[["failedIds"]])) {
      return(TRUE)
    }
    if (!is.null(status[["messages"]])) {
      print(status[["messages"]])
      return (FALSE)
    }
    Sys.sleep(pollingInterval)
    pollingInterval <- pollingInterval / 2  # Increase polling frequency
  }
  return(FALSE)
}

# Function to get results URL
getResultsURL <- function(redirectURL) {
  if (grepl("/idmapping/results/", 
            redirectURL, 
            fixed = TRUE)) {
    url <- gsub("/idmapping/results/", "/idmapping/stream/", redirectURL)
  } else {
    url <- gsub("/results/", "/results/stream/", redirectURL)
  }
  return(url)
}

# Function to process genes
processGenes <- function(genes, from, to) {
  genes_string <- paste(genes, collapse = ",")
  files <- list(
    ids = genes_string,
    from = from,
    to = to
  )
  r <- POST(url = "https://rest.uniprot.org/idmapping/run", 
            body = files, 
            encode = "multipart", 
            accept_json())
  submission <- content(r, as = "parsed")
  if (isJobReady(submission[["jobId"]])) {
    url <- paste("https://rest.uniprot.org/idmapping/details/", 
                 submission[["jobId"]], 
                 sep = "")
    r <- GET(url = url, accept_json())
    details <- content(r, as = "parsed")
    url <- getResultsURL(details[["redirectURL"]])
    url <- paste(url, "?format=tsv", sep = "")
    r <- GET(url = url, accept_json())
    resultsTable <- read_tsv(content(r), 
                             col_names = TRUE)
    
    # Filter results to keep only the first match for each input gene
    filtered_results <- resultsTable |>
      group_by(From) |>
      slice(which.min(ifelse(Organism == "Homo sapiens (Human)", 0, 1)))
    
    return(filtered_results)
  }
}

# Function to process all genes in a directory
processGenes <- function(genes, from, to) {
  genes_string <- paste(genes, collapse = ",")
  files <- list(
    ids = genes_string,
    from = from,
    to = to
  )
  r <- POST(url = "https://rest.uniprot.org/idmapping/run", 
            body = files, 
            encode = "multipart", 
            accept_json())
  submission <- content(r, as = "parsed")
  if (isJobReady(submission[["jobId"]])) {
    url <- paste("https://rest.uniprot.org/idmapping/details/", 
                 submission[["jobId"]], 
                 sep = "")
    r <- GET(url = url, accept_json())
    details <- content(r, as = "parsed")
    url <- getResultsURL(details[["redirectURL"]])
    url <- paste(url, "?format=tsv", sep = "")
    r <- GET(url = url, accept_json())
    resultsTable <- read_tsv(content(r), col_names = TRUE)
    
    # Filter results to keep only the first match for each input gene
    filtered_results <- resultsTable %>%
      group_by(From) %>%
      slice(which.min(ifelse(Organism == "Homo sapiens (Human)", 0, 1)))
    
    # Get the genes that were not found
    not_found_genes <- setdiff(genes, filtered_results$From)
    
    # Create a dataframe for not found genes with blank fields
    not_found_df <- data.frame(From = not_found_genes, 
                               Entry = rep(NA_character_, 
                                           length(not_found_genes)), 
                               stringsAsFactors = FALSE)
    
    # Combine the filtered results with the not found genes
    combined_results <- bind_rows(filtered_results, not_found_df)
    
    return(combined_results)
  }
}

processAllGenesInDirectory <- function(directory_path, 
                                       from = "Gene_Name", 
                                       to = c("UniProtKB-Swiss-Prot")) {
  files <- list.files(directory_path, pattern = "\\.txt$", full.names = TRUE)
  
  results_directory <- file.path(directory_path, "results")
  if (!file.exists(results_directory)) {
    dir.create(results_directory)
  }
  
  for (file in files) {
    lines <- readLines(file)
    headers <- lines[grepl("^#", lines)]
    lines <- lines[!grepl("^#", lines)]
    genes_list <- unlist(strsplit(lines, split = "\t"))
    results <- processGenes(genes_list, from, to)
    
    # Extract file name without extension
    file_name <- gsub("\\.txt$", "", basename(file))
    
    # Create results file name
    results_file <- file.path(results_directory, 
                              paste0(file_name, "_uniprot.txt"))
    
    # Write results to file
    writeLines(headers, results_file)
    write.table(results, file = results_file, 
                sep = "\t", 
                row.names = FALSE, 
                col.names = TRUE, 
                append = TRUE)
    
    cat("Processed", 
        length(genes_list), 
        "genes from file", 
        file, "\n")
  }
}

directory <- getwd()
processAllGenesInDirectory(directory)
