  ##############################################
  # Assignment: DNA Sequence Analysis 
  # File: chr1_GL383518v1_alt.fa
  # Author: Puja
  ##############################################
  
  # --- Part 1 ---
  # --- Step 1: Setup and File Verification ---
  file_name <- "chr1_GL383518v1_alt.fa"
  
  cat("Checking working directory...\n")
  cat("Current working directory:", getwd(), "\n")
  
  # Verify file exists
  if (!file.exists(file_name)) {
    cat("ERROR: File not found at path:", file_name, "\n")
    cat("Ensure file is saved in the same directory shown above.\n")
    stop("Program stopped: Input FASTA file missing.")
  } else {
    cat("File found successfully.\n")
  }
  
  # --- Step 2: Read FASTA File ---
  sequence_lines <- readLines(file_name, warn = FALSE)
  
  # Remove header lines (those beginning with ">")
  sequence <- paste(sequence_lines[!grepl("^>", sequence_lines)], collapse = "")
  sequence <- gsub("\\s+", "", sequence)  # Remove whitespace/newlines
  
  # Convert to uppercase to handle case sensitivity
  sequence <- toupper(sequence)
  
  # --- Step 3: Validate Sequence ---
  if (nchar(sequence) == 0) {
    stop("Sequence appears empty or unreadable. Check FASTA file format.")
  } else {
    cat("Sequence read successfully.\n")
    cat("Sequence length:", nchar(sequence), "bases.\n")
  }
  
  # --- Step 4: Print Required Bases ---
  if (nchar(sequence) < 758) {
    stop("QC Error: Sequence too short to access position 758.")
  }
  
  base10 <- substr(sequence, 10, 10)
  base758 <- substr(sequence, 758, 758)
  
  cat("1.a 10th base:", base10, "\n")
  cat("1.b 758th base:", base758, "\n")
  
  ##############################################
  # --- Part 2 ---
  ##############################################
  
  cat("Reading FASTA file...\n")
  sequence_lines <- readLines(file_name, warn = FALSE)
  
  sequence <- paste(sequence_lines[!grepl("^>", sequence_lines)], collapse = "")
  sequence <- gsub("\\s+", "", sequence)
  
  # Convert to uppercase to handle case sensitivity
  sequence <- toupper(sequence)
  
  if (nchar(sequence) == 0) {
    stop("Sequence appears empty or improperly formatted.")
  } else {
    cat("Sequence read successfully.\n")
    cat("Sequence length:", nchar(sequence), "bases.\n")
  }
  
  cat("Creating reverse complement...\n")
  complement <- chartr("ATGC", "TACG", sequence)
  rev_complement <- paste(rev(strsplit(complement, NULL)[[1]]), collapse = "")
  
  if (nchar(rev_complement) != nchar(sequence)) {
    stop("Reverse complement length mismatch.")
  } else {
    cat("Reverse complement generated successfully.\n")
  }
  
  if (nchar(rev_complement) < 800) {
    stop("Sequence shorter than 800 bases. Cannot extract required range.")
  }
  
  base79 <- substr(rev_complement, 79, 79)
  segment_500_800 <- substr(rev_complement, 500, 800)
  
  cat("2.a 79th base of reverse complement:", base79, "\n")
  cat("2.b 500th–800th bases of reverse complement:\n", segment_500_800, "\n")
  
  ##############################################
  # --- Part 3 ---
  ##############################################
  
  if (!require(stringr)) {
    install.packages("stringr", dependencies = TRUE)
    library(stringr)
  }
  
  file_name <- "chr1_GL383518v1_alt.fa"
  
  cat("Checking working directory...\n")
  cat("Current working directory:", getwd(), "\n")
  
  if (!file.exists(file_name)) {
    cat("ERROR: File not found at path:", file_name, "\n")
    cat("Ensure FASTA file is located in the above directory.\n")
    stop("Program stopped: Missing input file.")
  } else {
    cat("File located successfully.\n")
  }
  
  cat("Reading FASTA file...\n")
  sequence_lines <- readLines(file_name, warn = FALSE)
  
  sequence <- paste(sequence_lines[!grepl("^>", sequence_lines)], collapse = "")
  sequence <- gsub("\\s+", "", sequence)
  
  # Convert to uppercase to handle case sensitivity
  sequence <- toupper(sequence)
  
  if (nchar(sequence) == 0) {
    stop("Sequence appears empty or incorrectly formatted.")
  } else {
    cat("Sequence successfully read.\n")
    cat("Sequence length:", nchar(sequence), "bases.\n")
  }
  
  sequence_length <- nchar(sequence)
  num_kb <- ceiling(sequence_length / 1000)
  
  cat("Dividing sequence into", num_kb, "kilobases...\n")
  
  if (num_kb < 1) {
    stop("Invalid sequence length. Cannot divide into kilobases.")
  }
  
  kb_counts <- list()
  
  cat("Processing nucleotide counts per kilobase...\n")
  
  for (i in 1:num_kb) {
    start <- (i - 1) * 1000 + 1
    end <- min(i * 1000, sequence_length)
    kb_seq <- substr(sequence, start, end)
    
    # --- QC: Check if segment extracted correctly ---
    if (nchar(kb_seq) == 0) {
      warning(paste("Empty segment detected at kilobase", i))
    }
    
    # Count nucleotides (case-insensitive already handled by toupper earlier)
    kb_counts[[paste0("kb_", i)]] <- list(
      A = stringr::str_count(kb_seq, "A"),
      C = stringr::str_count(kb_seq, "C"),
      G = stringr::str_count(kb_seq, "G"),
      T = stringr::str_count(kb_seq, "T")
    )
  }
  
  cat("Successfully generated list of nucleotide counts.\n")
  cat("Number of kilobases processed:", length(kb_counts), "\n")
  
  # --- Step 6: Print First 5 Lists for QC ---
  cat("\nPreview of first 5 kilobase nucleotide counts:\n")
  print(kb_counts[1:5])
  
  ##############################################
  # --- Part 4 ---
  ##############################################
  
  # --- Step 1: Load Required Packages ---
  if (!require(dplyr)) {
    install.packages("dplyr", dependencies = TRUE)
    library(dplyr)
  }
  
  ##############################################
  # --- Part 4.a: Create a Data Frame for the First 1000 bp ---
  ##############################################
  
  cat("Verifying that 'kb_counts' list exists...\n")
  
  if (!exists("kb_counts")) {
    stop("ERROR: The list 'kb_counts' does not exist. Please run Part 3 first.")
  } else if (!is.list(kb_counts) || length(kb_counts) == 0) {
    stop("ERROR: 'kb_counts' must be a non-empty list.")
  } else {
    cat("'kb_counts' is available and valid.\n")
  }
  
  # Select first 5 kilobases (or fewer if the list is smaller)
  num_to_show <- min(5, length(kb_counts))
  
  # Convert to a data frame for easy viewing
  df_first5 <- do.call(rbind, lapply(kb_counts[1:num_to_show], as.data.frame))
  df_first5 <- as.data.frame(df_first5)
  rownames(df_first5) <- paste0("kb_", 1:num_to_show)
  
  cat("\nFirst 5 kilobase rows from kb_counts:\n")
  print(df_first5)
  
  ##############################################
  # --- Part 4.b: Repeat for Each Kilobase ---
  ##############################################
  
  cat("Combining all kilobase data into one data frame...\n")
  
  # Convert list to data frame
  df_all_kb <- do.call(rbind, lapply(kb_counts, as.data.frame))
  df_all_kb <- as.data.frame(df_all_kb)
  rownames(df_all_kb) <- paste0("kb_", seq_len(nrow(df_all_kb)))
  
  # QC Check: Verify structure
  if (nrow(df_all_kb) == 0 || ncol(df_all_kb) != 4) {
    stop("Data frame not constructed correctly. Check list contents.")
  } else {
    cat("All kilobases successfully combined into data frame.\n")
    cat("Rows (kilobases):", nrow(df_all_kb), " | Columns (A, C, G, T):", ncol(df_all_kb), "\n")
  }
  
  cat("\nPreview of first 5 rows from combined data frame:\n")
  print(head(df_all_kb, 5))
  
  ##############################################
  # --- Part 4.c: Calculate Row Sums ---
  ##############################################
  
  cat("Calculating total nucleotide counts per kilobase...\n")
  
  # Add a new column for total count (sum of A, C, G, and T)
  df_all_kb$Sum <- rowSums(df_all_kb[, c("A", "C", "G", "T")])
  
  # Expected total = 1000 per kilobase (except possibly the last)
  expected_sum <- 1000
  invalid_rows <- which(df_all_kb$Sum != expected_sum)
  
  # QC: Check for kilobases with unexpected totals
  if (length(invalid_rows) == 0) {
    cat("All kilobases have expected total (1000 bases).\n")
  } else {
    cat("Some kilobases have unexpected totals:\n")
    print(invalid_rows)
  }
  
  cat("\nPreview of Final Data Frame (First 5 Rows):\n")
  print(head(df_all_kb, 5))
  
  cat("\nTotal Kilobases Processed:", nrow(df_all_kb), "\n")
  cat("Number of kilobases with exact 1000-base totals:", nrow(df_all_kb) - length(invalid_rows), "\n")
  cat("Number of kilobases with different totals:", length(invalid_rows), "\n")
  

    # What is the expected sum for each list?
    # ------------------------------------------------
    # Each kilobase represents 1000 DNA bases (A, C, G, and T combined).
    # Therefore, the expected sum for each list (row in df) = 1000.
     
    # Are there any lists whose sums are not equal to the expected value?
    # ------------------------------------------------
    # In most cases, every row will have a total of exactly 1000 bases.
    # However, if the final kilobase segment is shorter than 1000 bases
    # (because the total sequence length is not an exact multiple of 1000),
    # then the last list will have a smaller sum.
    # This is normal and expected.
     
    # General explanation for differences between expected and observed results:
    # ------------------------------------------------
    # Biological DNA sequences often have total lengths that are not perfectly
    # divisible by 1000. When we split the sequence into 1-kb chunks, the final
    # chunk might contain fewer than 1000 nucleotides.
    #
    # This causes the last row (or last few rows, depending on sequence length) to have a smaller total count of nucleotides.
    # 
    # Example:
    #   Suppose total sequence length = 3050 bases.
    #   Kilobase chunks → 1000, 1000, 1000, 50.
    #   The last chunk’s sum = 50 (less than expected 1000).
    #
    # Therefore, the small discrepancies are due to sequence length, not errors in reading or processing.
    ##############################################
