library(sangeranalyseR)
library(xlsx)
library(MLSTar)
library(Biostrings)
library(dplyr)
library(stringr)
library(tibble)

createSampleTable <- function(dir_path = ".",
                           ab1files = list.files(path = "C:/Users/jaime/Desktop/TFM_Legionella/secuencias",
                                                 full.names = TRUE,
                                                 pattern = "\\.ab1$") ){
  genes <- c("fla", "pil", "asd", "mip", "momp", "pro", "neu")
  ab1names <- basename(ab1files)
  ab1route <- normalizePath(ab1files, winslash = "/")

  sample_table <- data.frame(
    ab1_file = character(),
    id_serv = character(),
    sample = character(),
    genes = character(),
    dir = character(),
    stringsAsFactors = FALSE
  )

  for (i in ab1names){
    prefix <- paste0("(.+?)(?=", paste(genes, collapse = "|"), ")")
    id_sample <- regmatches(i, regexpr(prefix, i, perl = TRUE))
    if (length(id_sample) == 0) id_sample <- NA
    id_sample <- gsub("[-_]+$", "", id_sample)
    id_sample <- ifelse(id_sample == "", NA, id_sample)

    if (!is.na(id_sample) && grepl("\\+", id_sample)) {
      id_serv <- sub("\\+.*", "", id_sample)
      sample <- sub(".*\\+", "", id_sample)
    } else {
      id_serv <- id_sample
      sample <- NA
    }

    regex_gen <- paste(genes, collapse = "|")
    gen <- regmatches(i, regexpr(regex_gen, i, perl = TRUE))
    gen <- ifelse(gen == "", NA, gen)
    if (length(gen) == 0) gen <- NA


    dir_match <- regmatches(i, regexpr(paste0("(", regex_gen, ")([FR])"), i, perl = TRUE))
    if (length(dir_match) == 0) {
      dir_match <- regmatches(i, regexpr("M13([FR])", i, perl = TRUE))
    }
    if (length(dir_match) == 0) {
      dir_val <- NA
    } else {
      dir_val <- sub(".*(F|R).*", "\\1", dir_match)
    }

    idx <- which(ab1names == i)

    sample_table <- rbind(
      sample_table,
      data.frame(
        ab1_file = ab1route[idx],
        id_serv = id_serv,
        sample = sample,
        genes = gen,
        dir = dir_val,
        stringsAsFactors = FALSE
      )
    )
  }


  output_name <- ("sample_table_ab1.xlsx")
  output_path <- file.path(dir_path, output_name)

  write.xlsx(sample_table, output_path, col.names = TRUE, row.names = FALSE)

  return(sample_table)
}



ab1tofasta <- function(abDir = ".",
                       newDir = "./fastas",
                       table_samples = sample_names){
  ab1 = list.files(path = "abDir", full.names = FALSE, pattern = "\\.ab1$")

  if (!dir.exists(newDir)) {
    dir.create(path = newDir, recursive = FALSE, showWarnings = FALSE)
    cat("Directorio", newDir, "creado", "\n")
  }
  else{
    cat("Directorio", newDir, "ya existe, continuando...", "\n")
  }

  for (i in 1:nrow(table_samples)){
    sanger_read <- ""
    orien <- ""
    file <- table_samples$ab1_file[i]

    dir_col <- table_samples$dir[i]

    if (dir_col == "F"){
      orien <- "Forward Read"}
    else if (dir_col == "R"){
      orien <- "Reverse Read"}
    else{
      cat("ATENCIÓN: Archivo", basename(ab1), "no tiene marcador F+ o R+. Saltando...\n")
    }

    sanger_read <- SangerRead(readFeature = orien,
                              readFileName = file,
                              geneticCode = GENETIC_CODE,
                              TrimmingMethod = "M1",
                              M1TrimmingCutoff = 0.0001,
                              showTrimmed = TRUE)
    fa <- writeFasta(sanger_read,
                     outputDir = newDir,
                     compress = FALSE)
    cat("Archivo FASTA creado:", basename(fa), "\n")
  }
  num_ab1 <- length(list.files(path = abDir, pattern = "\\.ab1$", full.names = TRUE))
  num_fa <- length(list.files(path = newDir, pattern = "\\.fa$", full.names = TRUE))

  if (num_ab1 == num_fa) {
    message("El número de archivos .ab1 (", num_ab1, ") coincide con el número de archivos .fa (", num_fa, ").")
  } else {
    stop("El número de archivos .ab1 (", num_ab1, ") no coincide con el número de archivos .fa (", num_fa, "). Verifica los directorios.")
  }

  return(paste("Las", num_fa, "secuencias .fa están en el directorio", newDir))
}


sample_multifasta <- function(sample_table, col_sample = "sample", col_route= "fa_file", dir_path = ".", out_dir = "."){
  if (!dir.exists(out_dir)) {
    dir.create(path = out_dir, recursive = FALSE, showWarnings = FALSE)
    cat("Directorio", out_dir, "creado", "\n")
  }
  else{
    cat("Directorio", out_dir, "ya existe, continuando...", "\n")
  }
  grouped_samples <- sample_table %>%
    group_by(!!sym(col_sample)) %>%
    summarise(
      fa_file = list(!!sym(col_route)),
      .groups = "drop"
    )

  for (i in 1:nrow(grouped_samples)){
    id <- grouped_samples[[col_sample]][i]
    routes <- unlist(grouped_samples$fa_file[i])
    seq_sample <- list()

    for (route in routes){
      if (file.exists(route)){
        seq_1 <- readDNAStringSet(filepath = route)
        seq_sample <- c(seq_sample, seq_1)
      } else {
        warning("Archivo no encontrado en la ruta", route)
      }
    }

    if (length(seq_sample) > 0 ) {
      all_seq <- do.call(c, seq_sample)

      out_file <- paste0(out_dir,"/",id,".fasta")
      writeXStringSet(
        x = all_seq,
        filepath = out_file,
        format = "fasta"
      )
      message(paste("Archivo creado con éxito para", id, "en", out_file))
    } else {
      warning(paste("No se encontraron secuencias válidas para la muestra:", id))
    }
  }
}


fastaBD <- function(ewgli, gen, output_dir = "."){

  sequ <- readDNAStringSet(ewgli)
  new_name <- paste0(gen, "_", names(sequ))
  names(sequ) <- new_name
  output_filename <- paste0(gen, ".fas")
  output_path <- file.path(output_dir, output_filename)
  writeXStringSet(sequ, output_path, format = "fasta")

  file_seq <- paste0(output_dir, "/", gen,".fas")
  lines <- readLines(file_seq)
  lines_clean <- gsub("\t", "\n", lines)
  writeLines(lines_clean, file_seq)


  cat(output_filename, "guardado en", output_path, "\n")
}

SBT <- function(files, db , scheme, dir_path = "."){
  MLST <-doMLSTw(infiles = files,
                 org = 'test',
                 schemeFastas = db,
                 write = "new",
                 schemeProfile = scheme)
  MLST$result
  write.csv(MLST$result, file = file.path(dir_path, "MLSTar_results.csv"))
  return(MLST)
}

