#Cargamos los paquetes necesarios
library(sangeranalyseR)
library(xlsx)
library(MLSTar)
library(Biostrings)
library(dplyr)
library(tools)
library(rBLAST)

#Seleccionar todos los archivos .ab1 del directorio
ab1files <- list.files(path = "C:/Users/jaime/Desktop/TFM_Legionella/secuencias",
                       full.names = TRUE, pattern = "\\.ab1$")

#Creamos una tabla con rutas, ids, gen, dirección de las secuencias
sample_names <- createSampleTable(dir_path = "C:/Users/jaime/Desktop/TFM_Legionella/secuencias")
head(sample_names)


#Creamos la función para convertir los archivos .ab1 a .fa
ab1tofasta()


#Añadimos una segunda columna con las rutas absolutas a los archivos fasta
sample_names["fa_file"] <- paste0(dirname(sample_names$ab1_file), "/fastas/",
                                   stringr::str_replace(basename(sample_names$ab1_file),
                                                        pattern = ".ab1",
                                                        replacement = ".fa"))
sample_names <- sample_names %>% relocate(fa_file, .after = 1)
head(sample_names)





#Reescribir nombres de los alelos
alelos_archivos <- list.files(path = "./references",
                      full.names = TRUE, pattern = "\\.txt$")


#Reescribir el archivo de alelos descritos para que pueda ser una entrada en makeblastdb
fastaBD(ewgli = "C:/Users/jaime/Desktop/TFM_Legionella/secuencias/references/asd ewgli.txt",
              gen = "asd",
              output_dir = "C:/Users/jaime/Desktop/TFM_Legionella/secuencias/references/references")

fastaBD(ewgli = "C:/Users/jaime/Desktop/TFM_Legionella/secuencias/references/references/fla ewgli.txt",
              gen = "fla",
              output_dir = "C:/Users/jaime/Desktop/TFM_Legionella/secuencias/references/references")

fastaBD(ewgli = "C:/Users/jaime/Desktop/TFM_Legionella/secuencias/references/references/mip ewgli.txt",
              gen = "mip",
              output_dir = "C:/Users/jaime/Desktop/TFM_Legionella/secuencias/references/references")

fastaBD(ewgli = "C:/Users/jaime/Desktop/TFM_Legionella/secuencias/references/references/momps ewgli.txt",
              gen = "momp",
              output_dir = "C:/Users/jaime/Desktop/TFM_Legionella/secuencias/references/references")

fastaBD(ewgli = "./references/neu ewgli.txt",
              gen = "neu",
              output_dir = "C:/Users/jaime/Desktop/TFM_Legionella/secuencias/references/references")


fastaBD(ewgli = "C:/Users/jaime/Desktop/TFM_Legionella/secuencias/references/references/pil ewgli.txt",
              gen = "pil",
              output_dir = "C:/Users/jaime/Desktop/TFM_Legionella/secuencias/references/references")

fastaBD(ewgli = "C:/Users/jaime/Desktop/TFM_Legionella/secuencias/references/references/pro ewgli.txt",
              gen = "pro",
              output_dir = "C:/Users/jaime/Desktop/TFM_Legionella/secuencias/references/references")





#Aplicar MLSTar version Windows a cada archivo .fasta creado

fafiles <- list.files(path = "./references",
                       full.names = TRUE, pattern = "\\.fas$")
fafiles

secuencias_20422 <- subset(sample_names, sample == "20422")$fa_file
scheme <- read.table("C:/Users/jaime/Desktop/TFM_Legionella/secuencias/scheme/st_profiles.txt",
                     header = TRUE)
head(scheme)

#Juntar los archivos por directorio
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

sample_multifasta(sample_table = sample_names, out_dir = "./samples_multifasta")

#Hacer gen a gen

SBT <- function(files, db , scheme, dir_path = "."){
  MLST <-doMLSTw(infiles = files,
          org = 'test',
          schemeFastas = db,
          write = "new",
          schemeProfile = scheme)
  MLST$result
  c <- data.frame(MLST$result)
  write.csv(MLST$result, file = file.path(dir_path, "MLSTar_results.csv"))

    }
  }
  return(MLST)
}

unlink("./alleles_test_1/", recursive = TRUE)


files = list.files("C:/Users/jaime/Desktop/TFM_Legionella/secuencias/fastas",
                   full.names = TRUE, pattern = "\\.fa$")
fafiles <- list.files("C:/Users/jaime/Desktop/TFM_Legionella/secuencias/references/",
                      full.names = TRUE, pattern = "\\.fas$")

sample_names
muestra <- "./fastas/AF4437+20412mipF+M13.fa"
muestra

SBT(files = files,
    db = fafiles,
    scheme = "C:/Users/jaime/Desktop/TFM_Legionella/secuencias/scheme/st_profiles.txt")




library(rBLAST)
route <- normalizePath(".")
route
db_mip <- "./references/mip"
db_mip
seq <- readDNAStringSet("./samples_multifasta/20412.fasta")
c1 <- predict(db_mip, seq)
eval <- 0.01

cmd<-paste0("blastn -word_size 50 -ungapped -dust no -query ",muestra,
            " -db ",db_mip,
            " -outfmt 1",
            " -max_target_seqs 1",
            " -evalue ",eval)

cmd
cmd
system(cmd)





SBT <- function(files, db , scheme, dir_path = "."){
  MLST <-doMLSTw(infiles = files,
                 org = 'test',
                 schemeFastas = db,
                 write = "all",
                 schemeProfile = scheme)
  mlstresu <- MLST$result
  write.csv(MLST$result, file = file.path(dir_path, "MLSTar_results.csv"))
  return(MLST)
}
unlink("./alleles_test_1/", recursive = TRUE)

mlstresu <- SBT(files = files,
                db = fafiles,
                scheme = "C:/Users/jaime/Desktop/TFM_Legionella/secuencias/scheme/st_profiles.txt")

resu <- mlstresu$result
resu

ibrary(dplyr)
library(stringr)
library(tibble)

# if rownames contain the IDs:
df <- resu %>% rownames_to_column("id")

# Extract the 5-digit number (e.g. 20422 or 20421)
df <- df %>%
  mutate(
    group = str_extract(id, "\\d{5}"),
    gen = {
      genes <- c("fla", "pil", "asd", "mip", "momp", "pro", "neu")
      regex_gen <- paste(genes, collapse = "|")
      g <- regmatches(id, regexpr(regex_gen, id, perl = TRUE))
      g <- ifelse(length(g) == 0 | g == "", NA, g)
      g
    }
  )%>%
  ungroup()

df

# Merge rows by that number
merged <- df %>%
  group_by(group) %>%
  summarise(
    across(
      all_of(setdiff(names(df), c("id","group"))),
      ~ {
        vals <- na.omit(.x)
        if(length(vals) == 0) "No encontrado" else first(vals)
      },
      .names = "{.col}"
    ),
    .groups = "drop"
  )
merged


#
createSchemes <- function(dir_path = ".", ref_genes){
  out_dir <- "./schemes"
  ref_dir <- paste0(dir_path,"references/")
  print(ref_dir)
  if (!dir.exists(out_dir)) {
    dir.create(path = out_dir, recursive = FALSE, showWarnings = FALSE)
    cat("Directorio", out_dir, "creado", "\n")
  }
  else{
    cat("Directorio", out_dir, "ya existe, continuando...", "\n")
  }

  for (i in ref_genes){
    namefas <- paste0(i,".fas")
    print(namefas)
    ref_route<- list.files(path = ref_dir, full.names = TRUE, pattern = namefas)
    print(ref_route)
    fasta <- readDNAStringSet(ref_route)
    head_gene <- names(fasta)
    num_gene <- str_split_i(head_gene, "_", 2)
    scheme <- data.frame(ST = num_gene)
    scheme[[i]] <- num_gene
    print(scheme)
    write.table(scheme, file = paste0(i, "_profile"))
  }
}

createSchemes(dir_path = "C:/Users/jaime/Desktop/TFM_Legionella/secuencias/", "asd")

fd <- readDNAStringSet("C:/Users/jaime/Desktop/TFM_Legionella/secuencias/references/asd.fas")
PolyPeakParser()





if (length(sbt_results) > 0) {
  final_results_df <- do.call(rbind, sbt_results)
  sbt_results(final_results_df)
} else {
  showNotification("Ningún análisis SBT pudo completarse con éxito.", type = "error")
  sbt_results(NULL)
}


library(rBLAST)
pro_db <- blast(db = "C:/Users/jaime/Desktop/TFM_Legionella/secuencias/references/pil", type = "blastn")
db <- "C:/Users/jaime/Desktop/TFM_Legionella/secuencias/references/pil"
seq <- "C:/Users/jaime/Desktop/TFM_Legionella/secuencias/fastas/AJ6519+20410pilF+M13F.fa"
predict(pro_db, seq)

fafiles <- list.files("C:/Users/jaime/Desktop/TFM_Legionella/secuencias/references/",
                      full.names = TRUE, pattern = "\\pro.fas")
pairwiseAlignment(pattern = db, subject = seq, type = "local")
fafiles

comando <- paste(
  "blastn -query",
  seq,
  "-db",
  db,
  " -outfmt  \"0 pident gaps \"",
  "-max_target_seqs 1 "
)
comando
system(comando)
Subject <- seq_along(db)
pairs = data.frame(Pattern = 1, Subject)
pairs

aligned <- AlignPairs(pattern=seq,
           subject=db,
           pairs= pairs,
           type="sequences",
           processors=NULL)
db$
head(aligned[[3]], n = 1)
library(seq2R)
ruta <- "C:/Users/jaime/Desktop/TFM_Legionella/secuencias/AF4429+20410proF+M13.ab1"
dat_ab <- read.abif(ruta)
plotabif(dat_ab)
vmatchPattern(db, seq)



db <- "C:/Users/jaime/Desktop/TFM_Legionella/secuencias/references/mip"
seq <- "C:/Users/jaime/Desktop/TFM_Legionella/secuencias/fastas/AF4437+20412mipF+M13.fa"

comando <- paste(
  "blastn -query",
  seq,
  "-db",
  db,
  " -outfmt  0",
  "-max_target_seqs 1 "
)

system(comando)
doMLST()
a1 <- readsangerseq("C:/Users/jaime/Desktop/TFM_Legionella/secuencias/AF4441+20400proR+M13.ab1")
a1
chromatogram(b1, width = 500,
             height = NA,
)
b1 <- makeBaseCalls(a1, ratio = 0.33)

seq1 <- "C:/Users/jaime/Desktop/TFM_Legionella/secuencias/AJ6519+20410pilF+M13F.ab1"

sread <- SangerRead(readFeature           = "Forward Read",
                    readFileName          = seq1,
                    geneticCode           = GENETIC_CODE,
                    TrimmingMethod        = "M1",
                    M1TrimmingCutoff      = 1,
                    M2CutoffQualityScore  = NULL,
                    M2SlidingWindowSize   = NULL,
                    baseNumPerRow         = 100,
                    heightPerRow          = 200,
                    signalRatioCutoff     = 0.2,
                    showTrimmed           = TRUE)
qualityBasePlot(sread)

seq2 <- readsangerseq("C:/Users/jaime/Desktop/TFM_Legionella/secuencias/AF4437+20412mipF+M13.ab1")
seq21 <- makeBaseCalls(seq2)
seq21
PolyPeakParser()
