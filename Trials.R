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



###### Seleccionar por muestra o por gen
levels(as.factor(sample_names$gen))
subset(sample_names, sample == "20450")
subset(sample_names, sample %in% c("20450", "20421"))
subset(sample_names, sample == "20450" & gen == "pro")$fa_file


#Reescribir nombres de los alelos
alelos_archivos <- list.files(path = "./references",
                      full.names = TRUE, pattern = "\\.txt$")


#Reescribir el archivo de alelos descritos para que pueda ser una entrada en makeblastdb
fastaBD(ewgli = "./references/asd ewgli.txt",
              gen = "asd",
              output_dir = "./references")

fastaBD(ewgli = "./references/fla ewgli.txt",
              gen = "fla",
              output_dir = "./references")

fastaBD(ewgli = "./references/mip ewgli.txt",
              gen = "mip",
              output_dir = "./references")

fastaBD(ewgli = "./references/momps ewgli.txt",
              gen = "momp",
              output_dir = "./references")

fastaBD(ewgli = "./references/neu ewgli.txt",
              gen = "neu",
              output_dir = "./references")

fastaBD(ewgli = "./references/neuAH ewgli.txt",
              gen = "neuAH",
              output_dir = "./references")

fastaBD(ewgli = "./references/pil ewgli.txt",
              gen = "pil",
              output_dir = "./references")

fastaBD(ewgli = "./references/pro ewgli.txt",
              gen = "pro",
              output_dir = "./references")





#Aplicar MLSTar version Windows a cada archivo .fasta creado

fafiles <- list.files(path = "./references",
                       full.names = TRUE, pattern = "\\.fas$")
fafiles

secuencias_20422 <- subset(sample_names, sample == "20422")$fa_file
secuencias_20410
scheme <- read.table("C:/Users/jaime/Desktop/TFM_Legionella/secuencias/scheme/st_profiles.txt",
                     header = TRUE)
head(scheme)

#Hacer gen a gen
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

SBT(files = secuencias_20422,
    db = fafiles,
    scheme = "C:/Users/jaime/Desktop/TFM_Legionella/secuencias/scheme/st_profiles.txt")

doMLSTw(infiles = secuencias_20422,
        org = 'test',
        schemeFastas = fafiles,
        write = "new",
        schemeProfile = "C:/Users/jaime/Desktop/TFM_Legionella/secuencias/scheme/st_profiles.txt")

tipificacion(seq = "./fastas/AF4440+20400flaR+M13.fa",
                 db = "./references/fla.fas")

makeblastdb("C:/Users/jaime/Desktop/TFM_Legionella/secuencias/neu ewgli.txt")

db_files <- list.files(path = ".", full.names = TRUE, pattern = "\\.txt$")
db_files
db <- makeblastdb("./references/fla.fas")

