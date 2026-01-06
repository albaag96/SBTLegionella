
library(shiny)
library(ggplot2)
library(plotly)
library(bslib)
library(shinyFiles)

fluidPage(
  titlePanel(
    img(src = "/bitmap.png", height = "100px")
    ),

  sidebarLayout(
    sidebarPanel(
      h5(strong("Create Table of Samples:")),
      shinyDirButton("dir", "Select Directory", "Select Directory:"),
      verbatimTextOutput("path"),

      actionButton("do_table", "Create Sample Table"),

      h5(strong("Convert ab1 samples to FASTA:")),
      actionButton("convert_fasta", "Convert to FASTA"),
    ),

    mainPanel(
      tabsetPanel(
        tabPanel("Available data",
                 fluidRow(
                   column(6,
                 pickerInput("sample", "Select sample:",
                             choices = c("All"))),
                 column(6,
                 pickerInput("gene", "Select gene:",
                             choices = c("All"),
                             selected = c("All"),
                             multiple = TRUE))
                 ),
                 DTOutput("table"),
                 actionButton("do_sbt", "SBT"),
                 dataTableOutput("sbt_table"),
                 downloadButton("download_sbt", "Download SBT")
                 ),
        tabPanel("Chromatograph",
                 actionButton("do_chromatograph", "Chromatograph"),
                 plotOutput("plot_chrom", height = "600px", width = "2000px")),
        tabPanel(
          "Alignment",
          actionButton("do_alignment", "Alignment"),
          htmlOutput("plot_alignment")),
        tabPanel("About",
                 HTML("
    <div style='font-family: Arial, sans-serif; line-height: 1.5;'>

      <h2 style='font-weight: bold; color: #2C3E50;'>Instruction Manual</h2>
      <ol>
        <li>Click 'Select Directory' to open a file browser, choose the folder containing the '.ab1' sequences, and then click 'Select'.</li>
        <li>Click 'Create Sample Table' to display the files and generate an Excel file.</li>
        <li>Click 'Convert to FASTA' to convert ABIF files into FASTA files in a new 'fastas' folder.</li>
        <li>Select the file(s) to type and click 'SBT'. You can filter by sample and gene.</li>
        <li>A table will appear below showing samples by row and genes by column with assigned alleles. You can download it by clicking 'Download SBT'.</li>
        <li>The table shows the allele number if the sequence matches a known allele exactly. Sequences with ≥90% identity and ≥90% coverage to a known allele are labeled as 'u+number'. If no allele meets these criteria, the cell is left blank.</li>
        <li>If an allele could not be assigned, potential base-calling errors can be checked by selecting the sequence and clicking 'Chromatograph' in the 'Chromatograph' tab.</li>
        <li>You can also view the alignment with the best-matching allele by selecting the sequence and clicking 'Alignment' in the 'Alignment' tab.</li>
      </ol>

      <hr style='border:1px solid #ccc;'>

      <h2 style='font-weight: bold; color: #2C3E50;'>About this Application</h2>
      <p>
        This application was developed as a Master’s thesis project in collaboration with the
        Clinical and Environmental Infectious Disease Study Group (CEID) at the Germans Trias i Pujol Research Institute (IGTP),
        by Alba Arranz García, José Francisco Sánchez Herrero and Noemí Párraga Niño.
      </p>

      <hr style='border:1px solid #ccc;'>

      <h2 style='font-weight: bold; color: #2C3E50;'>Acknowledgments</h2>
      <p>
        This application uses the <strong>sangeranalyseR</strong> and <strong>MLSTar</strong> packages. We thank the authors for making these tools available:
      </p>
      <ul>
        <li><strong>MLSTar:</strong> Ferrés I, Iraola G. MLSTar: automatic multilocus sequence typing of bacterial genomes in R. PeerJ. 2018;6
        <a href='https://doi.org/10.7717/peerj.5098' target='_blank'>https://doi.org/10.7717/peerj.5098</a>.</li>
        <li><strong>sangeranalyseR:</strong> Chao K, Barton K, Palmer S, Lanfear R (2021). “sangeranalyseR: simple and interactive analysis of Sanger sequencing data in R.” Genome Biology and Evolution.
        <a href='https://doi.org/10.1093/gbe/evab028' target='_blank'>https://doi.org/10.1093/gbe/evab028</a>.</li>
      </ul>

    </div>
  "))
      )
    )
  )
)


