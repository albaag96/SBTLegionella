
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
                 actionButton("write_sbt", "Save Results")
                 ),
        tabPanel("Chromatograph",
                 plotOutput("plot_chrom"),
                 actionButton("do_chromatograph", "Chromatograph")),
        tabPanel("Alignment",
                 htmlOutput("plot_alignment"))
      )
    )
  )
)


