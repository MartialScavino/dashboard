if (!require("rstudioapi")) install.packages("rstudioapi")

# Setting working directory
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(shinydashboard)
library(shiny)
library(shinysky)
library(shinyWidgets)
library(shinyjs)
library(shinyFiles)
library(Seurat)
library(ggplot2)
library(viridis)
library(cowplot)
library(DT)
library(babelgene)
library(UCell)
library(stringr)

runApp('.')
