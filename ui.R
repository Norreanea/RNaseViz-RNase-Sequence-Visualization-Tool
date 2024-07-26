

#install.packages("shinydashboard")
#install.packages("shinyBS")
#install.packages("shiny.semantic")
#Load libraries
library(shiny)
library(shinydashboard)
library(seqvisr)
library(bslib)
library(shinyBS)
library(shinyjs)
library(shinyMobile)
library(shiny.semantic)

#name == "Saccharomyces_cerevisiae"  ~ 'SC',
#       name == "Drosophila_melanogaster"  ~ 'DM',
#       name == "Caenorhabditis_elegans"  ~ 'CE',
#       name == "Danio_rerio"  ~ 'DR',
#       name == "Xenopus_tropicalis"  ~ 'XT',
#       name == "Homo_sapiens"  ~ 'HS',
#       name == "Rattus_norvegicus"  ~ 'RT',
#       name == "Mus_musculus"  ~ 'MM',
#       TRUE ~ 'Other'),
ui <- function() {
  semanticPage(
    title = "RNase Sequence Visualization (RNaseViz)",
    h1("RNase Sequence Visualization (RNaseViz)"),
    div(
      class = "ui stackable grid",  # stackable grid for responsiveness
      div(
        class = "six wide column",  # defines the width of the sidebar
        div(
          selectInput("rnase", "Choose RNase:",
                      choices = c("RNASEH2A", "RNASEH2B", "RNASEH2C", "AGO2", "RNASET2", "DICER1", "DIS3L2", "ELAC2", "PRORP", "PARN", "RNASEH1")),
          actionButton("visualizeBtn", "Visualize", class = "large primary active button", style = "font-size: 18px;"),
          
          # Separated elements with grey background
          div(
            class = "ui segment",
            style = "background-color: #f0f0f0; padding: 15px; border-radius: 5px;",
            h3(" "),
            # selectInput("referenceOrganism", "Select organism:", 
            #             choices = c("Homo sapiens", "Mus musculus", "Rattus norvegicus", "Danio rerio", "Xenopus tropicalis",
            #                         "Caenorhabditis elegans", "Drosophila melanogaster","Saccharomyces cerevisiae")),
            selectInput("referenceOrganism", "Select organism:", 
                        choices = c("Homo sapiens", "Mus musculus", "Rattus norvegicus", "Danio rerio", "Xenopus tropicalis",
                                    "Caenorhabditis elegans", "Drosophila melanogaster","Saccharomyces cerevisiae"), 
                        selected = "Homo sapiens"),
            h3(" "),
            numericInput("mutationPos", "Enter position in reference sequence for selected organism:", 1, min = 1),
            actionButton("highlightMutation", "Highlight position", class = "large secondary button", style = "font-size: 18px;"),
            htmlOutput("alignedPosInfo")
          ),
          
          uiOutput("rnaseInfo"),
          uiOutput("rnaseImageContainer"),
          uiOutput("BioRender"),
          style = "padding-bottom: 20px;"
        ),
        style = "overflow: auto; max-height: 100vh;"
      ),
      div(
        class = "ten wide column",  # defines the width of the main panel
        uiOutput("mainPanelSegments"),
        style = "overflow: auto; max-height: 100vh;"
      )
    ),
    style = "height: 100vh; padding: 1em; text-align: justify;",  # ensures the page takes the full height
    tags$head(
      tags$style(HTML("
      /* Layout */
        .ui.stackable.grid {
          display: flex;
          
          flex-direction: column;
          height: 100%;
        }
      /* Buttons */
        ui.button, .ui.buttons .button {
          font-size: 20px !important; 
          font-weight: bold !important; 
        }
        /* Labels in forms */
        div[class*='field'] > label {
          font-size: 20px;
          font-weight: bold;
        }
        /* Labels UI frameworks */
        .ui.form .field > label {
          font-size: 20px;  
          font-weight: bold;  
        }

        /* rnaseInfo */
        #rnaseInfo {
          font-size: 18px; 
        }
        /* alignedPosInfo */
        #alignedPosInfo {
          font-size: 16px;
        }

        /* main title */
        h1, .header {
          text-align: center;
          font-size: 24px; 
        }

        /* rnaseImageContainer */
        #rnaseImageContainer {
          max-width: auto;
          overflow: hidden;
          display: flex;
          align-items: center; /* Centers the image vertically */
          justify-content: center; /* Centers the image horizontally */
        }

        /* rnaseImage */
        #rnaseImage {
         max-width: 100%;
          max-height: calc(600px - 20px);
          height: auto; /* Maintain aspect ratio */
          width: auto; /* Adjust width based on height */
        }

        /* Custom grey background for segment */
        .grey-background {
          background-color: #f0f0f0;
          padding: 15px;
          border-radius: 5px;
        }
        
        
      "))
    )
  )
}


#max-width: 100%;
#max-height: calc(600px - 20px);
