

#install.packages("shinydashboard")
#install.packages("shinyBS")
#install.packages("shiny.semantic")

library(shiny)
library(shinydashboard)
library(seqvisr)
library(bslib)
library(shinyBS)
library(shinyjs)
library(shinyMobile)
library(shiny.semantic)

ui <- function() {
  semanticPage(
    title = "RNase Sequence Visualization (RNaseViz)",
    h1("RNase Sequence Visualization (RNaseViz)"),
    div(
      class = "ui stackable grid",  # stackable grid for responsiveness
      div(
        class = "six wide column",  # This defines the width of the sidebar
        div(
          selectInput("rnase", "Choose RNase:",
                      choices = c("RNASEH1", "RNASEH2A", "RNASEH2B", "RNASEH2C", "AGO2", "DICER1", "ELAC2", "DIS3L2", "RNASET2", "PARN", "PRORP")),
          actionButton("visualizeBtn", "Visualize", class = "large primary active button",style = "font-size: 18px;"),
          numericInput("mutationPos", "Enter position in human reference sequence:", 1, min = 1),
          actionButton("highlightMutation", "Highlight position", class = "large secondary button",style = "font-size: 18px;"),
          textOutput("alignedPosInfo"),
          uiOutput("rnaseInfo"),
          uiOutput("rnaseImageContainer"),
          style = "padding-bottom: 20px;"
        ),
        style = "overflow: auto; max-height: 100vh;"
      ),
      div(
        class = "ten wide column",  # This defines the width of the main panel
        uiOutput("mainPanelSegments")
      )
    ),
    style = "height: 100vh; padding: 1em; text-align: justify;",  # Ensures the page takes the full height
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
          fint-size: 16px;
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
        
      "))
    )
  )
}
#max-width: 100%;
#max-height: calc(600px - 20px);
