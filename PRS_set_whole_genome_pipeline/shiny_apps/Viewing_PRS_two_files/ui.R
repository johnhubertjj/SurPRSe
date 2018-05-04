#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(scales)
library(data.table)
library(plyr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(DT)
library(shiny)
library(stringr)
library(grid)
# Define UI for application that draws a histogram
shinyUI( fluidPage(
          # Sidebar layout with input and output definitions ----
          
          headerPanel("PRS gene-set analysis viewer"),
          
          sidebarLayout(
            
            ## Write shiny UI across 4 parameters in the data table
            sidebarPanel(
              fileInput("file1", "Choose an input file",
                        multiple = F),
              radioButtons("step_regress", label = "Add stepwise regression for DSM grouping?",
                            c("Yes", "No"),
                           selected = "No"),
              uiOutput("Significance_threshold"
                       ),
              checkboxGroupInput("Gene_regions", label = "Length of Gene regions:",
                                 choices = c("extended","normal","Full")),
              uiOutput("DSM"),
              uiOutput("geneset"),
              sliderInput('plotHeight', 'Bar which does nothing, use if bored', 
                          min = 100, max = 2000, value = 1000)
            ),
            
            mainPanel(
              tabsetPanel(id = "tabs",
                tabPanel("Plots",
                        conditionalPanel(
                          condition = "input.tabs == 'Plots' & input.step_regress == 'No'",
                          plotOutput('PvalPlot'),
                         plotOutput('Beta_plot'),
                         plotOutput('R2_plot')),
                        conditionalPanel(condition = "input.tabs == 'Plots' & input.step_regress == 'Yes'",
                                         plotOutput('PvalPlot_step'),
                                         plotOutput('Beta_plot_step'),
                                         plotOutput('R2_plot_step'))),
                tabPanel("Table", dataTableOutput("summary_table"))
              )
            )
          )
        )
) 
