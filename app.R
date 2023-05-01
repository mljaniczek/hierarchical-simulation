#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)
library(lme4)
library(lmerTest)
library(gtsummary)

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Multi-level Data (Families within Regions)"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            numericInput("nregions", "Number of Regions", 2),
            numericInput("randintvar", "Region Random Intercept Variance", 3),
            numericInput("randintvar2", "Family Random Intercept Variance", 3),
            numericInput("intercept", "Overall Intercept", 0),
            numericInput("groupeffect", "Variable Effect", 1),
            numericInput("residvar", "Fixed variance", 1),
            numericInput("npatients", "Number of subjects per Region", 50),
            numericInput("nfamilies", "Number of Families per Region", 2)
        ),

        # Show a plot of the generated distribution
        mainPanel(
          tabsetPanel(type = "tabs",
                      # h4("Male and Female Genotype Percentage per Generation"),
                      tabPanel("Plots", plotOutput("distPlot1"), plotOutput("distPlot2"))
        )
    )
))

# Define server logic required to draw a histogram
server <- function(input, output) {
    output$distPlot1 <- renderPlot({
      set.seed(1219)
      randint_reg = rnorm(input$nregions, sd = sqrt(input$randintvar))
      
      # make random size of patients per region - maybe do this later
      # minpatients = 20
      # maxpatients = 50
      # npatients = as.integer(minpatients + ceiling(runif(n)*maxpatients-minpatients))
      
      # keep same number of patients per region for now
      npatients = rep(input$npatients, input$nregions)
      region = seq(1,input$nregions)
      nregions2 = input$nregions
      sampsize <- sum(npatients)
      data1 <- as.data.frame(cbind(region, randint_reg, npatients))
      data2 <- slice(data1, rep(seq(nregions2), npatients))
      # make groups of correlated data within region called "families"
      # keeping just 2 families for now for simplicity
      nfamilies = input$nfamilies
      # make random intercept for each family
      randint2 = data.frame(
        familyid = c(1:(nregions2*nfamilies)),
        randintfam = rnorm(nregions2*nfamilies, sd = sqrt(input$randintvar2)))
      
      
      # now collect the data together
      simple <- data2 %>%
        mutate(familyid = rep(seq(nregions2*nfamilies), c(rep(sampsize/(nregions2*nfamilies),nregions2*nfamilies)))) %>%
        left_join(randint2, by = "familyid") %>%
        mutate(familyid = as.factor(familyid))
      finaldat <- simple %>%
            # make normal predictor
            mutate(x1 = (rnorm(sampsize)),
                   # make binary demographic variable (male/female)
                   female = rep(c(0,1), times = sampsize/2)) %>%
            mutate(
              # now make linear predictor model
              z = randint_reg + randintfam + x1 + input$groupeffect*female + (rnorm(sampsize, sd = sqrt(input$residvar))), # add random effects, fixed effects, and random noise
              # trying out 2 ways to make binary outcome
              pr1 = exp(z)/(1+exp(z)),
              runis = runif(sampsize, 0, 1),
              y1 = ifelse(runis < pr1, 1, 0)) # pulled from Ken's blog
      ggplot(finaldat, aes(x=as.factor(y1), y=x1, color = familyid)) +
        geom_boxplot()+
        #ylim(0, NA) +
        theme_minimal()
    })
    
    set.seed(1219)
    # pulled from Ken's blog
    
    randint_reg = reactive({rnorm(input$nregions, sd = sqrt(input$randintvar))})
    
    # make random size of patients per region - maybe do this later
    # minpatients = 20
    # maxpatients = 50
    # npatients = as.integer(minpatients + ceiling(runif(n)*maxpatients-minpatients))
    
    # keep same number of patients per region for now
    npatients = reactive({rep(input$npatients, input$nregions)})
    region = reactive({seq(1,input$nregions)})
    sampsize <- reactive({sum(npatients())})
    data1 <- reactive({as.data.frame(cbind(region(), randint_reg(), npatients()))})
    data2 <- reactive({slice(data1(), rep(seq(input$nregions), npatients()))})
    # make groups of correlated data within region called "families"
    # keeping just 2 families for now for simplicity
    # make random intercept for each family
    randint2 = reactive({data.frame(
      familyid = c(1:(input$nregions*input$nfamilies)),
      randintfam = rnorm(input$nregions*input$nfamilies, sd = sqrt(input$randintvar2)))})
    
    
    # now collect the data together
    simple <- reactive({data2() %>%
        mutate(familyid = rep(seq(input$nregions*input$nfamilies), c(rep(sampsize()/(input$nregions*input$nfamilies),input$nregions*input$nfamilies)))) %>%
        left_join(randint2(), by = "familyid") %>%
        mutate(familyid = as.factor(familyid))})
    finaldat <- reactive({simple() %>%
        # make normal predictor
        mutate(x1 = (rnorm(sampsize())),
               # make binary demographic variable (male/female)
               female = rep(c(0,1), times = sampsize()/2)) %>%
        mutate(
          # now make linear predictor model
          z = randint_reg() + randintfam + x1 + input$groupeffect*female + (rnorm(sampsize(), sd = sqrt(input$residvar))), # add random effects, fixed effects, and random noise
          # trying out 2 ways to make binary outcome
          pr1 = exp(z)/(1+exp(z)),
          runis = runif(sampsize(), 0, 1),
          y1 = ifelse(runis < pr1, 1, 0))
    })
    
    output$distPlot2 <- renderPlot({
      print(ggplot(finaldat(), aes(x=as.factor(y1), y=x1, color = familyid)) +
        geom_boxplot()+
        #ylim(0, NA) +
        theme_minimal())
    })
}

# Run the application 
shinyApp(ui = ui, server = server)



#old
# 
# # make random intercept for each region
# randint_reg = reactive({norm(input$nregions, sd = sqrt(input$randintvar))})

# # keep same number of patients per region for now
# npatients = reactive({rep(input$npatients, input$nregions)})
# region = reactive({seq(1,input$nregions)})
# nregions2 = reactive({input$nregions})
# randintvar2 = reactive({input$randintvar2})
# sampsize = sum(npatients())
# data1 = as.data.frame(cbind(region(), randint_reg(), npatients()))
# data2 = slice(data1, rep(seq(nregions2()), npatients()))
# # make groups of correlated data within region called "families"
# # keeping just 2 families for now for simplicity
# nfamilies = reactive({input$nfamilies})
# # make random intercept for each family
# randint2 = data.frame(
#   familyid = c(1:(nregions2()*nfamilies())),
#   randintfam = rnorm(nregions2()*nfamilies(), sd = sqrt(randintvar2())))
# 
# 
# # now collect the data together
# simple = data2 %>%
#   mutate(familyid = rep(seq(nregions2()*nfamilies()), c(rep(sampsize()/(nregions2()*nfamilies()),nregions2()*nfamilies())))) %>%
#   left_join(randint2(), by = "familyid") %>%
#   mutate(familyid = as.factor(familyid))
# 
# finaldat = simple %>%
#   # make normal predictor
#   mutate(x1 = (rnorm(sampsize())),
#          # make binary demographic variable (male/female)
#          female = rep(c(0,1), times = sampsize()/2)) %>%
#   mutate(
#     # now make linear predictor model
#     z = randint_reg + randintfam + x1 + input$groupeffect*female + (rnorm(sampsize, sd = sqrt(input$residvar))), # add random effects, fixed effects, and random noise
#     # trying out 2 ways to make binary outcome
#     pr1 = exp(z)/(1+exp(z)),
#     runis = runif(sampsize, 0, 1),
#     y1 = ifelse(runis < pr1, 1, 0))
