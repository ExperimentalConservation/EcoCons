
####################################
#GLV model for ecology course
# https://stefanoallesina.github.io/Sao_Paulo_School/intro.html
####################################

library(shiny)
library(deSolve)
library(shinyMatrix)
library(tidyverse)

################################################################################
## define some functions
################################################################################

# Generalized Lotka-Volterra model
GLV <- function(t, x, parameters){
  with(as.list(c(x, parameters)), {
    x[x < 10^-8] <- 0 # prevent numerical problems
    dxdt <- x * (r + A %*% x)
    list(dxdt)
  })
}
# function to plot output
plot_ODE_output <- function(out){
  #out <- as.data.frame(out)
  colnames(out) <- c("time", c("Belpharisma", "Paramecium", "Spirostomum", "Didinium", "Stenostomum"))
  out <- as_tibble(out) %>% gather(species, density, -time)
  return(out)
}

# general function to integrate GLV
integrate_GLV <- function(r, A, x0, maxtime = 100, steptime = 0.5){
  times <- seq(0, maxtime, by = steptime)
  parameters <- list(r = r, A = A)
  # solve numerically
  out <- ode(y = x0, times = times, 
             func = GLV, parms = parameters, 
             method = "ode45")
  # plot and make into tidy form
  out <- plot_ODE_output(out)
  return(out)
}

################################################################################
# ##initial values
start_matrix <- -matrix(runif(25), 5, 5, 
       dimnames = list(c("Belph", "Para", "Spiro", "Didi", "Steno"), c("Belph", "Para", "Spiro", "Didi", "Steno")))

 start_r_matrix <- matrix(runif(5), dimnames = list(c("Belpharisma", "Paramecium", "Spirostomum", "Didinium", "Stenostomum")))
 
 start_k_matrix <- matrix(round(runif(5, 10, 500)), dimnames = list(c("Belpharisma", "Paramecium", "Spirostomum", "Didinium", "Stenostomum")))

################################################################################
# Define UI for application
ui <- fluidPage(
  
  # Application title
  titlePanel("GLV model"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      style = "height: 90vh; overflow-y: auto;",
      br(),
      matrixInput("rmat",
                  "Growth rates of species:",
                  value = start_r_matrix,
                  rows = list(extend = F, delete = F, names=T),
                  cols = list(extend = F, delete = F, names=T),
                  class = "numeric",
                  lazy=F),
      matrixInput("kmat",
                  "Carrying capacity for each species (per ml):",
                  value = start_k_matrix,
                  rows = list(extend = F, delete = F, names=T),
                  cols = list(extend = F, delete = F, names=T),
                  class = "numeric",
                  lazy=F),
      matrixInput("matrix1",
                  "Interaction matrix:",
                  value = start_matrix,
                  rows = list(extend = F, delete = F, names=T),
                  cols = list(extend = F, delete = F, names=T),
                  class = "numeric",
                  lazy=F
      ),
      sliderInput("time",
                  "Length of simulation (days):",
                  min = 1,
                  max = 50,
                  value = 7),
      sliderInput("as1",
                  "Starting abundance (Belpharisma):",
                  min = 50,
                  max = 500,
                  value = 50),
      sliderInput("as2",
                  "Starting abundance (Paramecium):",
                  min = 50,
                  max = 500,
                  value = 50),
      sliderInput("as3",
                  "Starting abundance (Spirostomum):",
                  min = 50,
                  max = 200,
                  value = 50),
      sliderInput("as4",
                  "Starting abundance (Didinium):",
                  min = 5,
                  max = 50,
                  value = 5),
      sliderInput("as5",
                  "Starting abundance (Stenostomum):",
                  min = 5,
                  max = 20,
                  value = 5),
      #,
      uiOutput("Empty")
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("dynplot"),
      tableOutput("table")
    )
  )
)

####
# Define server logic required to run the functions
server <- function(input, output, session) {
  
  data <- reactiveValues(plt = NULL)
  
    output$dynplot <- renderPlot({
    ##the matrixxw
     A_4 <- input$matrix1
     diag(A_4)[which(input$rmat>0)] <- -input$rmat[which(input$rmat>0)]/input$kmat[which(input$rmat>0)]
     diag(A_4)[which(input$rmat<0)] <- input$rmat[which(input$rmat<0)]/input$kmat[which(input$rmat<0)]
    # #A_4 <- start_matrix
    # A_4 <- -matrix(runif(25)/100, 5, 5, 
    #                       dimnames = list(paste("Species",1:5,sep=" "), 
    #                                       paste("Species",1:5,sep=" ")))
    # # ##make sure values are negative
    # # #A_4 <- A_4
    # # ##generate some starting growth values
     r_4 <- input$rmat
    # # ##and some initial abundances
     x0_4 <- c(input$as1, input$as2, input$as3, input$as4, input$as5)/75

    # # ##run the GLV
     plt<-integrate_GLV(r_4, A_4, x0_4, maxtime = input$time)
     ##reshape for displaying
     disp_plt<-plt %>% filter(time==max(plt$time)) %>% 
        select(c(species, density)) %>%
        pivot_wider(names_from=species, values_from = density)
     rownames(disp_plt)<-"Abundance"
     data$plt <- disp_plt
    # # 
    # # ##plot the output
    ggplot(data = plt) +
      aes(x = time, y = density*75, colour = species) +
      geom_line() +
      theme_classic()+
      xlab("Day") +
      ylab("Landscape abundane")
    
    })
    

    output$table=renderTable({
      data$plt
    })
    
  #})
  
}

# Run the application 
shinyApp(ui = ui, server = server)
