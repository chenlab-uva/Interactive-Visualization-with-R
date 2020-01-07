library(ggplot2)
library(shiny)
library(dplyr)

prefix="ukb"
seg_name <- paste0(prefix, ".roh")
segments_name <- paste0(prefix, ".rohseg.gz")
all_seg_name <- paste0(prefix, "allsegs.txt")
if( !(file.exists(seg_name) & file.exists(segments_name) & file.exists(all_seg_name)) ) stop("Missing RoH files")
all_seg <- read.table(all_seg_name, header = TRUE)
all_seg <- subset(all_seg, select = c(Chr, StartMB, StopMB))
segments <- read.table(segments_name, header = TRUE)
segments <- subset(segments, select = c(FID, ID, Chr, StartMB, StopMB))
roh <- read.table(seg_name, header = TRUE)
roh_info <- roh[roh$F_ROH > 2^-4.5, c("FID","ID","F_ROH_X","F_ROH")]
roh_info$FID <- as.character(roh_info$FID)
roh_info$ID <- as.character(roh_info$ID)

ui <- fluidPage(
  titlePanel("Filter"),
  sidebarLayout(position = "left",
                sidebarPanel(id = "sidebar",
                             sliderInput("F_ROH_range", "F_ROH_Range:", min = min(roh_info$F_ROH), max = max(roh_info$F_ROH),
                                         value = c(min(roh_info$F_ROH),max(roh_info$F_ROH))),
                             sliderInput("F_ROH_X_range", "F_ROH_X_Range:", min = min(roh_info$F_ROH_X), max = max(roh_info$F_ROH_X),
                                         value = c(min(roh_info$F_ROH_X),max(roh_info$F_ROH_X)))),
                mainPanel(id="main",
                          tabPanel("Plots", 
                                   fluidRow(
                                     splitLayout(style = "border: 1px solid silver:", 
                                                 cellWidths = c("50%", "50%"),
                                                 #PLOT1
                                                 plotOutput(outputId = "plot1", click = "plot_click",height = "600px"),
                                                 #PLOT2
                                                 plotOutput(outputId = "plot2", height = "600px", width = "100%")
                                                 )),
                                   fluidRow(
                                     column(width = 5, verbatimTextOutput("click_info"),verbatimTextOutput("last_infor"))
                                     )
                                   )
                          )
                )
)

server <- function(input, output) {
  output$plot1 <- renderPlot({
    target.data <- filter(roh_info, F_ROH >= input$F_ROH_range[1] & F_ROH <= input$F_ROH_range[2],
                          F_ROH_X >= input$F_ROH_X_range[1] & F_ROH_X <= input$F_ROH_X_range[2])
    plot(target.data$F_ROH_X, target.data$F_ROH,xlab = "F_ROH_X", ylab="F_ROH", main = paste0("F_ROH vs F_ROH_x in ", prefix))
  })
  
  output$click_info <- renderPrint({
    if(!is.null(input$plot_click)){
      #print(c(input$plot_click$x, input$plot_click$y))
      min.index <- which.min(abs(roh_info$F_ROH-input$plot_click$y)^2 + abs(roh_info$F_ROH_X-input$plot_click$x)^2)
      name <- roh_info[min.index,"ID"]
      #print(name)
      if (!(abs(roh_info[min.index,"F_ROH_X"]-input$plot_click$x) <= 0.01 & abs(roh_info[min.index,"F_ROH"]-input$plot_click$y) <= 0.01)) {
        k <- NULL
      } else {
        k <- subset(segments,ID==name)
        print(k) 
      }
      output$plot2 <- renderPlot({
        validate(
          need(nrow(k) > 0, "Please select a related pair")
        )
        theme_set(theme_bw(base_size = 18))
        f_roh <- roh_info[roh_info$ID==name,"F_ROH"]
        fid <- as.character(k[1,1])
        id <- as.character(k[1,2])
        g <- ggplot() +
          geom_rect(data = all_seg, aes(xmin = StartMB, xmax = StopMB, ymin = 0, max = 0.9), fill = 'white', color = "black", size = 0.85) +
          geom_rect(data = k, aes(xmin = StartMB, xmax = StopMB, ymin = 0, ymax = 0.9), fill = "red") +
          geom_rect(data = all_seg, aes(xmin = StartMB, xmax = StopMB, ymin = 0, max = 0.9), color = "black", alpha = 0, size = 0.85) +
          facet_grid(Chr ~ .) + scale_x_continuous(expand  = c(0, 0), limits = c(0, NA)) +
          labs(x = "Position (Mb)", y = "", title = bquote(paste('Run of Homozygosity for ', .(id), ' from FAM ', .(fid), ' in ', .(prefix), ' (F'['ROH']*' = ', .(f_roh), ')'))) +
          theme(legend.position = "none",
                panel.background = element_rect(fill = 'grey80', color = 'grey80'), panel.border = element_blank(),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.title=element_text(size = 12))
        print(g)
        
      })
    }
  })
}

shinyApp(ui, server)
