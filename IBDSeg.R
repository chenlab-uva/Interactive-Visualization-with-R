# IBDseg with UK Biobank data
library(shiny)
library(dplyr)
library(ggplot2)

prefix <- "ukb"

#' file names
seg_name <- paste0(prefix, ".seg")
segments_name <- paste0(prefix, ".segments.gz")
all_seg_name <- paste0(prefix, "allsegs.txt")
if( !(file.exists(seg_name) & file.exists(segments_name) & file.exists(all_seg_name)) ) stop("Missing IBDSeg files")

#' load inference info
individuals_all <- read.table(seg_name, header = TRUE, stringsAsFactors = FALSE)


# segments considered
all_seg <- read.table(all_seg_name, header = TRUE)
all_seg <- subset(all_seg, select = c(Chr, StartMB, StopMB))
all_seg <- all_seg[all_seg[,"Chr"]!=23,]
colnames(all_seg) <- c("chr", "start", "end")

# king segments gz file
segments <- read.table(segments_name, header = TRUE)
segments <- subset(segments, select = c(ID1, ID2, IBDType, Chr, StartMB, StopMB))
colnames(segments) <- c("ID1", "ID2", "IBDType", "chr", "start", "end")



# provides the user interface for the app.
ui <- fluidPage(
  
  # App title ----
  titlePanel("Options"),
  sidebarLayout(position = "left",
                sidebarPanel(
                  #id = "sidebar", paste("Sidebar: Often reserved for Inputs"), 
                  selectInput("InfType", paste("Choose a InfType:"),
                              choices =c(Choose='',unique(individuals_all$InfType))),
                  sliderInput("IBD1_Seg_Range", "Length Proportion of IBD1 Segments:",min = 0, max = 1,value = c(0,1),step = 0.001),
                  sliderInput("IBD2_Seg_Range", "Length Proportion of IBD2 Segments:",min = 0, max = 1,value = c(0,1)),step = 0.001),
                # Sidebar layout with input and output definitions ----
                mainPanel(
                  tabsetPanel(
                    tabPanel("List",dataTableOutput(outputId = "dt1")),
                    tabPanel("Plot", 
                             plotOutput(outputId = "plot1",click = "plot_click", width = "70%"), verbatimTextOutput("click_info"),
                             plotOutput(outputId = "plot2",height = "680px", width = "70%")
                             ),
                    tabPanel("Table",
                             dataTableOutput(outputId = "dt2")
                             )
                    )
                )
  )
)


# server
server <- function(input, output) {
  output$dt1 <- renderDataTable({filter(individuals_all, InfType==input$InfType) %>% select("ID1","ID2","IBD1Seg","IBD2Seg","PropIBD","InfType")}, 
                                options=list(searching = FALSE))
  
  output$dt2 <- renderDataTable({
    individuals_all %>%
      group_by(InfType) %>% 
      summarize(n())},
    options=list(searching = FALSE))
  
  output$plot1 <- renderPlot({
    individuals_all_select <- filter(individuals_all,IBD1Seg >= input$IBD1_Seg_Range[1] & IBD1Seg <= input$IBD1_Seg_Range[2],
                                     IBD2Seg >= input$IBD2_Seg_Range[1] & IBD2Seg <= input$IBD2_Seg_Range[2])
    
    d0 <- individuals_all_select$IBD2Seg>0.7
    d1.PO <- (!d0) & individuals_all_select$IBD1Seg+individuals_all_select$IBD2Seg>0.96 | (individuals_all_select$IBD1Seg+individuals_all_select$IBD2Seg>0.9 & individuals_all_select$IBD2Seg<0.08)
    d1.FS <- (!d0) & (!d1.PO) & individuals_all_select$PropIBD>0.35355 & individuals_all_select$IBD2Seg>=0.08
    d2 <- individuals_all_select$PropIBD>0.17678 & individuals_all_select$IBD1Seg+individuals_all_select$IBD2Seg<=0.9 & (!d1.FS)
    d3 <- individuals_all_select$PropIBD>0.08839 & individuals_all_select$PropIBD<=0.17678
    d4 <- individuals_all_select$PropIBD>0.04419 & individuals_all_select$PropIBD<=0.08839
    dU <- individuals_all_select$PropIBD>0 & individuals_all_select$PropIBD<=0.04419
    plot(individuals_all_select$IBD1Seg[dU], individuals_all_select$IBD2Seg[dU], type="p", col = "black", cex.lab=1.2,
         xlim=c(min(individuals_all_select$IBD1Seg), max(individuals_all_select$IBD1Seg)),
         ylim=c(min(individuals_all_select$IBD2Seg), max(individuals_all_select$IBD2Seg)),
         main = paste("IBD Segments In Inferred",prefix,"Relatives"),
         xlab=expression(paste("Length Proportion of IBD1 Segments (", pi[1], ")",sep="")),
         ylab=expression(paste("Length Proportion of IBD2 Segments (", pi[2], ")",sep="")))
    points(individuals_all_select$IBD1Seg[d0], individuals_all_select$IBD2Seg[d0], col="purple")
    points(individuals_all_select$IBD1Seg[d1.PO], individuals_all_select$IBD2Seg[d1.PO], col="red")
    points(individuals_all_select$IBD1Seg[d1.FS], individuals_all_select$IBD2Seg[d1.FS], col="green")
    points(individuals_all_select$IBD1Seg[d2], individuals_all_select$IBD2Seg[d2], col="blue")
    points(individuals_all_select$IBD1Seg[d3], individuals_all_select$IBD2Seg[d3], col="magenta")
    points(individuals_all_select$IBD1Seg[d4], individuals_all_select$IBD2Seg[d4], col="gold")
    abline(h = 0.08, col = "green", lty = 3, lwd = 2)
    abline(a = 0.96, b = -1, col = "red", lty = 3, lwd = 2)
    abline(a = 0.3535534, b = -0.5, col = "green", lty = 3, lwd = 2)
    abline(a = 0.1767767, b = -0.5, col = "blue", lty = 3, lwd = 2)
    abline(a = 0.08838835, b = -0.5, col = "magenta", lty = 3, lwd = 2)
    abline(a = 0.04419, b = -0.5, col = "gold", lty = 3, lwd = 2)
    allcolors <- c("purple", "red", "green", "blue", "magenta", "gold", "black")
    legend("topright", c("Inferred MZ", "Inferred PO", "Inferred FS", "Inferred 2nd", "Inferred 3rd", "Inferred 4th", "Inferred UN"),
           col=allcolors, text.col = allcolors, pch = 19, cex = 1.2)
  })
  output$click_info <- renderPrint({
    if(!is.null(input$plot_click)){
      point.index <- which.min((individuals_all$IBD1Seg-input$plot_click$x)^2+(individuals_all$IBD2Seg-input$plot_click$y)^2) 
      if (!(abs(individuals_all[point.index, "IBD1Seg"]-input$plot_click$x) <= 0.01 & abs(individuals_all[point.index,"IBD2Seg"]-input$plot_click$y) <= 0.01)) {
        target.data <- NULL}
      else{
        segments$IBDType <- factor(segments$IBDType, levels = c("IBD0", "IBD1", "IBD2"))
        target.data <- filter(segments, ID1==individuals_all[point.index,"ID1"], ID2==individuals_all[point.index,"ID2"])
      }
      output$plot2 <- renderPlot({
        validate(
          need(nrow(target.data) > 0, "Please select a related pair")
        )
        theme_set(theme_bw(base_size = 18))
        g <- ggplot() +
          geom_rect(data = all_seg, aes(xmin = start, xmax = end, ymin = 0, max = 0.9), fill = 'white', color = "black", size = 0.85) + 
          geom_rect(data = target.data , aes(xmin = start, xmax = end, ymin = 0, ymax = 0.9, fill = IBDType)) + 
          geom_rect(data = all_seg, aes(xmin = start, xmax = end, ymin = 0, max = 0.9), color = "black", alpha = 0, size = 0.85) + 
          scale_fill_manual(values = c("IBD0" = "white", "IBD1" = "dodgerblue2", "IBD2" = "firebrick2"), drop = FALSE) + 
          facet_grid(chr ~ .) +
          labs(x = "Position (Mb)", y = "", title=paste(prefix, "IBDSeg", target.data$ID1, target.data$ID2, sep = "_"))+
          theme(
            legend.position = "bottom", legend.key = element_rect(color = "black"),
            panel.background = element_rect(fill = 'grey80', color = 'grey80'), panel.border = element_blank(),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            axis.text.y = element_blank(), axis.ticks.y = element_blank()
          )
        print(g)
      })
    }
  })
}

# Run the application 
shinyApp(ui = ui, server = server)

