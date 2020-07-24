
library(dplyr)


sc = .GlobalEnv$.sc_oj
all.pathways = .GlobalEnv$.pathways

GM_state <- function(cds,start_node){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$cluster)[,start_node]
    return(as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))]))
  }else {
    return (1)
  }
}




ui = fluidPage(

  # App title
  titlePanel("PAS results browser"),

  ## clustering
  fluidRow(
    sidebarLayout(
      sidebarPanel(
        selectInput(inputId = 'dimens', label = 'dimensional reduction method:',
                    choices = c('umap','tsne'),
                    selected = 'umap'),
      ),
      mainPanel(
        h3("clustering analysis"),
      )
    )
  ),


  fluidRow(
    sidebarLayout(
      sidebarPanel(
        sliderInput(inputId = 'resol', label = 'resolution for clustering:',
                    min = 0, max = 1,
                    value = 0.3)
      ),
      mainPanel(
        plotOutput(outputId = 'dimensionalReduction')
      )
    )
  ),



  fluidRow(
    sidebarLayout(
      sidebarPanel(
        selectInput(inputId = 'feature', label = 'pathways:',
                    choices = all.pathways)
      ),
      mainPanel(
        plotOutput(outputId = 'featureplot')
      )
    )
  ),

  ## differential analysis

  fluidRow(
    sidebarLayout(
      sidebarPanel(
        h3("Differential anlysis"),
        checkboxInput(inputId = "show_marker",
                      label = "Go!",
                      value = F),
      ),
      mainPanel(
        h3("signature pathways")

      )
    )
  ),

  ## differential analysis
  fluidRow(
    sidebarLayout(
      sidebarPanel(

        sliderInput(inputId = 'n_marker', label = 'number of signatures for each cluster:',
                    min = 1, max = 20,
                    value = 5),
        sliderInput(inputId = 'logFC', label = 'threshold of logFC:',
                    min = 0, max = 1,
                    value = 0.1),

        sliderInput(inputId = 'pct', label = 'threshold of pct:',
                    min = 0, max = 0.5,
                    value = 0.1)
      ),
      mainPanel(
        plotOutput(outputId = 'heatmap'),
        DT::dataTableOutput(outputId = "featuretable")

      )
    )
  ),

  ## Trajectory analysis
  fluidRow(
    sidebarLayout(
      sidebarPanel(
        h3("Trajectory anlysis"),
        checkboxInput(inputId = "show_traject",
                      label = "Go!",
                      value = F),
      ),
      mainPanel(
        h3("pseudo time")
      )
    )
  ),

  fluidRow(
    sidebarLayout(
      sidebarPanel(
        textInput(inputId = 'start_cluster', label = 'start cluster','0'),
      ),
      mainPanel(
        fluidRow(
          splitLayout(cellWidths = c("50%", "50%"),
                      plotOutput("tarjectory"), plotOutput("pseudotime"))
        )
      )
    )
  )

)




server = function(input, output){

  ## update seurat object if resolution changed
  scPCA = reactive({
    req(input$resol)
    sc = Seurat::FindNeighbors(sc)
    sc = Seurat::FindClusters(sc,resolution=input$resol)
  })




  ## dimensional reduction plot   UMAP/TSNE

  output$dimensionalReduction = renderPlot(
    Seurat::DimPlot(scPCA(),pt.size = 2,
                    reduction = input$dimens)
  )

  ## feature plot

  output$featureplot = renderPlot(
    Seurat::FeaturePlot(scPCA(),features = input$feature,
                        reduction = input$dimens,cols = c("lightgrey", "red"))
  )

  ## finding markers based on logFC and pct
  markers = reactive({
    req(input$logFC, input$pct)
    all.markers = Seurat::FindAllMarkers(scPCA(),
                                         logfc.threshold = input$logFC,
                                         min.pct = input$pct)
  })


  ## update monocle object based on clusters from seurat
  monocle_cds = reactive({
    req(input$resol)
    ss = scPCA()
    pdata  = data.frame(Idents(ss),
                        row.names = colnames(GetAssayData(ss)))
    colnames(pdata) = c('cluster')
    fdata <- data.frame(gene_short_name = row.names(GetAssayData(ss)),
                        row.names = row.names(GetAssayData(ss)))
    cds <- monocle::newCellDataSet(GetAssayData(ss),
                          phenoData = new('AnnotatedDataFrame', data = pdata),
                          featureData = new('AnnotatedDataFrame', data = fdata),
                          expressionFamily = monocle::tobit())
    cds <- monocle::reduceDimension(cds, max_components = 2, method = 'DDRTree')
    cds <- monocle::orderCells(cds)
  })


  ## update pseudo time based on start node
  monocle_time = reactive({
    req(input$start_cluster)
    cds <- monocle_cds()
    cds <- monocle::orderCells(cds, root_state = GM_state(cds,input$start_cluster))
  })




  ## heatmap of signature pathways
  output$heatmap <- renderPlot(
    if(input$show_marker){
      all.markers = markers()
      top = all.markers %>% group_by(cluster) %>%
        top_n(n = input$n_marker,wt = avg_logFC) %>%
        as.data.frame
      Seurat::DoHeatmap(scPCA(),
                        features = top$gene) +
        NoLegend()+
        scale_fill_gradientn(colors =
                               colorRampPalette(c("forestgreen", "gray90", "orange"))(300))

    }
  )


  ## display signature table
  output$featuretable <- DT::renderDataTable(
    if(input$show_marker){
      all.markers = markers()
      top = all.markers %>% group_by(cluster) %>%
        top_n(n = input$n_marker,wt = avg_logFC) %>%
        as.data.frame
      colnames(top)[ncol(top)] = 'pathway'
      DT::datatable(data = top,
                    options = list(pageLength = 5),
                    rownames = F)
    }
  )



  ## trajectory from monocle
  output$tarjectory <- renderPlot(
    if(input$show_traject){
      monocle::plot_cell_trajectory(monocle_cds(), color_by = "cluster")
    }
  )

  ## pseudo time
  output$pseudotime <- renderPlot(
    if(input$show_traject){
      monocle::plot_cell_trajectory(monocle_time(), color_by = "Pseudotime")
    }
  )

}


shinyApp(ui, server)
