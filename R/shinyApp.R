library(shiny)
library(dplyr)
library(Seurat)
library(RColorBrewer)
source('./visualization.R')

getVarib = function(sc){
  tryCatch({
    sc = Seurat::FindVariableFeatures(sc,selection.method = 'vst',verbose=F)
    return('vst')
  },error=function(e){
    tryCatch({
      sc = Seurat::FindVariableFeatures(sc,selection.method= 'disp',verbose=F)
      return('disp')
    },error=function(e){
      sc = Seurat::FindVariableFeatures(sc,selection.method= 'mvp',verbose=F)
      return('mvp')
    })
  })
}

score = readRDS('../data/test.rds')
sc = Seurat::CreateSeuratObject(score)
sc = Seurat::ScaleData(sc)
sc = Seurat::FindVariableFeatures(sc,selection.method = getVarib(sc),verbose=F)
sc = Seurat::RunPCA(sc,verbose=F,npcs = 10)
sc = Seurat::RunUMAP(sc,dims = 1:10, n.components = 2,
                     verbose = F)
sc = Seurat::RunTSNE(sc,dims = 1:10, verbose = F)
sc = Seurat::FindNeighbors(sc,dims = 1:10)
sc = Seurat::FindClusters(sc,resolution=0.3)

ui = fluidPage(

  # App title
  titlePanel("PAS results browser"),

  # Sidebar layout with a input and output definitions
  sidebarLayout(
    sidebarPanel(

      selectInput(inputId = 'dimens', label = 'dimensional reduction method:',
                  choices = c('umap','tsne'),
                  selected = 'umap'),

      sliderInput(inputId = 'resol', label = 'resolution for clustering:',
                  min = 0, max = 1,
                  value = 0.3),
      br(),
      br(),

      # Show marker pathway
      checkboxInput(inputId = "show_marker",
                    label = "Show marker pathwats",
                    value = F),

      sliderInput(inputId = 'n_marker', label = 'number of marker pathways:',
                  min = 1, max = 20,
                  value = 5),

      sliderInput(inputId = 'logFC', label = 'threshold of logFC:',
                  min = 0, max = 1,
                  value = 0.25),

      sliderInput(inputId = 'pct', label = 'threshold of pct:',
                  min = 0, max = 0.5,
                  value = 0.1),
    ),

    mainPanel(
      plotOutput(outputId = 'dimensionalReduction'),
      br(),    # a little bit of visual separation

      # Show data table
      plotOutput(outputId = 'heatmap')
    )
  )
)




server = function(input, output){
  scPCA = reactive({
    req(input$resol)
    sc = Seurat::FindNeighbors(sc,dims = 1:10)
    sc = Seurat::FindClusters(sc,resolution=input$resol)
  })

  markers = reactive({
    req(input$logFC, input$pct)
    all.markers = Seurat::FindAllMarkers(scPCA(),
                                logfc.threshold = input$logFC,
                                min.pct = input$pct)
  })

  output$dimensionalReduction = renderPlot(
    Seurat::DimPlot(scPCA(),pt.size = 2,reduction = input$dimens)
  )


  output$heatmap <- renderPlot(
    if(input$show_marker){
      all.markers = markers()
      top = all.markers %>% group_by(cluster) %>% top_n(n = input$n_marker,
                                                        wt = avg_logFC)
      Seurat::DoHeatmap(scPCA(),
                        features = top$gene) +
        NoLegend()+
        scale_fill_gradientn(colors =
                               colorRampPalette(c("forestgreen", "gray90", "orange"))(300))
    }
  )
}


shinyApp(ui, server)
