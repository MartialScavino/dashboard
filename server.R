server <- function(input, output, session) {
  
  # Reactive value for dataset
  val <- reactiveValues()
  
  
  # Loading dataset
  observeEvent(input$data,{
    # Check for file type
    if (input$data$type != "" & !(str_sub(input$data$name,-3,-1) == "rds")){ 
      alert("Enter a rds file")
      return(0)
    }
    
    df <- readRDS(input$data$datapath)

    if (!("percent.mt" %in% names(df@meta.data))){
      df[["percent.mt"]] <- PercentageFeatureSet(df, pattern = "^MT-")
    }  
    
    val$data <- df
    
    # Show Metadata
    output$dataset <- renderDataTable(val$data@meta.data, extensions = 'Buttons', 
                                      options = list(dom = 'Bfrtip', fixedColumns = TRUE,
                                                     buttons = c('copy', 'csv', 'excel')))
    })
 
  observeEvent(input$filtered, {

    addClass(id = "UpdateAnimateLoad", class = "loading dots")
    disable("filtered")

    counts <- ReadMtx(input$matrixfile$datapath,
            cells = input$barcodesfile$datapath,
            features = input$featuresfile$datapath,
            feature.column = 1)

    df <- CreateSeuratObject(counts = counts)

    if (!("percent.mt" %in% names(df@meta.data)))
      df[["percent.mt"]] <- PercentageFeatureSet(df, pattern = "^MT-")

    val$data <- df

    # Show Metadata
    output$dataset <- renderDataTable(val$data@meta.data, extensions = 'Buttons',
                                      options = list(dom = 'Bfrtip', fixedColumns = TRUE,
                                                     buttons = c('copy', 'csv', 'excel')))

    enable("filtered")
    removeClass(id = "UpdateAnimateLoad", class = "loading dots")
  })
    
  
  # Slider featre slide bar
  output$sliderFeature <- renderUI(sliderInput("features",
                                               "Choissisez le nombre de gène 
                                               exprimé voulu dans chaque cellule",
                                               min = 0,
                                               max = max(val$data$nFeature_RNA),
                                               value = c(600, max(val$data$nFeature_RNA)), 
                                               step = 50))
  
  # Number of cell removed text
  output$texte <- renderText({
    
    nb_cell_bad <- length(rownames(val$data@meta.data[which(val$data$percent.mt > input$mt | 
                                                              val$data$nFeature_RNA < input$features[1] | 
                                                              val$data$nFeature_RNA > input$features[2]),]))
    nb_cell_tot <- length(rownames(val$data@meta.data))
    
    paste("On enlève", nb_cell_bad, "cellules du jeu de donnée soit", round(100*nb_cell_bad/nb_cell_tot, 2), "% des cellules")
    
  })
  
  output$scatter_QC_MT <- renderPlot({
    
    
    ggplot(val$data@meta.data, aes(nCount_RNA, percent.mt, color = nFeature_RNA)) +
      geom_point(data = val$data@meta.data[val$data$percent.mt > input$mt,], alpha = 0.2) + 
      geom_point(data = val$data@meta.data[val$data$percent.mt <= input$mt,],) + 
      geom_hline(yintercept = input$mt, col = "#FEB078") + 
      scale_color_viridis(option = "B") + theme_light()
    
    
  })
    
  output$scatter_QC_Feature <- renderPlot({
    
    ggplot(val$data@meta.data, aes(nCount_RNA, nFeature_RNA, color = percent.mt)) +
      geom_point(data = val$data@meta.data[val$data$nFeature_RNA < input$features[1],], alpha = 0.2) + 
      geom_point(data = val$data@meta.data[val$data$nFeature_RNA > input$features[2],], alpha = 0.2) + 
      geom_point(data = val$data@meta.data[val$data$nFeature_RNA <= input$features[2] & val$data$nFeature_RNA >= input$features[1],]) + 
      scale_color_viridis() +
      geom_hline(yintercept = input$features, col = "#832681") + theme_light()
    
  })
  
  
  output$hist_QC_MT <- renderPlot({
    
    ggplot(val$data@meta.data, aes(x = percent.mt)) + 
      geom_histogram(aes(y = after_stat(density)), bins = 200, fill = "#FEB078", ) + 
      geom_density(color = "#f8765c") + 
      geom_vline(xintercept = input$mt, col = "#800000") + theme_light()
  
    })
  
  output$hist_QC_Feature <- renderPlot({
    
    ggplot(val$data@meta.data, aes(x = nFeature_RNA)) + 
      geom_histogram(aes(y = after_stat(density)), bins = 200, fill = "#DDA0DD", ) + 
      geom_density(color = "#8B0000") + 
      geom_vline(xintercept = input$features, col = "#832681") + theme_light()
    
  })
  
  observeEvent(input$Trim, {
    
    addClass(id = "UpdateAnimateTrim", class = "loading dots")
    disable("Trim")
    
    val$data <- subset(val$data,
                       subset = nFeature_RNA > input$features[1] & nFeature_RNA < input$features[2] & percent.mt < input$mt)
    
    enable("Trim")
    removeClass(id = "UpdateAnimateTrim", class = "loading dots")
    
  })
  
  
  observeEvent(input$LaunchPreprocessing, {
    
    addClass(id = "UpdateAnimatePreprocessing", class = "loading dots")
    disable("LaunchPreprocessing")
    
    val$data <- NormalizeData(val$data, normalization.method = input$NormMethod, 
                              scale.factor = input$ScaleFactor)
    val$data <- FindVariableFeatures(val$data, selection.method = input$VariableMethod,
                                     nfeatures = input$nfeatures)
    
    s.genes <- cc.genes$s.genes
    g2m.genes <- cc.genes$g2m.genes
    
    # Convert S genes
    s.genes_mus.musculus <- orthologs(genes = s.genes, species = "mouse")
    
    # Convert G2M genes
    g2m.genes_mus.musculus <- orthologs(genes = g2m.genes, species = "mouse")
    
    if (input$species == "Human")
      val$data <- CellCycleScoring(val$data, s.features = s.genes, g2m.features = g2m.genes)
    
    else
      val$data <- CellCycleScoring(val$data, s.features = s.genes_mus.musculus$symbol, g2m.features = g2m.genes_mus.musculus$symbol)
    
    
    dict <- list("All genes" = rownames(val$data), "Variable genes" = VariableFeatures(val$data))
    
    if (input$regressing)
      val$data <- ScaleData(val$data, features = as.vector(unlist(dict[input$FeatureScale])), 
                            vars.to.regress = c("S.Score", "G2M.Score"))

    else
      val$data <- ScaleData(val$data, features = as.vector(unlist(dict[input$FeatureScale])))
    
    
    enable("LaunchPreprocessing")
    removeClass(id = "UpdateAnimatePreprocessing", class = "loading dots")
    
  })
  
  
  output$variablefeatureplot <- renderPlot({
    
    if (length(VariableFeatures(val$data)) == 0)
      return(0)
    
    top20 <- head(VariableFeatures(val$data), 20)
    p1 <- VariableFeaturePlot(val$data)
    LabelPoints(plot = p1, points = top20, repel = TRUE, xnudge = 0, ynudge = 0)
    
    })
  
  
  output$groupbypca <- renderUI(selectInput("groupbypca", "group cells by",
                                         choices = names(val$data@meta.data)))
  
  observeEvent(input$dopca, {
    
    addClass(id = "UpdateAnimatePCA", class = "loading dots")
    disable("dopca")
    
    dict <- list("All genes" = rownames(val$data), "Variable genes" = VariableFeatures(val$data))
    
    val$data <- RunPCA(val$data, features = as.vector(unlist(dict[input$featurePCA])))
    
    enable("dopca")
    removeClass(id = "UpdateAnimatePCA", class = "loading dots")
  })
  
  
  
  output$PCA <- renderPlot({
    if (!("pca" %in% names(val$data@reductions)))
      return(0)
    
    if(length(table(val$data@meta.data[,input$groupbypca])) > 50)
      FeaturePlot(val$data, reduction = "pca", features = input$groupbypca)
    
    else
      DimPlot(val$data, reduction = "pca", group.by = input$groupbypca)
    })
  
  output$elbow <- renderPlot({
    if (!("pca" %in% names(val$data@reductions)))
      return(0)
    ElbowPlot(val$data, ndims = input$ndimElbow)
    })
  
  output$dimheatmap <- renderPlot({
    if (!("pca" %in% names(val$data@reductions)))
      return(0)
    DimHeatmap(val$data, dims = 1:input$ndimheatmap, balanced = TRUE, nfeatures = 20)
    })
  
  
  
  output$groupbyumap <- renderUI(selectInput("groupbyumap", "group cells by",
                                             choices = names(val$data@meta.data)))
  
  observeEvent(input$doumap,{
    
    addClass(id = "UpdateAnimateUMAP", class = "loading dots")
    disable("doumap")
      
               val$data <- RunUMAP(val$data, dims = 1:input$ndimumap)
               
    enable("doumap")
    removeClass(id = "UpdateAnimateUMAP", class = "loading dots")
               })
  
  output$UMAP <- renderPlot({
    if(!("umap" %in% names(val$data@reductions)))
      return(0)
    
    if(length(table(val$data@meta.data[,input$groupbyumap])) > 50)
      FeaturePlot(val$data, reduction = "umap", features = input$groupbyumap)
    
    else
      DimPlot(val$data, reduction = "umap", group.by = input$groupbyumap)
    
  })
  
  output$groupbycluster <- renderUI(selectInput("groupbycluster", "group cells by",
                                             choices = names(val$data@meta.data),
                                             selected = "seurat_clusters"))
  
  output$sliderndimclusters <- renderUI(
    sliderInput("ndimscluster", "Number of dimension", 
                min = 1, 
                max = length(val$data@reductions$pca), 
                value = 30, 
                step = 1)
  )
  
  observeEvent(input$docluster, {
    addClass(id = "UpdateAnimateCluster", class = "loading dots")
    disable("docluster")
    
    val$data <- FindNeighbors(val$data, k.param = input$kparam, dims = 1:input$ndimscluster)
    val$data <- FindClusters(val$data, resolution = input$resolution, algorithm = as.integer(input$algocluster))
    
    enable("docluster")
    removeClass(id = "UpdateAnimateCluster", class = "loading dots")
  })
    
  output$UMAPCluster <- renderPlot({
    if(length(table(val$data@meta.data[,input$groupbycluster])) > 50)
      FeaturePlot(val$data, reduction = "umap", features = input$groupbycluster)
    
    else
      DimPlot(val$data, reduction = "umap", group.by = input$groupbycluster)
  })
  
  output$savedata <- downloadHandler(
    filename = "seurat_object.rds",
    content = function(file) {
      addClass(id = "UpdateAnimateSave", class = "loading dots")
      disable("savedata")
      
      saveRDS(val$data, file)
      
      enable("savedata")
      removeClass(id = "UpdateAnimateSave", class = "loading dots")
    }
  )
  
  
  ## Feature Plot
  output$typeaheadFeature <- renderUI({
    
    genes <- data.frame(list(gene = rownames(val$data)))
    
    textInput.typeahead("gene", placeholder = "Enter a gene",
                        local = genes, 
                        valueKey = "gene",
                        tokens =  c(1:length(genes$gene)),
                        template = HTML("<p class='repo-language'>{{info}}</p> <p class='repo-gene'>{{gene}}</p>")
    )
  
    })
  
  liste = reactiveValues(df = data.frame(list(gene_name = character(0), show_button = character(0))))
  
  # Ajoute un gène à la liste lorsqu'on écrit un nouveau gène
  observeEvent(input$gene,{
    
    if (input$gene != "" & input$gene %in% liste$df$gene_name == FALSE){
      liste$df[nrow(liste$df) + 1,] <- c(input$gene, TRUE)
    }
    
    else if (input$gene != "" & input$gene %in% liste$df$gene_name == TRUE){
      if (liste$df[which(liste$df$gene_name == input$gene), "show_button"] == FALSE){
        
        liste$df[which(liste$df$gene_name == input$gene), "show_button"] <- TRUE
        
      }
    }
  })
  
  
  # Créer un bouton pour chaque élément dans la liste de gène signature
  observeEvent(liste$df$gene_name, output$list_button <- renderUI({
    
    if (length(liste$df$gene_name) > 0){
      tagList(
        lapply(1:length(liste$df$gene_name), function(i){
          
          if (liste$df$show_button[i] == TRUE){
            
            id1 <- paste0('slider_',liste$df$gene_name[i])
            actionButton(id1, liste$df$gene_name[i], "info", 
                         onclick = "Shiny.onInputChange('myclick', {id : this.id, val : this})")
          }
        })
      )
    }
    
    else{
      
      output$list_button <- renderUI(actionButton('bla', "No gene selected", "inverse"))
      
    }
  }))
  
  
  # Reset list quand on appuie sur le bouton
  observeEvent(input$reset, {
    
    liste$df <- data.frame(list(gene_name = character(0), show_button = character(0)))
    
  })
  
  # Make button removable when clicking on it
  observeEvent(input$myclick, {
    
    id <- strsplit(input$myclick$id, "_")[[1]][2]
    liste$df[which(liste$df$gene_name == id), "show_button"] <- FALSE
    
  })
  
  
  # Calcule le plot quand le bouton est sélectionné 
  observeEvent(input$dofeature,{
    
    output$featuresplot <- renderPlot({
      
      isolate({
        
        liste_gene_plot <- liste$df[which(liste$df$show_button == TRUE), "gene_name"]
        
        if (length(liste_gene_plot) == 1){
          
          p <- FeaturePlot(val$data, features = liste_gene_plot)
          return(p)
        }
        
        else if (length(liste_gene_plot) > 1){
          
          addClass(id = "UpdateAnimateFeature", class = "loading dots")
          disable("dofeature")
          
          val$data <- AddModuleScore_UCell(obj = val$data, features = list(Signature = liste_gene_plot))
          p <- FeaturePlot(object = val$data, features = "Signature_UCell") + ggtitle("Signature")
          
          enable("dofeature")
          removeClass(id = "UpdateAnimateFeature", class = "loading dots")
          return(p)
          
        }
        
        else { return() }
        
      })
    })
  })
  
  
}
