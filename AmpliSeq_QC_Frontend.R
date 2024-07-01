library(shiny)
library(shinyFiles)
library(bslib)
library(tidyverse)
library(ggridges)
library(ggstance)
library(ape)
library(phangorn)
library(ggtree)
#library(ggimage)
library(ggnewscale)
library(cowplot)

x.theme.axis.rotate.angle <- theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
legend.size <- theme(legend.key.size = unit(0.55,"line"))
standard.textsize <- 11
text.size.within <- (5/14)*(standard.textsize-2)
theme.text.size <- theme(text = element_text(size = standard.textsize))


# Define UI
ui <- fluidPage(
  titlePanel("Treponema AmpliSeq QC"),
  
  sidebarLayout(
    sidebarPanel(
      style = "position:fixed;width:22%;",
      shinyDirButton("directory", "Choose a Directory", "Please select a directory"),
      hr(),
      uiOutput("checkboxes"),
      hr(),
    ),
    
    mainPanel(
      tags$h3("Selected Directory"),
      verbatimTextOutput("selectedDirectory"), 
      tags$h3("Selected Files"),
      textOutput("selectedFiles"),
      hr(),
      tags$h3("Pipeline Status"),
      plotOutput(outputId = "p.pipeline.status"),
      hr(),
      tags$h3("Sequencing QC"), 
      
      conditionalPanel(
        condition = "output.folderSelected == true",
        p("For more detailed MultiQC Report, follow link"),
        uiOutput("multiqcLink")
      ),
      
      hr(),
      plotOutput(outputId = "p.readlengths"),
      hr(),
      plotOutput(outputId = "p.ontargetmapping"),
      hr(),
      tags$h4("Amplicon Coverage"), 
      navset_card_underline(
        nav_panel("Median", plotOutput("p.median.cov")),
        nav_panel("Minimum", plotOutput("p.min.cov")),
        nav_panel("%<10x", plotOutput("p.10x.cov"))
      ),
      hr(),
      tags$h3("Sample Relatedness"), 
      textOutput("SNPcount"),
      p("Note, if ≤10 SNPs are identified, consider investigating further."),
      tags$h4("NJ Phylogeny for Run"), 
      plotOutput(outputId = "p.NJ.tree"),
      hr(),
      tags$h4("Contextualised NJ Tree"), 
      column(
        width = 10, uiOutput("p.contextual.tree.ui")
      ),
      #tags$h4("Contextualised NJ Tree"), 
      #plotOutput(outputId = "p.contextual.tree"),
      hr(),
      tags$h3("Resistance/Lineage Summary"),
      plotOutput("p.Lineage.Resistance.bars"),
      hr(),
      tags$h3("Sample Report"),
      dataTableOutput("t.resistancetable"),
      hr()
      
    )
  )
)



# Define Server
server <- function(input, output, session) {
  # Reactive variable to hold selected files
  selectedFiles <- reactive({
    input$selected_files
  })
  
  # Set up shinyFiles to browse directories
  shinyDirChoose(input, "directory", roots = c(home = "~"), filetypes = c("", "txt"))
  
  # Reactive value to store the previously selected directory
  previousDirectory <- reactiveVal(NULL)
  
  # Flag to indicate whether a folder is selected
  folderSelected <- reactiveVal(FALSE)
  
  observe({
    cat("Directory selection changed\n")
    
    # Check if a directory is chosen
    if (is.null(input$directory) || length(input$directory) == 0) {
      cat("No directory selected or directory input is empty\n")
      output$selectedDirectory <- renderText("No directory selected")
      folderSelected(FALSE)  # Set the flag to FALSE when no directory is selected
      return()
    }
    
    # Get the path of the selected directory
    directoryPath <- try(parseDirPath(c(home = "~"), input$directory), silent = TRUE)
    cat("Parsed directory path:", directoryPath, "\n")
    
    # Check if directoryPath is valid
    if (inherits(directoryPath, "try-error") || is.null(directoryPath) || length(directoryPath) == 0) {
      cat("Directory path is invalid or empty\n")
      output$selectedDirectory <- renderText("Invalid directory path")
      return()
    }
    
    # Display the selected directory path
    output$selectedDirectory <- renderText({
      paste(directoryPath)
    })
    
    # Append the fixed subdirectory path
    subdirectoryPath <- file.path(directoryPath, "mapped_reads/")
    cat("Subdirectory path:", subdirectoryPath, "\n")
    
    # Check if subdirectoryPath exists
    if (!isTRUE(dir.exists(subdirectoryPath))) {
      cat("Subdirectory does not exist\n")
      output$checkboxes <- renderUI({
        h4("Subdirectory does not exist")
      })
      folderSelected(TRUE)  # Set the flag to TRUE since a valid directory is selected
      return()
    }
    
    # List files in the subdirectory and sort alphabetically
    files <- sort(list.files(subdirectoryPath))
    cat("Files in subdirectory:", paste(files, collapse = ", "), "\n")
    
    # Apply regex to trim filenames (remove '_sorted.bam' extension)
    trimmedFiles <- gsub("\\_sorted\\.bam$", "", files)
    
    # Create checkboxes dynamically
    current.checkbox_list <- checkboxGroupInput("selected_files", "\nSelect Files:", choices = trimmedFiles, selected = trimmedFiles)
    output$checkboxes <- renderUI({
      current.checkbox_list
    })
    #})
    
    # Check if the directory has changed
    if (!is.null(previousDirectory()) && previousDirectory() == directoryPath) {
      cat("Directory has not changed, skipping resource path addition\n")
      folderSelected(TRUE)  # Set the flag to TRUE since a valid directory is selected
      return()
    }
    
    # Update the previous directory
    previousDirectory(directoryPath)
    
    # Serve the directory containing the MultiQC report
    cat("Directory path for resource:", directoryPath, "\n")
    if (dir.exists(directoryPath)) {
      addResourcePath("multiqc", directoryPath)  # Map URL path to local directory
      cat("Resource path added for:", directoryPath, "\n")
    }
    folderSelected(TRUE)  # Set the flag to TRUE when a valid directory is selected and processed
  })
  
  # Display selected files
  output$selectedFiles <- renderText({
    if (is.null(selectedFiles()) || length(selectedFiles()) == 0) {
      return("No files selected")
    }
    paste("", paste(selectedFiles(), collapse = ", "))
  })
  
  # Reactive expression to capture MultiQC report path
  multiQC <- reactive({
    req(input$directory)  # Ensure the directory input is not null
    multiQCPath <- file.path(parseDirPath(c(home = "~"), input$directory), "multiqc/multiqc_report.html")
    cat("MultiQC html path:", multiQCPath, "\n")
    multiQCPath
  })
  
  # Render a link to the MultiQC report
  output$multiqcLink <- renderUI({
    req(input$directory)  # Ensure the directory input is not null
    multiQCPath <- multiQC()
    cat("MultiQC report link being generated\n")
    if (!file.exists(multiQCPath)) {
      cat("MultiQC report not found at:", multiQCPath, "\n")
      return(h4("MultiQC report not found"))
    }
    reportUrl <- file.path("multiqc", "multiqc/multiqc_report.html")  # Construct the URL
    cat("MultiQC report URL:", reportUrl, "\n")
    tags$a(href = reportUrl, target = "_blank", "Open MultiQC Report")
  })
  
  # Output to indicate whether a folder is selected
  output$folderSelected <- reactive({
    folderSelected()
  })
  outputOptions(output, "folderSelected", suspendWhenHidden = FALSE)
  
  ####################
  # Make some plots
  
  # First check if different processes have finished
  pipelineStatus <- reactive({
    req(input$selected_files)
    
    pipeline.status.df <- data.frame(Samples=input$selected_files, raw.mapping="yes")
    
    # check if reads are present
    fastq.dir <- file.path(parseDirPath(c(home = "~"), input$directory), "fastqs")
    cat("\ninternal fastq path:", fastq.dir, "\n")
    if (dir.exists(fastq.dir)){
      pipeline.status.df <- pipeline.status.df %>% left_join( data.frame(Samples=(sort(gsub("\\.fastq\\.gz","", list.files(fastq.dir)))), fastq="yes"), by="Samples") %>%
        replace_na(list(fastq = "no"))
    }
    else {
      pipeline.status.df <-  pipeline.status.df %>% left_join( data.frame(Samples=input$selected_files, fastq = "no"), by="Samples")
    }
    # Are post filter QC files available - readlengths
    readlength.dir <- file.path(parseDirPath(c(home = "~"), input$directory), "qc/post_filter_qc/readlengths/") 
    if (dir.exists(readlength.dir)){
      pipeline.status.df <- pipeline.status.df %>% left_join( data.frame(Samples=(sort(gsub("\\.read\\-lengths\\.tsv","", list.files(readlength.dir)))), post.qc.readlengths="yes"), by="Samples") %>%
        replace_na(list(post.qc.readlengths = "no")) }
    else {
      pipeline.status.df <-  pipeline.status.df %>% left_join( data.frame(Samples=input$selected_files, post.qc.readlengths = "no"), by="Samples")
    }
    # Are coverage summary files available
    cov.dir <- file.path(parseDirPath(c(home = "~"), input$directory), "qc/post_filter_qc/coverage/coverage_summary/") 
    if (dir.exists(cov.dir)){
      pipeline.status.df <- pipeline.status.df %>% left_join( data.frame(Samples=(sort(gsub("\\_coverage\\_summary\\.tsv","", list.files(cov.dir)))), coverage.summary="yes"), by="Samples") %>%
        replace_na(list(coverage.summary = "no"))
    }
    else {
      pipeline.status.df <- pipeline.status.df %>% left_join( data.frame(Samples=input$selected_files, coverage.summary = "no"), by="Samples")
    }
    # whether on target stats are available
    ontarget.file <- file.path(parseDirPath(c(home = "~"), input$directory), "qc/post_filter_qc/on_and_off_target_stats.csv")
    if (file.exists(ontarget.file)){
      Nextflow.mapping.stats1 <- read.csv(ontarget.file)
      pipeline.status.df <- pipeline.status.df %>% left_join( data.frame(Samples=Nextflow.mapping.stats1$Name, mapping.stats="yes"), by="Samples") %>%
        replace_na(list(mapping.stats = "no"))
    }
    else {
      pipeline.status.df <- pipeline.status.df %>% left_join( data.frame(Samples=input$selected_files, mapping.stats="no"), by="Samples")
    }
    # whether variants were called
    vars.dir <- file.path(parseDirPath(c(home = "~"), input$directory), "variants")
    if (dir.exists(vars.dir)){
      pipeline.status.df <- pipeline.status.df %>% left_join( data.frame(Samples=(sort(gsub("_clair3.gvcf.gz","", list.files(vars.dir)))), variants.file="yes"), by="Samples") %>%
        replace_na(list(variants.file = "no"))
    }
    else {
      pipeline.status.df <- pipeline.status.df %>% left_join( data.frame(Samples=input$selected_files, variants.file = "no"), by="Samples")
    }
    # whether a consensus fasta was made
    consensus.dir <- file.path(parseDirPath(c(home = "~"), input$directory), "curated_consensus")
    if (dir.exists(consensus.dir)){
      pipeline.status.df <- pipeline.status.df %>% left_join( data.frame(Samples=gsub("\\.fasta","",sort(list.files(consensus.dir)[grep("multi_locus",list.files(consensus.dir), invert=T)])), fasta.file="yes"), by="Samples") %>%
        replace_na(list(fasta.file = "no"))
    }
    else {
      pipeline.status.df <- pipeline.status.df %>% left_join( data.frame(Samples=input$selected_files, fasta.file = "no"), by="Samples")
    }
    
    cat("\nshow pipeline.status.df:")
    print(pipeline.status.df)
    pipeline.status.df
  })
  
  pipelineStatusMelt <- reactive({
    req(input$selected_files)
    req(pipelineStatus)
    # melt to long form for plotting
    pipeline.status.df.melt <- pipelineStatus() %>%
      pivot_longer(-Samples, names_to="Process", values_to="process.done") %>%
      mutate(Process=factor(Process, levels=c("fastq", "raw.mapping","mapping.stats","post.qc.readlengths","coverage.summary","variants.file","fasta.file"))) %>%
      mutate(Process.done=factor(process.done, levels=c("yes","no")))
    
    cat("show pipeline.status.df.melt:")
    print(pipeline.status.df.melt)
    pipeline.status.df.melt
    
  })
  
  # Now plot pipeline status with emojis
  output$p.pipeline.status <- renderPlot({
    # make plot
    p.pipeline.status <- pipelineStatusMelt() %>%
      ggplot(aes(y=Samples, x=Process, fill=Process.done)) +
      geom_tile(color='grey95', alpha=0.5, size=1.5) +
      theme_bw() + theme.text.size + legend.size + x.theme.axis.rotate.angle +
      scale_x_discrete(expand = c(0, 0)) +
      scale_fill_manual(values=c("green3", "red1"), breaks=c("yes","no")) +
      #geom_emoji(aes(image = ifelse(process.done=="yes", '1f600', '1f622'))) + # this looks fun, but takes ages to load
      theme(legend.position='top')
    p.pipeline.status
  })
  
  
  # Reactive data for read length distributions
  collatedLengths <- reactive({
    req(input$selected_files)
    
    readlength.directory <- file.path(parseDirPath(c(home = "~"), input$directory), "qc/post_filter_qc/readlengths/")
    cat("\nChecking read length distributions\n")
    cat("readlength path:", readlength.directory, "\n")
    
    collated.lengths <- NULL
    for (current.sample in input$selected_files) {
      current.length.sample <- paste0(readlength.directory, current.sample, ".read-lengths.tsv")
      if (file.exists(current.length.sample)){
        current.lengths <- read.table(current.length.sample, col.names = c("read.length", "count"))
        current.lengths$sample <- current.sample
      }
      else{
        cat("\nFile", current.length.sample, "does not exist.\n")
      }
      collated.lengths <- rbind(collated.lengths, current.lengths)
    }
    collated.lengths
  })
  
  output$p.readlengths <- renderPlot({
    p.readlengths <- collatedLengths() %>%
      arrange(sample) %>%
      group_by(sample) %>%
      mutate(total.reads = sum(count)) %>%
      ggplot(aes(x = read.length, y = sample, height = count, fill = sample)) +
      geom_density_ridges(stat = 'identity', scale = 1, linewidth = 0.35, alpha = 0.9) +
      theme_bw() + x.theme.axis.rotate.angle + #theme.text.size + legend.size +
      coord_cartesian(xlim=c(400,850)) +
      labs(y = "Sample", x = "Read Length") + theme(legend.position = 'none') +
      labs(title = "Read Length Distributions")
    p.readlengths
  })
  
  # On target mapping (no loop here)
  Nextflow.mapping.stats1.file <- reactive({
    file.path(parseDirPath(c(home = "~"), input$directory), "qc/post_filter_qc/on_and_off_target_stats.csv")
  })
  
  Nextflow.mapping.stats1 <- reactive({
    read.csv(Nextflow.mapping.stats1.file())
  })
  
  output$p.ontargetmapping <- renderPlot({
    p.on.target.mapping <- Nextflow.mapping.stats1() %>%
      arrange(Name) %>%
      select(Name, On.target.count, Off.target.count, On.target.percentage) %>%
      filter(Name %in% selectedFiles()) %>%
      pivot_longer(-c(Name, On.target.percentage), names_to = "Reads", values_to = "Count") %>%
      mutate(On.target.percentage = ifelse(Reads == "On.target.count", On.target.percentage, NA)) %>%
      ggplot(aes(y = Name, x = Count, fill = Reads)) +
      geom_barh(stat = 'identity', position = 'stack', width = 0.6) +
      theme_bw() + theme.text.size + legend.size + x.theme.axis.rotate.angle +
      geom_text(aes(y = Name, x = Count + 1500, label = On.target.percentage), angle = 0, size = text.size.within, hjust = 0) +
      coord_cartesian(x = c(0, max(Nextflow.mapping.stats1()$On.target.count) + 3000)) +
      scale_fill_brewer(palette = 'Dark2') + 
      scale_x_continuous(breaks = pretty, labels = scales::comma) +
      labs(y = "Sample", x = "Read Count", title = "Reads mapping to target regions")
    p.on.target.mapping
  })
  
  # Reactive data for per sample amplicon coverage
  collatedCov <- reactive({
    req(input$selected_files)
    
    samplecov.directory <- file.path(parseDirPath(c(home = "~"), input$directory), "qc/post_filter_qc/coverage/coverage_summary/")
    cat("\nChecking Sample Coverage\n")
    cat("samplecov path:", samplecov.directory, "\n")
    
    collated.cov <- NULL
    for (current.sample in input$selected_files) {
      current.cov.sample <- paste0(samplecov.directory, current.sample, "_coverage_summary.tsv")
      cat("\n current_samplecov file:", current.cov.sample)
      if (file.exists(current.cov.sample))
        current.cov <- read.table(current.cov.sample, header = TRUE)
      else
        cat("\nFile", current.cov.sample, "does not exist\n")
      collated.cov <- rbind(collated.cov, current.cov)
    }
    collated.cov
  })
  
  output$p.median.cov <- renderPlot({
    p.per_region.mediancov.heatmap <- collatedCov() %>% 
      select(sample, name, depth_median, start) %>%
      arrange(sample, start) %>%
      mutate(name = factor(name, levels=unique(name))) %>%
      ggplot(aes(x = name, y = sample, fill = depth_median)) +
      geom_tile(color = 'grey95') +
      theme_bw() + x.theme.axis.rotate.angle + theme.text.size + legend.size +
      scale_fill_viridis_b(option = 'D', trans = "log10", breaks = c(10, 25, 50), direction = -1, na.value = 'grey95') +
      labs(x = "Amplicon", y = "Sample", fill = "Median\nCoverage") +
      labs(title = "Median Coverage (X) per sample & amplicon")
    p.per_region.mediancov.heatmap
  })
  
  output$p.min.cov <- renderPlot({
    p.per_region.mincov.heatmap <- collatedCov() %>% 
      select(sample, name, depth_min, start) %>%
      arrange(sample, start) %>%
      mutate(name = factor(name, levels=unique(name))) %>%
      ggplot(aes(x = name, y = sample, fill = depth_min)) +
      geom_tile(color = 'grey95') +
      theme_bw() + x.theme.axis.rotate.angle + theme.text.size + legend.size +
      scale_fill_viridis_b(option = 'D', trans = "log10", breaks = c(10, 25, 50), direction = -1, na.value = 'grey95') +
      labs(x = "Amplicon", y = "Sample", fill = "Minimum\nCoverage") +
      labs(title = "Minimum Coverage (X) per sample & amplicon")
    p.per_region.mincov.heatmap
  })
  
  output$p.10x.cov <- renderPlot({
    p.per_region.10xcov.heatmap <- collatedCov() %>% 
      select(sample, name, cov_perc_10.0x, start) %>%
      mutate(cov_perc_below_10.0x = 100 - cov_perc_10.0x) %>%
      arrange(sample, start) %>%
      mutate(name = factor(name, levels=unique(name))) %>%
      ggplot(aes(x = name, y = sample, fill = cov_perc_below_10.0x)) +
      geom_tile(color = 'grey95') +
      theme_bw() + x.theme.axis.rotate.angle + theme.text.size + legend.size +
      scale_fill_viridis_b(option = 'D', trans = "log10", direction = 1, na.value = 'grey95') +
      labs(x = "Amplicon", y = "Sample", fill = "% Sites\n<10x") +
      labs(title = "% sites in amplicon <10X coverage")
    p.per_region.10xcov.heatmap
  })
  
  # Prepare and plot a basic phylogeny of SNPs
  
  
  collatedPhylo <- reactive({
    req(input$selected_files)
    
    # Read in concatenated SNPs (in fasta format)
    # Specify path for data
    multi_locus.filepath <- file.path(parseDirPath(c(home = "~"), input$directory), "curated_consensus/") 
    multi_locus.files <- paste0(multi_locus.filepath, input$selected_files, "_multi_locus.fasta", sep="")
    
    cat("\nChecking multi-locus directory\n")
    cat("Multi-locus file path:", multi_locus.filepath, "\n")
    cat("Multi-locus file list:", multi_locus.files, "\n")
    
    # read fasta files into a list
    multi_locus.sequences <- lapply(multi_locus.files, read.dna, format = "fasta")
    # convert into an alignment
    multi_locus.sequences_alignment <- do.call("rbind", multi_locus.sequences)
    # Update fasta headers
    AmpliSeq.full.multi_locus.fasta <- updateLabel(multi_locus.sequences_alignment, labels(multi_locus.sequences_alignment), gsub("\\_multi\\_locus\\ joined.+$","",labels(multi_locus.sequences_alignment)))
    
    # Subset tip labels to the selection
    AmpliSeq.full.multi_locus.selected <- AmpliSeq.full.multi_locus.fasta[labels(AmpliSeq.full.multi_locus.fasta) %in% input$selected_files,]
    
    # Convert to a phyDat object, calculate a distance matrix using phangorn, then infer an NJ tree 
    AmpliSeq.full.multi_locus.selected.phydat <- phyDat(AmpliSeq.full.multi_locus.selected, type = "DNA", levels = NULL)
    cat("\nMaking NJ tree\n")
    AmpliSeq.full.multi_locus.selected.dna_dist <- dist.ml(AmpliSeq.full.multi_locus.selected.phydat, model="JC69")
    AmpliSeq.full.multi_locus.selected.NJ <- NJ(AmpliSeq.full.multi_locus.selected.dna_dist)
    midpoint(AmpliSeq.full.multi_locus.selected.NJ)
    #cat("\nMaking ML tree\n")
    #AmpliSeq.full.multi_locus.selected.phydat_fitGTR <- pml_bb(AmpliSeq.full.multi_locus.selected.phydat, model="GTR+G(4)+I")
    #midpoint(AmpliSeq.full.multi_locus.selected.phydat_fitGTR$tree)
  })
  
  output$p.NJ.tree <- renderPlot({
    # Plot tree
    options(ignore.negative.edge=TRUE)
    p.AmpliSeq.full.selectedNJ <- ggtree(collatedPhylo() ) +
      geom_tiplab(size=text.size.within) +
      #coord_cartesian(xlim=c(0,max(ggtree(collatedPhylo() )$data$x)+2)) +
      coord_cartesian(xlim=c(0,max(ggtree(collatedPhylo() )$data$x)+0.0001)) + 
      geom_treescale(fontsize=text.size.within)
    p.AmpliSeq.full.selectedNJ
  })
  
  # Reactive element to get the number of variant positions in the alignment used for making the tree
  SNPslength <- reactive({
    req(input$selected_files)
    # Read in concatenated SNPs (in fasta format)
    SNP.concat.directory <- file.path(parseDirPath(c(home = "~"), input$directory), "snp_aln/")
    cat("\nChecking number of sites in SNP alignment\n")
    cat("SNP file path:", SNP.concat.directory, "\n")
    
    AmpliSeq.full.fasta.file <- paste0(SNP.concat.directory, "merged.fasta.snp.aln")
    cat("SNP file:", AmpliSeq.full.fasta.file)
    AmpliSeq.full.fasta.dnaBin <- read.dna(AmpliSeq.full.fasta.file, 'fasta')
    length(as.character(AmpliSeq.full.fasta.dnaBin)[1,]) # all sequences are the same length, so just get the length of the first
  })
  
  # Display selected files
  output$SNPcount <- renderText({
    if (is.null(selectedFiles()) || length(selectedFiles()) == 0) {
      return("No files selected")
    }
    paste("There were ", SNPslength(), " SNPs identified in the dataset.", collapse = "")
  })
  
  
  ## Now add contextual data
  
  # Specify contextual data
  AmpliSeq_contextual.fasta.file <- "/Users/mb29/Treponema/Treponema_Discriminatory_Sites__MinION/nextflow_pipeline_example_run_20240510/MAGUS_context_treemer0.4.multilocus.concat.aln"
  cat("\nContextual fasta sequences:",AmpliSeq_contextual.fasta.file,"\n")
  
  # Read in contextual sequence data
  ContextualSeqs <- reactive({
    (read.dna(AmpliSeq_contextual.fasta.file, 'fasta'))
  })
  
  
  # read in contextual data
  ContextualisedTree <- reactive({
    req(input$selected_files)
    req(ContextualSeqs())
    
    # Read in concatenated SNPs (in fasta format)
    # Specify path for data
    multi_locus.filepath <- file.path(parseDirPath(c(home = "~"), input$directory), "curated_consensus/") 
    multi_locus.files <- paste0(multi_locus.filepath, input$selected_files, "_multi_locus.fasta", sep="")
    
    cat("\nChecking multi-locus directory\n")
    cat("\nMulti-locus file path:", multi_locus.filepath, "\n")
    cat("Multi-locus file list:", multi_locus.files, "\n")
    
    # read fasta files into a list
    multi_locus.sequences <- lapply(multi_locus.files, read.dna, format = "fasta")
    # convert into an alignment
    multi_locus.sequences_alignment <- do.call("rbind", multi_locus.sequences)
    # Update fasta headers
    AmpliSeq.full.multi_locus.fasta <- updateLabel(multi_locus.sequences_alignment, labels(multi_locus.sequences_alignment), gsub("\\_multi\\_locus\\ joined.+$","",labels(multi_locus.sequences_alignment)))
    
    # Subset tip labels to the selection
    AmpliSeq.full.multi_locus.selected <- AmpliSeq.full.multi_locus.fasta[labels(AmpliSeq.full.multi_locus.fasta) %in% input$selected_files,]
    
    # Combine contextual and current sequence data
    cat("\nCombining new run seqs with contextual seqs\n")
    AmpliSeq_contextual_and_selected.dnabin <- rbind(AmpliSeq.full.multi_locus.selected, ContextualSeqs())
    
    # Convert to a phyDat object, calculate a distance matrix using phangorn, then infer an NJ tree 
    AmpliSeq_contextual_and_selected.phydat <- phyDat(AmpliSeq_contextual_and_selected.dnabin, type = "DNA", levels = NULL)
    
    # import current and contextual data, then make a NJ tree and infer lineages
    # Fit data using NJ and return tree
    cat("\nCalculating Tree for Current+Contextual sequences\n")
    AmpliSeq_contextual_and_selected_fitNJ <- dist.ml(AmpliSeq_contextual_and_selected.phydat, model="JC69")
    AmpliSeq_contextual_and_selected.phydat.NJ <- NJ(AmpliSeq_contextual_and_selected_fitNJ)
    AmpliSeq_contextual_and_selected_tree <- midpoint(AmpliSeq_contextual_and_selected.phydat.NJ)
    cat("\nOutput contextual tree\n")
    AmpliSeq_contextual_and_selected_tree
  })
  
  InferredLineages <- reactive({
    req(input$selected_files)
    req(ContextualisedTree())
    cat("\nExtracting metadata from contextual sequence headers\n")
    
    # contextual.metadata <- data.frame(sample=labels(ContextualSeqs())) %>%
    #   mutate(Lineage=gsub("^.+__","",sample)) %>%
    #   mutate(Country= gsub("^.+__","", gsub("__SS14","",gsub("__Nichols","",sample))))
    
    contextual.metadata <- data.frame(sample=ContextualisedTree()$tip.label) %>%
      filter(grepl('Nichols|SS14', sample)) %>%
      mutate(Lineage=gsub("^.+__","",sample)) %>%
      mutate(Country= gsub("^.+__","", gsub("__SS14","",gsub("__Nichols","",sample))))
    
    cat("\nShow contextual metadata:\n")
    print(contextual.metadata)
    #contextual.metadata
    
    # Infer Lineages (Nichols/SS14) for novel samples using MRCA/Descendents phylogenetic method
    inferred_Nichols.list <- data.frame(sample= ContextualisedTree()$tip.label[phangorn::Descendants(ContextualisedTree(), phangorn::mrca.phylo(ContextualisedTree(), filter(contextual.metadata, Lineage=="Nichols") %>% pull(sample)))[[1]] ],
                                        Lineage="Nichols")
    inferred_SS14.list <- data.frame(sample= ContextualisedTree()$tip.label[phangorn::Descendants(ContextualisedTree(), phangorn::mrca.phylo(ContextualisedTree(), filter(contextual.metadata, Lineage=="SS14") %>% pull(sample)))[[1]] ],
                                     Lineage="SS14")
    inferred_lineages <- data.frame(rbind(inferred_Nichols.list, inferred_SS14.list))
    
    cat("\nInferred Lineages:\n")
    print(inferred_lineages)
    
    inferred_lineages
  })
  

  output$p.contextual.tree <- renderPlot({
    req(input$selected_files)
    req(ContextualSeqs())
    req(ContextualisedTree())
    cat("\nCreating initial contextual tree object\n")
    # Plot tree
    p.AmpliSeq_contextual_and_selected_tree <- ggtree(ContextualisedTree(), ladderise='right') +
      #coord_cartesian(xlim=c(0,max(ggtree(fitGTR.tree)$data$x)+2)) +
      coord_cartesian(xlim=c(0,max(ggtree(ContextualisedTree(), ladderise='right')$data$x)+0.00030)) +
      geom_treescale(fontsize=text.size.within)
    cat("\nAdding coloured tips (inferred from current dataset)\n")
    p.AmpliSeq_contextual_and_selected_tree <- p.AmpliSeq_contextual_and_selected_tree %<+% data.frame(seq=input$selected_files, study="current") +
      geom_tiplab(aes(color = factor(study)), size=text.size.within, align = T, offset=0.000015) + 
      scale_color_manual(breaks=c("current"), values=c("green4"), na.value = "grey5", name="Current\nSequencing\nRun") +
      new_scale_color()
    cat("\nNow trying to add inferred Lineage data\n")
     p.AmpliSeq_contextual_and_selected_tree <- p.AmpliSeq_contextual_and_selected_tree %<+% InferredLineages() +
       geom_tippoint(aes(color=factor(Lineage)), alpha=0.75, size=4) +
       scale_color_manual(breaks=c("Nichols","SS14"), values=c("royalblue2", "indianred1"), name="Lineage")

    p.AmpliSeq_contextual_and_selected_tree
    
  }, height = 700, width = 550 )
  
  
  output$p.contextual.tree.ui <- renderUI({
    plotOutput("p.contextual.tree", height = 700)
  })
  
  
  ResistanceTable <- reactive({
    req(input$selected_files)

    variants.filepath <- file.path(parseDirPath(c(home = "~"), input$directory), "variants/merged_gvcf/")
    # Function to get the most recent file with suffix _merged.tsv
    get_most_recent_file <- function(directory, suffix) {
      files <- list.files(directory, pattern = paste0(".*", suffix, "$"), full.names = TRUE)
      if (length(files) == 0) {
        return(NULL)
      }
      files_info <- file.info(files)
      most_recent_file <- rownames(files_info)[which.max(files_info$mtime)]
      return(most_recent_file)
    }
    latest_variants.file <- get_most_recent_file(variants.filepath, "_merged.tsv")
    latest_variants <- read.csv(latest_variants.file, sep=" ", col.names = c("Reference","POS","REF.allele", "ALT.alleles","SampleID", "GT", "GT_allele"), header = F)
    # Filter to only include selected samples 
    latest_variants.selected <- latest_variants %>% filter(SampleID %in% input$selected_files)
    
    # Clean up and summarise
    latest_variants.selected <- latest_variants.selected %>%
      filter(POS %in% c(235246)) %>%
      mutate(ResistanceSite="A2058") %>%
      select(SampleID, ResistanceSite, GT_allele) %>%
      mutate(Resistant = ifelse(GT_allele=="G", "Resistant", "Sensitive")) %>%
      arrange(SampleID) %>%
      # combine with lineage information inferred earlier
      left_join(InferredLineages(), by=c("SampleID"="sample"))
    latest_variants.selected
    cat("Compile 23S variants into a table")
    print(latest_variants.selected)
  })
  
  output$p.Lineage.Resistance.bars <- renderPlot({
    req(input$selected_files)
    req(ResistanceTable())
    # Prepare macrolide resistance bar plot
    p.macrolide.Res.bar <- ResistanceTable() %>%
      mutate(total.samples=n()) %>%
      group_by(Resistant) %>%
      mutate(Res.Count=n(), perc.Resistant=round((Res.Count/total.samples)*100,2)) %>%
      distinct(Res.Count, perc.Resistant) %>%
      ggplot(aes(x=Resistant, y=Res.Count, fill=Resistant)) +
      geom_bar(stat='identity', width=0.6) +
      theme_minimal() + 
      x.theme.axis.rotate.angle + theme.text.size + legend.size +
      geom_text(aes(x=Resistant, y=Res.Count+1, label = paste(perc.Resistant,"%")), size=text.size.within, inherit.aes = F) +
      scale_fill_manual(values=c("grey5", "grey85"), breaks=c("Resistant","Sensitive")) +
      labs(y="Sample Count", x="A2058G Macrolide Resistance") + theme(legend.position='none')
    # Now prepare Lineage summary plot
    p.Lineage.bar <- ResistanceTable() %>%
      mutate(total.samples=n()) %>%
      group_by(Lineage) %>%
      mutate(Lineage.Count=n(), perc.Lineage=round((Lineage.Count/total.samples)*100,2)) %>%
      distinct(Lineage.Count, perc.Lineage) %>%
      ggplot(aes(x=Lineage, y=Lineage.Count, fill=Lineage)) +
      geom_bar(stat='identity', width=0.6) +
      #theme_bw() + 
      theme_minimal() +
      x.theme.axis.rotate.angle + theme.text.size + legend.size +
      geom_text(aes(x=Lineage, y=Lineage.Count+1, label = paste(perc.Lineage,"%")), size=text.size.within, inherit.aes = F) +
      scale_fill_manual(breaks=c("Nichols","SS14"), values=c("royalblue2", "indianred1"), name="Lineage") +
      labs(y="Sample Count", x="Lineage") + theme(legend.position='none')
    # Now make combined figure using cowplot
    plot_grid(p.macrolide.Res.bar, p.Lineage.bar, ncol=2, labels=c("Macrolide Resistance", "Lineage"), label_size = 11, scale=0.95)
  })
    
    
    
  


  output$t.resistancetable <- renderDataTable({
    req(input$selected_files)
    ResistanceTable()
  })
  
  
  
  
  
}


# Run the application 
shinyApp(ui = ui, server = server)

