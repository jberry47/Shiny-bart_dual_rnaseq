library(reshape2)
library(grid)
library(stringr)
library(scales)
library(gridExtra)
library(shiny)
library(shinyBS)
library(shinydashboard)
library(DT)
library(shinyjs)
library(plyr)
library(ggplot2)
library(openxlsx)

sorg_data <- list(
  ntj2_mi = read.table("data/cuffdiff/cuffdiff_classicfpkm_dispersionpercondition_dualRNAseq_ntj2_m.x_Sbicolor_454_v3.0.1_xhan_rnaseq1_s.annot_190118/gene_exp.diff",header=T,sep="\t",stringsAsFactors = F),
  bs_mi = read.table("data/cuffdiff/cuffdiff_classicfpkm_dispersionpercondition_dualRNAseq_bs_m.x_Sbicolor_454_v3.0.1_xhan_rnaseq1_s.annot_190118/gene_exp.diff",header=T,sep="\t",stringsAsFactors = F),
  grassl_mi = read.table("data/cuffdiff/cuffdiff_classicfpkm_dispersionpercondition_dualRNAseq_grassl_m.x_Sbicolor_454_v3.0.1_xhan_rnaseq1_s.annot_190118/gene_exp.diff",header=T,sep = "\t",stringsAsFactors = F),
  mock_ntj2_bs = read.table("data/cuffdiff/cuffdiff_classicfpkm_dispersionpercondition_dualRNAseq_bs.mock_ntj2.mock_Sbicolor_454_v3.0.1_xhan_rnaseq1_s.annot_190118/gene_exp.diff",header=T,sep = "\t",stringsAsFactors = F),
  mock_ntj2_grassl = read.table("data/cuffdiff_classicfpkm_dispersionpercondition_dualRNAseq_grassl.mock_ntj2.mock_Sbicolor_454_v3.0.1_xhan_rnaseq1_s.annot_190118/gene_exp.diff",header=T,sep = "\t",stringsAsFactors = F),
  mock_bs_grassl = read.table("data/cuffdiff/cuffdiff_classicfpkm_dispersionpercondition_dualRNAseq_bs.mock_grassl.mock_Sbicolor_454_v3.0.1_xhan_rnaseq1_s.annot_190118/gene_exp.diff",header=T,sep = "\t",stringsAsFactors = F),
  inf_ntj2_bs = read.table("data/cuffdiff/cuffdiff_classicfpkm_dispersionpercondition_dualRNAseq_ntj2.xanh_bs.xanh_Sbicolor_454_v3.0.1_xhan_rnaseq1_s.annot_190118/gene_exp.diff",header=T,sep = "\t",stringsAsFactors = F),
  inf_ntj2_grassl = read.table("data/cuffdiff/cuffdiff_classicfpkm_dispersionpercondition_dualRNAseq_ntj2.xanh_grassl.xanh_Sbicolor_454_v3.0.1_xhan_rnaseq1_s.annot_190118/gene_exp.diff",header=T,sep = "\t",stringsAsFactors = F),
  inf_bs_grassl = read.table("data/cuffdiff/cuffdiff_classicfpkm_dispersionpercondition_dualRNAseq_grassl.xanh_bs.xanh_Sbicolor_454_v3.0.1_xhan_rnaseq1_s.annot_190118/gene_exp.diff",header=T,sep = "\t",stringsAsFactors = F)
)

xanh_data <- list(
  xanh_bs = read.table("data/cuffdiff/cuffdiff_classicfpkm_dispersionpercondition_dualRNAseq_xanh_bs.xanh_Sbicolor_454_v3.0.1_xhan_rnaseq1_xanh.annot_190118/isoform_exp.diff",header=T,sep = "\t",stringsAsFactors = F),
  xanh_ntj2 = read.table("data/cuffdiff/cuffdiff_classicfpkm_dispersionpercondition_dualRNAseq_xanh_ntj2.xanh_Sbicolor_454_v3.0.1_xhan_rnaseq1_xanh.annot_190118/isoform_exp.diff",header=T,sep = "\t",stringsAsFactors = F),
  xanh_grassl = read.table("data/cuffdiff/cuffdiff_classicfpkm_dispersionpercondition_dualRNAseq_xanh_grassl.xanh_Sbicolor_454_v3.0.1_xhan_rnaseq1_xanh.annot_190118/isoform_exp.diff",header=T,sep = "\t",stringsAsFactors = F)
)

sorg_annot <- read.table("data/sorg_annot.tsv",header=T,stringsAsFactors = F,sep="\t",quote = "",fill=T)
xanh_annot <- read.csv("data/xanh_new anno.csv",stringsAsFactors = F)

load("data/ballgown/sorg_ballgown.RData")
xanh_fpkm <- read.csv("data/ballgown/xanh_ballgown.csv",stringsAsFactors = F)

ui <- dashboardPage(skin="black", title="Dual-RNAseq",
  dashboardHeader(
    title = tagList(
      tags$span(
        class = "logo-mini", "DRS"
      ),
      tags$span(
        class = "logo-lg", "Dual-RNAseq"
      )
    ),
    titleWidth = 450
  ),
  dashboardSidebar(
    tags$script(HTML("$('body').addClass('sidebar-mini');")),
    width = 150,
    sidebarMenu(
      shinyjs::useShinyjs(),
      id="tabs",
      hidden(menuSubItem("dummy", tabName = "dummy")),
      menuItem("Overview", tabName = "overview"),
      menuItem("Sorghum Reads",tabName = "sorg_filt",icon = icon("sitemap")),
      menuItem("Xanh Reads",tabName = "xanh_filt",icon = icon("sitemap"))
    )
  ),
  dashboardBody(
    tags$style(HTML("
            .tabbable > .nav > li[class=active]    > a {background-color: #444444; color:white}
            .multicol{
              -webkit-column-count: 4; /* Chrome, Safari, Opera */
              -moz-column-count: 4; /* Firefox */
              column-count: 4;
            }
            .twocol{
              -webkit-column-count: 2; /* Chrome, Safari, Opera */
              -moz-column-count: 2; /* Firefox */
              column-count: 2;
            }
            .warning { 
              color: red;
            }"
    )),
    uiOutput("login"),
    uiOutput("body")
  )
)

server = function(input, output, session){
  output$login <- renderUI({
    showModal(modalDialog(easyClose = TRUE,footer = NULL,
      passwordInput("password", "Password:"),
      actionButton("go", "Go")
    ))
  })
  
  cc_auth <- reactiveValues(data="NO")
  
  observeEvent(input$go,{
    isolate({
      pass <- input$password
    })
    if(pass == "test"){
      cc_auth$data <- "YES"
      #      obs1$suspend()
      removeModal()
      updateTabItems(session, "tabs", "overview")
    }else{
      cc_auth$data <- "NO"
    }
  })
  
  
  output$body <- renderUI({
    tags$head(tags$style("#container * {display: inline;}"))
    if(cc_auth$data == "YES"){
      tabItems(
        tabItem(tabName = "overview",
          fluidRow(
            box(style = "overflow-y:scroll",width=10,title = "Welcome",solidHeader = T,status = 'success',collapsible = TRUE,
              p("test")
            )
          )
        ),
        tabItem(tabName = "sorg_filt",
          fluidRow(
            box(style = "overflow-y:scroll",width=12,title = "Select data to filter",solidHeader = T,status = 'success',collapsible = F,
              column(width=12,
                column(width=3,
                  radioButtons("ntj2_mi_sel", "NTJ2: Mock vs Inf",c("Significant" = "sig", "NS"="ns","Either"="either"),selected = 'either'),
                  radioButtons("mock_ntj2_bs_sel", "Mock: NTJ2 vs BS",c("Significant" = "sig", "NS"="ns","Either"="either"),selected = 'either'),
                  radioButtons("inf_ntj2_bs_sel", "Inf: NTJ2 vs BS",c("Significant" = "sig", "NS"="ns","Either"="either"),selected = 'either')
                ),
                column(width=3,
                  radioButtons("bs_mi_sel", "BS: Mock vs Inf",c("Significant" = "sig", "NS"="ns","Either"="either"),selected = 'either'),
                  radioButtons("mock_ntj2_grassl_sel", "Mock: NTJ2 vs Grassl",c("Significant" = "sig", "NS"="ns","Either"="either"),selected = 'either'),
                  radioButtons("inf_ntj2_grassl_sel", "Inf: NTJ2 vs Grassl",c("Significant" = "sig", "NS"="ns","Either"="either"),selected = 'either')
                ),
                column(width=3,
                  radioButtons("grassl_mi_sel", "Grassl: Mock vs Inf",c("Significant" = "sig", "NS"="ns","Either"="either"),selected = 'either'),
                  radioButtons("mock_bs_grassl_sel", "Mock: BS vs grassl",c("Significant" = "sig", "NS"="ns","Either"="either"),selected = 'either'),
                  radioButtons("inf_bs_grassl_sel", "Inf: BS vs grassl",c("Significant" = "sig", "NS"="ns","Either"="either"),selected = 'either')
                ),
                column(width=3,
                  numericInput("fc_cut","Fold-change cutoff",value = 2,min = 0,max = Inf,step = 1),
                  numericInput("q_cut","Q-value cutoff",value = 0.01,min = 0,max = 1,step = 0.1)
                )
              ),
              actionButton("filter_genes","Filter Genes"),
              br(),br(),hr(),br(),
              dataTableOutput("dt_filtered_genes"),
              br(),
              uiOutput("filtered_genes_name_ui"),
              uiOutput("filtered_genes_download_ui"),
              br(),br(),
              tabsetPanel(
                tabPanel("Gene",
                  plotOutput("sel_gene_plot",width=800),
                  uiOutput("gene_boxplot_download_ui")
                ),
                tabPanel("Transcripts",
                  uiOutput("sel_transcript_ui"),
                  plotOutput("sel_transcript_boxplot",width=800),
                  uiOutput("transcript_boxplot_download_ui")
                )
              )

            )
          )
        ),
        tabItem(tabName = "xanh_filt",
          fluidRow(
            box(style = "overflow-y:scroll",width=12,title = "Select data to filter",solidHeader = T,status = 'success',collapsible = F,
              column(width=12,
                column(width=3,
                  radioButtons("xanh_bs_sel", "Xanh vs Xanh-BS",c("Significant" = "sig", "NS"="ns","Either"="either"),selected = 'either')
                ),
                column(width=3,
                  radioButtons("xanh_ntj2_sel", "Xanh vs Xanh-NTJ2",c("Significant" = "sig", "NS"="ns","Either"="either"),selected = 'either')
                ),
                column(width=3,
                  radioButtons("xanh_grassl_sel", "Xanh vs Xanh-Grassl",c("Significant" = "sig", "NS"="ns","Either"="either"),selected = 'either')
                ),
                column(width=3,
                  numericInput("xanh_fc_cut","Fold-change cutoff",value = 2,min = 0,max = Inf,step = 1),
                  numericInput("xanh_q_cut","Q-value cutoff",value = 0.01,min = 0,max = 1,step = 0.1)
                )
              ),
              actionButton("xanh_filter_genes","Filter Genes"),
              br(),br(),hr(),br(),
              dataTableOutput("xanh_dt_filtered_genes"),
              br(),
              uiOutput("xanh_filtered_genes_name_ui"),
              uiOutput("xanh_filtered_genes_download_ui"),
              br(),hr(),br(),
              uiOutput("xanh_transcript_ui"),
              plotOutput("xanh_gene_boxplot",width = 500),
              uiOutput("xanh_boxplot_download_ui")
            )
          )
        )
      )
    }
  })
  
  #********************************************************************************************************
  # SORGHUM READS
  #********************************************************************************************************
  filtered_genes <- eventReactive(input$filter_genes,{
    selections <- c(input$ntj2_mi_sel,input$bs_mi_sel,input$grassl_mi_sel,
                    input$mock_ntj2_bs_sel,input$mock_ntj2_grassl_sel,input$mock_bs_grassl_sel,
                    input$inf_ntj2_bs_sel,input$inf_ntj2_grassl_sel,input$inf_bs_grassl_sel)
    fc_cut <- input$fc_cut
    q_cut <- input$q_cut

    if(!all(selections == "either")){
      sig_sel <- c("ntj2_mi","bs_mi","grassl_mi","mock_ntj2_bs","mock_ntj2_grassl","mock_bs_grassl","inf_ntj2_bs","inf_ntj2_grassl","inf_bs_grassl")[which(selections == "sig")]
      ns_sel <- c("ntj2_mi","bs_mi","grassl_mi","mock_ntj2_bs","mock_ntj2_grassl","mock_bs_grassl","inf_ntj2_bs","inf_ntj2_grassl","inf_bs_grassl")[which(selections == "ns")]
      if(length(sig_sel)>0){
        sig <- unique(Reduce(intersect,lapply(sig_sel,function(i) with(sorg_data[[i]],gene[abs(log2.fold_change.) > fc_cut & q_value < q_cut]))))
      }else{sig <- NULL}
      
      if(length(ns_sel) >0){
        ns <- unique(Reduce(union,lapply(ns_sel,function(i) with(sorg_data[[i]],gene[abs(log2.fold_change.) > fc_cut & q_value < q_cut]))))
      }else{ns <- NULL}
      
      if(is.null(ns)){
        out <- sig
      }else if(is.null(sig)){
        out <- ns
      }else{
        out <- sig[!(sig %in% ns)]
      }
      sorg_fpkm[sorg_fpkm$gene_name %in% out,]
    }else{
      data.frame("Empty"="empty")
    }
  })
  
  output$dt_filtered_genes <- renderDataTable(server = T,{
    dat <- sorg_annot[sorg_annot$locusName %in% filtered_genes()$gene_name,]
    datatable(dat,rownames = F,selection="single")
  })
  
  output$filtered_genes_name_ui <- renderUI({
    if(length(input$dt_filtered_genes_rows_all)!=0){
      textInput("xlsx_filename","Filename to use on download","dual_rnaseq_filtered.xlsx",width=400)
    }
  })
  
  output$filtered_genes_download_ui <- renderUI({
    if(length(input$dt_filtered_genes_rows_all)!=0){
      downloadButton("filtered_genes_download","Download Genes (xlsx)")
    }
  })
  
  output$filtered_genes_download <- downloadHandler(
    filename = function() {input$xlsx_filename},
    content=function(file){
      my_list <- list(
        ballgown=filtered_genes(),
        annotation=sorg_annot[sorg_annot$locusName %in% filtered_genes()$gene_name,]
      )
      write.xlsx(my_list, file = file)
    }
  )
  
  gene_data <- eventReactive(input$dt_filtered_genes_row_last_clicked,{
    dat <- sorg_annot[sorg_annot$locusName %in% filtered_genes()$gene_name,]
    temp <- dat[input$dt_filtered_genes_row_last_clicked,"locusName"]
    my_gene <- melt(filtered_genes()[filtered_genes()$gene_name == temp,],id=c("chr","t_name","gene_id","gene_name"))
    my_gene
  })
  
  gene_boxplot <- eventReactive(input$dt_filtered_genes_row_last_clicked,{
    if(!is.null(input$dt_filtered_genes_row_last_clicked)){
      my_gene <- gene_data()
      meta <- strsplit(str_sub(my_gene$variable,28,40),"[.]")
      my_gene$genotype <- unlist(lapply(meta,function(i)i[1]))
      my_gene$inoc <- unlist(lapply(meta,function(i)i[2]))
      gene_name <- my_gene$gene_name[1]
      my_gene <- aggregate(data=my_gene,value~genotype+inoc,FUN="sum")
      my_gene$inoc <- gsub('[[:digit:]]+', '', my_gene$inoc)
      my_gene$genotype <- ordered(my_gene$genotype,levels=c("bs","grassl","ntj2"),labels=c("BS","Grassl","NTJ2"))
      my_gene$inoc <- ordered(my_gene$inoc,levels=c("mock","xanh"),labels=c("Mock","Xanh"))

      ggplot(my_gene,aes(inoc,value))+
        ggtitle(gene_name)+
        facet_grid(~genotype)+
        geom_boxplot()+
        ylab("FPKM")+
        xlab("")+
        theme_light()+
        theme(axis.text = element_text(size = 16),
          axis.title= element_text(size = 18),
          plot.title = element_text(size=22))+
        theme(plot.title = element_text(hjust = 0.5),
          strip.background=element_rect(fill="gray50"),
          strip.text.x=element_text(size=16,color="white"),
          strip.text.y=element_text(size=16,color="white"))
    }else{ggplot()}
  })
  
  output$sel_gene_plot <- renderPlot({
    gene_boxplot()
  })
  
  output$gene_boxplot_download_ui <- renderUI({
    if(!is.null(input$dt_filtered_genes_row_last_clicked)){
      downloadButton("gene_boxplot_download","Download Plot (png)")
    }
  })
  
  output$gene_boxplot_download <- downloadHandler(
    filename = function() {paste0(gene_data()$gene_name[1],"_boxplot.png")},
    content=function(file){
      ggsave(file,gene_boxplot(),device = "png",width = 6,height = 3.5,dpi = 300)
    })
  
  output$sel_transcript_ui <- renderUI({
    my_gene <- gene_data()
    my_transcripts <- sort(unique(my_gene$t_name))
    selectInput("which_trans","Select transcript", my_transcripts,my_transcripts[1],width = 200)
  })
  
  transcript_boxplot <- eventReactive(input$which_trans,{
    if(!is.null(input$dt_filtered_genes_row_last_clicked)){
      my_gene <- gene_data()
      my_gene <- my_gene[my_gene$t_name == input$which_trans,]
      meta <- strsplit(str_sub(my_gene$variable,28,40),"[.]")
      my_gene$genotype <- unlist(lapply(meta,function(i)i[1]))
      my_gene$inoc <- unlist(lapply(meta,function(i)i[2]))
      my_gene$inoc <- gsub('[[:digit:]]+', '', my_gene$inoc)
      my_gene$genotype <- ordered(my_gene$genotype,levels=c("bs","grassl","ntj2"),labels=c("BS","Grassl","NTJ2"))
      my_gene$inoc <- ordered(my_gene$inoc,levels=c("mock","xanh"),labels=c("Mock","Xanh"))
      
      ggplot(my_gene,aes(inoc,value))+
        ggtitle(input$which_trans)+
        facet_grid(~genotype)+
        geom_boxplot()+
        ylab("FPKM")+
        xlab("")+
        theme_light()+
        theme(axis.text = element_text(size = 16),
          axis.title= element_text(size = 18),
          plot.title = element_text(size=22))+
        theme(plot.title = element_text(hjust = 0.5),
          strip.background=element_rect(fill="gray50"),
          strip.text.x=element_text(size=16,color="white"),
          strip.text.y=element_text(size=16,color="white"))
    }else{ggplot()}
  })
  
  output$sel_transcript_boxplot <- renderPlot({
    transcript_boxplot()
  })
  
  output$transcript_boxplot_download_ui <- renderUI({
    if(!is.null(input$dt_filtered_genes_row_last_clicked)){
      downloadButton("transcript_boxplot_download","Download Plot (png)")
    }
  })
  
  output$transcript_boxplot_download <- downloadHandler(
    filename = function() {paste0(input$which_trans,"_boxplot.png")},
    content=function(file){
      ggsave(file,transcript_boxplot(),device = "png",width = 6,height = 3.5,dpi = 300)
    })
  
  #********************************************************************************************************
  # XANH READS
  #********************************************************************************************************
  filtered_genes_xanh <- eventReactive(input$xanh_filter_genes,{
    selections <- c(input$xanh_bs_sel,input$xanh_ntj2_sel,input$xanh_grassl_sel)
    fc_cut <- input$xanh_fc_cut
    q_cut <- input$xanh_q_cut
    
    if(!all(selections == "either")){
      sig_sel <- c("xanh_bs","xanh_ntj2","xanh_grassl")[which(selections == "sig")]
      ns_sel <- c("xanh_bs","xanh_ntj2","xanh_grassl")[which(selections == "ns")]
      if(length(sig_sel)>0){
        sig <- unique(Reduce(intersect,lapply(sig_sel,function(i) with(xanh_data[[i]],test_id[abs(log2.fold_change.) > fc_cut & q_value < q_cut]))))
      }else{sig <- NULL}
      
      if(length(ns_sel) >0){
        ns <- unique(Reduce(union,lapply(ns_sel,function(i) with(xanh_data[[i]],test_id[abs(log2.fold_change.) > fc_cut & q_value < q_cut]))))
      }else{ns <- NULL}
      
      if(is.null(ns)){
        out <- sig
      }else if(is.null(sig)){
        out <- ns
      }else{
        out <- sig[!(sig %in% ns)]
      }
      xanh_fpkm[xanh_fpkm$t_name %in% out,]
    }else{
      data.frame("Empty"="empty")
    }
  })
  
  output$xanh_dt_filtered_genes <- renderDataTable(server = T,{
    dat <- xanh_annot[xanh_annot$ID %in% filtered_genes_xanh()$t_name,]
    datatable(dat,rownames = F,selection="single")
  })
  
  output$xanh_filtered_genes_name_ui <- renderUI({
    if(length(input$xanh_dt_filtered_genes_rows_all)!=0){
      textInput("xanh_xlsx_filename","Filename to use on download","xanh_dual_rnaseq_filtered.xlsx",width=400)
    }
  })
  
  output$xanh_filtered_genes_download_ui <- renderUI({
    if(length(input$xanh_dt_filtered_genes_rows_all)!=0){
      downloadButton("xanh_filtered_genes_download","Download Genes (xlsx)")
    }
  })
  
  output$xanh_filtered_genes_download <- downloadHandler(
    filename = function() {input$xanh_xlsx_filename},
    content=function(file){
      my_list <- list(
        ballgown=filtered_genes_xanh(),
        annotation=xanh_annot[xanh_annot$ID %in% filtered_genes_xanh()$t_name,]
      )
      write.xlsx(my_list, file = file)
    }
  )
  
  xanh_gene_data <- eventReactive(input$xanh_dt_filtered_genes_row_last_clicked,{
    dat <- xanh_annot[xanh_annot$ID %in% filtered_genes_xanh()$t_name,]
    temp <- dat[input$xanh_dt_filtered_genes_row_last_clicked,"ID"]
    my_gene <- melt(filtered_genes_xanh()[filtered_genes_xanh()$t_name == temp,],id=c("chr","t_name","gene_id","gene_name"))
    my_gene
  })
  
  output$xanh_transcript_ui <- renderUI({
    my_gene <- xanh_gene_data()
    my_transcripts <- sort(unique(my_gene$t_name))
    selectInput("xanh_which_trans","Select transcript", my_transcripts,my_transcripts[1])
  })
  
  xanh_gene_boxplot <- eventReactive(input$xanh_which_trans,{
    if(!is.null(input$xanh_dt_filtered_genes_row_last_clicked)){
      my_gene <- xanh_gene_data()
      my_gene <- my_gene[my_gene$t_name == input$xanh_which_trans,]
      meta <- strsplit(str_sub(my_gene$variable,23,40),"[.]")
      my_gene$genotype <- unlist(lapply(meta,function(i)if(length(i)==1){"xanh"}else{i[2]}))
      my_gene$inoc <- unlist(lapply(meta,function(i)if(length(i)==1){i[1]}else{i[3]}))
      my_gene$inoc <- gsub('[[:digit:]]+', '', my_gene$inoc)
      my_gene <- my_gene[my_gene$inoc == "xanh",]
      my_gene$genotype <- sapply(my_gene$genotype,function(i)if(i != "xanh"){paste0("xanh+",i,collapse="")}else{"xanh"})
      my_gene$genotype <- ordered(my_gene$genotype,levels=c("xanh","xanh+bs","xanh+grassl","xanh+ntj2"),labels=c("Xanh","Xanh+BS","Xanh+Grassl","Xanh+NTJ2"))

      ggplot(my_gene,aes(genotype,value))+
        ggtitle(input$xanh_which_trans)+
        geom_boxplot()+
        ylab("FPKM")+
        xlab("")+
        theme_light()+
        theme(axis.text = element_text(size = 16),
          axis.title= element_text(size = 18),
          plot.title = element_text(size=22))+
        theme(plot.title = element_text(hjust = 0.5),
          strip.background=element_rect(fill="gray50"),
          strip.text.x=element_text(size=16,color="white"),
          strip.text.y=element_text(size=16,color="white"))+
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }else{ggplot()}
  })
  
  xanh_gene_name <- eventReactive(input$xanh_dt_filtered_genes_row_last_clicked,{
    filtered_genes_xanh()[input$xanh_dt_filtered_genes_row_last_clicked,"t_name"]
  })
  
  output$xanh_gene_boxplot <- renderPlot({
    xanh_gene_boxplot()
  })
  
  output$xanh_boxplot_download_ui <- renderUI({
    if(!is.null(input$xanh_dt_filtered_genes_row_last_clicked)){
      downloadButton("xanh_boxplot_download","Download Plot (png)")
    }
  })
  
  output$xanh_boxplot_download <- downloadHandler(
    filename = function() {paste0(input$xanh_which_trans,"_boxplot.png")},
    content=function(file){
      ggsave(file,xanh_gene_boxplot(),device = "png",width = 4,height = 5,dpi = 300)
    })
}

shinyApp(ui, server)
  