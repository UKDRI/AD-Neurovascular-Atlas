library(shiny)
library(shinyhelper)
library(data.table)
library(Matrix)
library(DT)
library(magrittr)
library(ggplot2)
library(ggrepel)
library(hdf5r)
library(ggdendro)
library(grid)
library(gridExtra)
library(qs)
library(tradeSeq)
library(ggsci)
library(SummarizedExperiment)
library(purrr)
library(dplyr)
library(ggpubr)
sc1conf = readRDS("./sc1conf.rds")
sc1meta = readRDS("./sc1meta.rds")
sc1gene = readRDS("./sc1gene.rds")
sc1dimr = readRDS("./sc1dimr.rds")
sc1def = readRDS("./sc1def.rds")
mast_de <- readr::read_csv("./mast_de_merged.csv")
all_pseudotime_sig_genes <- qs::qread("./tradeseq_sig_degs.qs")
# Get all unique genes
all_pseudotime_sig_genes <- all_pseudotime_sig_genes |>
  purrr::map(dplyr::select, hgnc_symbol, ensembl_gene_id) |>
  purrr::list_rbind() |>
  unique()
sce_subset <- qs::qread(
  "./slingshot_tradeseq_3k_case_and_control.qs"
)

cat(file = stderr(), "Finished reading in data objects")
counts <- assays(sce_subset)$counts
source("./shinyFunc.R")
pdf(file = NULL)


### Useful stuff
# Colour palette
cList = list(
  c(
    "grey85",
    "#FFF7EC",
    "#FEE8C8",
    "#FDD49E",
    "#FDBB84",
    "#FC8D59",
    "#EF6548",
    "#D7301F",
    "#B30000",
    "#7F0000"
  ),
  c(
    "#4575B4",
    "#74ADD1",
    "#ABD9E9",
    "#E0F3F8",
    "#FFFFBF",
    "#FEE090",
    "#FDAE61",
    "#F46D43",
    "#D73027"
  )[c(1, 1:9, 9)],
  c(
    "#FDE725",
    "#AADC32",
    "#5DC863",
    "#27AD81",
    "#21908C",
    "#2C728E",
    "#3B528B",
    "#472D7B",
    "#440154"
  )
)
names(cList) = c("White-Red", "Blue-Yellow-Red", "Yellow-Green-Purple")

# Panel sizes
pList = c("400px", "600px", "800px")
names(pList) = c("Small", "Medium", "Large")
pList2 = c("500px", "700px", "900px")
names(pList2) = c("Small", "Medium", "Large")
pList3 = c("600px", "800px", "1000px")
names(pList3) = c("Small", "Medium", "Large")
sList = c(18, 24, 30)
names(sList) = c("Small", "Medium", "Large")


### Start server code
shinyServer(function(input, output, session) {
  ### For all tags and Server-side selectize
  observe_helpers()
  optCrt = "{ option_create: function(data,escape) {return('<div class=\"create\"><strong>' + '</strong></div>');} }"
  observe({
    invalidateLater(30000) # ping every 30 seconds to keep connection alive
    cat(".")
  })

  ### Functions for tab A1
  getGsc1a1inp1 <- reactive({
    req(gsub("^Assay: ", "", input$sc1a1ass1))
    if (gsub("^Assay: ", "", input$sc1a1ass1) == "Cell Information") {
      res <- sc1conf$UI
      resDef <- sc1def$meta1
      resLen <- length(res)
    } else {
      res <- names(sc1gene[[gsub("^Assay: ", "", input$sc1a1ass1)]])
      resDef <- sc1def$gene1[[gsub("^Assay: ", "", input$sc1a1ass1)]]
      resLen <- 7
    }
    return(list(res, resDef, resLen))
  })
  observeEvent(gsub("^Assay: ", "", input$sc1a1ass1), {
    updateSelectizeInput(
      session,
      "sc1a1inp1",
      choices = getGsc1a1inp1()[[1]],
      server = TRUE,
      selected = getGsc1a1inp1()[[2]],
      options = list(
        maxOptions = getGsc1a1inp1()[[3]],
        create = TRUE,
        persist = TRUE,
        render = I(optCrt)
      )
    )
  })
  output$sc1a1sub1.ui <- renderUI({
    sub = strsplit(sc1conf[UI == input$sc1a1sub1]$fID, "\\|")[[1]]
    checkboxGroupInput(
      "sc1a1sub2",
      "Select which cells to show",
      inline = TRUE,
      choices = sub,
      selected = sub
    )
  })
  observeEvent(input$sc1a1sub1non, {
    sub = strsplit(sc1conf[UI == input$sc1a1sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(
      session,
      inputId = "sc1a1sub2",
      label = "Select which cells to show",
      choices = sub,
      selected = NULL,
      inline = TRUE
    )
  })
  observeEvent(input$sc1a1sub1all, {
    sub = strsplit(sc1conf[UI == input$sc1a1sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(
      session,
      inputId = "sc1a1sub2",
      label = "Select which cells to show",
      choices = sub,
      selected = sub,
      inline = TRUE
    )
  })

  sc1a1oup1xy <- reactiveValues(x = NULL, y = NULL)
  observe({
    brush <- input$sc1a1inp1.br
    if (!is.null(brush)) {
      sc1a1oup1xy$x <- c(brush$xmin, brush$xmax)
      sc1a1oup1xy$y <- c(brush$ymin, brush$ymax)
    } else {
      sc1a1oup1xy$x <- NULL
      sc1a1oup1xy$y <- NULL
    }
  })
  sc1a1oup1br <- reactive({
    sc2Ddimr(
      sc1conf,
      sc1meta,
      sc1dimr,
      input$sc1a1dr,
      input$sc1a1inp1,
      "sc1assay_",
      sc1gene,
      input$sc1a1ass1,
      input$sc1a1sub1,
      input$sc1a1sub2,
      input$sc1a1min1,
      input$sc1a1max1,
      input$sc1a1siz / 2,
      input$sc1a1ord1,
      cList[[input$sc1a1col1]],
      sList[input$sc1a1fsz] / 2,
      input$sc1a1asp,
      FALSE,
      FALSE
    )
  })
  output$sc1a1oup1.br <- renderPlot({
    sc1a1oup1br() + theme(legend.position = "none")
  })

  sc1a1oup1 <- reactive({
    sc2Ddimr(
      sc1conf,
      sc1meta,
      sc1dimr,
      input$sc1a1dr,
      input$sc1a1inp1,
      "sc1assay_",
      sc1gene,
      input$sc1a1ass1,
      input$sc1a1sub1,
      input$sc1a1sub2,
      input$sc1a1min1,
      input$sc1a1max1,
      input$sc1a1siz,
      input$sc1a1ord1,
      cList[[input$sc1a1col1]],
      sList[input$sc1a1fsz],
      input$sc1a1asp,
      input$sc1a1txt,
      input$sc1a1lab1
    )
  })
  output$sc1a1oup1 <- renderPlot({
    if (is.null(sc1a1oup1xy$x[1])) {
      sc1a1oup1() + theme(legend.position = "none")
    } else {
      sc1a1oup1() +
        theme(legend.position = "none") +
        scale_x_continuous(limits = sc1a1oup1xy$x, expand = c(0, 0)) +
        scale_y_continuous(limits = sc1a1oup1xy$y, expand = c(0, 0))
    }
  })
  output$sc1a1oup1.ui <- renderUI({
    plotOutput("sc1a1oup1", height = pList[input$sc1a1psz])
  })
  output$sc1a1oup1.dl <- downloadHandler(
    filename = function() {
      paste0("sc1", input$sc1a1dr, "_", input$sc1a1inp1, ".", input$sc1a1oup1.f)
    },
    content = function(file) {
      if (is.null(sc1a1oup1xy$x[1])) {
        ggsav(
          file,
          height = input$sc1a1oup1.h,
          width = input$sc1a1oup1.w,
          plot = sc1a1oup1() + theme(legend.position = "none")
        )
      } else {
        ggsav(
          file,
          height = input$sc1a1oup1.h,
          width = input$sc1a1oup1.w,
          plot = sc1a1oup1() +
            theme(legend.position = "none") +
            scale_x_continuous(limits = sc1a1oup1xy$x, expand = c(0, 0)) +
            scale_y_continuous(limits = sc1a1oup1xy$y, expand = c(0, 0))
        )
      }
    }
  )
  output$sc1a1oup3 <- renderPlot({
    grid.newpage()
    grid.draw(g_legend(sc1a1oup1()))
  })
  output$sc1a1oup3.ui <- renderUI({
    plotOutput(
      "sc1a1oup3",
      height = 72 *
        convertHeight(
          grobHeight(g_legend(sc1a1oup1())),
          unitTo = "in",
          valueOnly = TRUE
        ) +
        50
    )
  })
  output$sc1a1oup3.dl <- downloadHandler(
    filename = function() {
      paste0(
        "sc1",
        input$sc1a1dr,
        "_",
        input$sc1a1inp1,
        "_leg.",
        input$sc1a1oup3.f
      )
    },
    content = function(file) {
      grid.newpage()
      grid.draw(g_legend(sc1a1oup1()))
      ggsav(
        file,
        height = input$sc1a1oup3.h,
        width = input$sc1a1oup3.w,
        plot = grid.grab()
      )
    }
  )

  output$sc1a1.dt <- renderDataTable({
    ggData = sc2Dnum(
      sc1conf,
      sc1meta,
      sc1dimr,
      input$sc1a1dr,
      sc1a1oup1xy$x,
      sc1a1oup1xy$y,
      input$sc1a1inp1,
      "sc1assay_",
      sc1gene,
      input$sc1a1ass1,
      input$sc1a1sub1,
      input$sc1a1sub2,
      input$sc1a1splt
    )
    datatable(
      ggData,
      rownames = FALSE,
      extensions = "Buttons",
      options = list(
        pageLength = -1,
        dom = "tB",
        buttons = c("copy", "csv", "excel")
      )
    ) %>%
      formatRound(columns = c("pctZoom"), digits = 2)
  })

  ### Functions for tab A2
  getGsc1a2inp1 <- reactive({
    req(gsub("^Assay: ", "", input$sc1a2ass1))
    if (gsub("^Assay: ", "", input$sc1a2ass1) == "Cell Information") {
      res <- sc1conf$UI
      resDef <- sc1def$meta1
      resLen <- length(res)
    } else {
      res <- names(sc1gene[[gsub("^Assay: ", "", input$sc1a2ass1)]])
      resDef <- sc1def$gene1[[gsub("^Assay: ", "", input$sc1a2ass1)]]
      resLen <- 7
    }
    return(list(res, resDef, resLen))
  })
  observeEvent(gsub("^Assay: ", "", input$sc1a2ass1), {
    updateSelectizeInput(
      session,
      "sc1a2inp1",
      choices = getGsc1a2inp1()[[1]],
      server = TRUE,
      selected = getGsc1a2inp1()[[2]],
      options = list(
        maxOptions = getGsc1a2inp1()[[3]],
        create = TRUE,
        persist = TRUE,
        render = I(optCrt)
      )
    )
  })
  getGsc1a2inp2 <- reactive({
    req(gsub("^Assay: ", "", input$sc1a2ass2))
    if (gsub("^Assay: ", "", input$sc1a2ass2) == "Cell Information") {
      res <- sc1conf$UI
      resDef <- sc1def$meta1
      resLen <- length(res)
    } else {
      res <- names(sc1gene[[gsub("^Assay: ", "", input$sc1a2ass2)]])
      resDef <- sc1def$gene1[[gsub("^Assay: ", "", input$sc1a2ass2)]]
      resLen <- 7
    }
    return(list(res, resDef, resLen))
  })
  observeEvent(gsub("^Assay: ", "", input$sc1a2ass2), {
    updateSelectizeInput(
      session,
      "sc1a2inp2",
      choices = getGsc1a2inp2()[[1]],
      server = TRUE,
      selected = getGsc1a2inp2()[[2]],
      options = list(
        maxOptions = getGsc1a2inp2()[[3]],
        create = TRUE,
        persist = TRUE,
        render = I(optCrt)
      )
    )
  })
  output$sc1a2sub1.ui <- renderUI({
    sub = strsplit(sc1conf[UI == input$sc1a2sub1]$fID, "\\|")[[1]]
    checkboxGroupInput(
      "sc1a2sub2",
      "Select which cells to show",
      inline = TRUE,
      choices = sub,
      selected = sub
    )
  })
  observeEvent(input$sc1a2sub1non, {
    sub = strsplit(sc1conf[UI == input$sc1a2sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(
      session,
      inputId = "sc1a2sub2",
      label = "Select which cells to show",
      choices = sub,
      selected = NULL,
      inline = TRUE
    )
  })
  observeEvent(input$sc1a2sub1all, {
    sub = strsplit(sc1conf[UI == input$sc1a2sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(
      session,
      inputId = "sc1a2sub2",
      label = "Select which cells to show",
      choices = sub,
      selected = sub,
      inline = TRUE
    )
  })

  sc1a2oup1 <- reactive({
    sc2Ddimr(
      sc1conf,
      sc1meta,
      sc1dimr,
      input$sc1a2dr,
      input$sc1a2inp1,
      "sc1assay_",
      sc1gene,
      input$sc1a2ass1,
      input$sc1a2sub1,
      input$sc1a2sub2,
      input$sc1a2min1,
      input$sc1a2max1,
      input$sc1a2siz,
      input$sc1a2ord1,
      cList[[input$sc1a2col1]],
      sList[input$sc1a2fsz],
      input$sc1a2asp,
      input$sc1a2txt,
      input$sc1a2lab1
    )
  })
  output$sc1a2oup1 <- renderPlot({
    sc1a2oup1() + theme(legend.position = "none")
  })
  output$sc1a2oup1.ui <- renderUI({
    plotOutput("sc1a2oup1", height = pList[input$sc1a2psz])
  })
  output$sc1a2oup1.dl <- downloadHandler(
    filename = function() {
      paste0("sc1", input$sc1a2dr, "_", input$sc1a2inp1, ".", input$sc1a2oup1.f)
    },
    content = function(file) {
      ggsav(
        file,
        height = input$sc1a2oup1.h,
        width = input$sc1a2oup1.w,
        plot = sc1a2oup1() + theme(legend.position = "none")
      )
    }
  )
  output$sc1a2oup3 <- renderPlot({
    grid.newpage()
    grid.draw(g_legend(sc1a2oup1()))
  })
  output$sc1a2oup3.ui <- renderUI({
    plotOutput(
      "sc1a2oup3",
      height = 72 *
        convertHeight(
          grobHeight(g_legend(sc1a2oup1())),
          unitTo = "in",
          valueOnly = TRUE
        ) +
        50
    )
  })
  output$sc1a2oup3.dl <- downloadHandler(
    filename = function() {
      paste0(
        "sc1",
        input$sc1a2dr,
        "_",
        input$sc1a2inp1,
        "_leg.",
        input$sc1a2oup3.f
      )
    },
    content = function(file) {
      grid.newpage()
      grid.draw(g_legend(sc1a2oup1()))
      ggsav(
        file,
        height = input$sc1a2oup3.h,
        width = input$sc1a2oup3.w,
        plot = grid.grab()
      )
    }
  )

  sc1a2oup2 <- reactive({
    sc2Ddimr(
      sc1conf,
      sc1meta,
      sc1dimr,
      input$sc1a2dr,
      input$sc1a2inp2,
      "sc1assay_",
      sc1gene,
      input$sc1a2ass2,
      input$sc1a2sub1,
      input$sc1a2sub2,
      input$sc1a2min2,
      input$sc1a2max2,
      input$sc1a2siz,
      input$sc1a2ord2,
      cList[[input$sc1a2col2]],
      sList[input$sc1a2fsz],
      input$sc1a2asp,
      input$sc1a2txt,
      input$sc1a2lab2
    )
  })
  output$sc1a2oup2 <- renderPlot({
    sc1a2oup2() + theme(legend.position = "none")
  })
  output$sc1a2oup2.ui <- renderUI({
    plotOutput("sc1a2oup2", height = pList[input$sc1a2psz])
  })
  output$sc1a2oup2.dl <- downloadHandler(
    filename = function() {
      paste0("sc1", input$sc1a2dr, "_", input$sc1a2inp2, ".", input$sc1a2oup2.f)
    },
    content = function(file) {
      ggsav(
        file,
        height = input$sc1a2oup2.h,
        width = input$sc1a2oup2.w,
        plot = sc1a2oup2() + theme(legend.position = "none")
      )
    }
  )

  output$sc1a2oup4 <- renderPlot({
    grid.newpage()
    grid.draw(g_legend(sc1a2oup2()))
  })
  output$sc1a2oup4.ui <- renderUI({
    plotOutput(
      "sc1a2oup4",
      height = 72 *
        convertHeight(
          grobHeight(g_legend(sc1a2oup2())),
          unitTo = "in",
          valueOnly = TRUE
        ) +
        50
    )
  })
  output$sc1a2oup4.dl <- downloadHandler(
    filename = function() {
      paste0(
        "sc1",
        input$sc1a2dr,
        "_",
        input$sc1a2inp2,
        "_leg.",
        input$sc1a2oup4.f
      )
    },
    content = function(file) {
      grid.newpage()
      grid.draw(g_legend(sc1a2oup2()))
      ggsav(
        file,
        height = input$sc1a2oup4.h,
        width = input$sc1a2oup4.w,
        plot = grid.grab()
      )
    }
  )

  sc1a2oup5 <- reactive({
    sc2Dcomp(
      sc1conf,
      sc1meta,
      input$sc1a2inp1,
      input$sc1a2inp2,
      "sc1assay_",
      sc1gene,
      input$sc1a2ass1,
      input$sc1a2ass2,
      input$sc1a2sub1,
      input$sc1a2sub2,
      input$sc1a2min1,
      input$sc1a2max1,
      input$sc1a2min2,
      input$sc1a2max2,
      input$sc1a2siz,
      cList[[input$sc1a2col2]],
      sList[input$sc1a2fsz]
    )
  })
  output$sc1a2oup5 <- renderPlot({
    sc1a2oup5() + theme(legend.position = "none")
  })
  output$sc1a2oup5.ui <- renderUI({
    plotOutput("sc1a2oup5", height = pList[input$sc1a2psz])
  })
  output$sc1a2oup5.dl <- downloadHandler(
    filename = function() {
      paste0(
        "sc1compare_",
        input$sc1a2inp1,
        "_",
        input$sc1a2inp2,
        ".",
        input$sc1a2oup5.f
      )
    },
    content = function(file) {
      ggsav(
        file,
        height = input$sc1a2oup5.h,
        width = input$sc1a2oup5.w,
        plot = sc1a2oup5()
      )
    }
  )
  output$sc1a2.dt <- renderDataTable({
    ggData = sc2Dcnum(
      sc1conf,
      sc1meta,
      input$sc1a2inp1,
      input$sc1a2inp2,
      "sc1assay_",
      sc1gene,
      input$sc1a2ass1,
      input$sc1a2ass2,
      input$sc1a2sub1,
      input$sc1a2sub2,
      input$sc1a2min1,
      input$sc1a2max1,
      input$sc1a2min2,
      input$sc1a2max2,
      input$sc1a2cut
    )
    datatable(
      ggData,
      rownames = FALSE,
      extensions = "Buttons",
      options = list(
        pageLength = -1,
        dom = "tB",
        buttons = c("copy", "csv", "excel")
      )
    )
  })

  ### Functions for tab A3
  getGsc1a3inp1 <- reactive({
    req(input$sc1a3ass1)
    res <- names(sc1gene[[input$sc1a3ass1]])
    resDef1 <- sc1def$gene1[[input$sc1a3ass1]]
    resDef2 <- sc1def$gene2[[input$sc1a3ass1]]
    return(list(res, resDef1, resDef2))
  })
  observeEvent(input$sc1a3ass1, {
    updateSelectizeInput(
      session,
      "sc1a3inp1",
      choices = getGsc1a3inp1()[[1]],
      server = TRUE,
      selected = getGsc1a3inp1()[[2]],
      options = list(
        maxOptions = 7,
        create = TRUE,
        persist = TRUE,
        render = I(optCrt)
      )
    )
  })
  observeEvent(input$sc1a3ass1, {
    updateSelectizeInput(
      session,
      "sc1a3inp2",
      choices = getGsc1a3inp1()[[1]],
      server = TRUE,
      selected = getGsc1a3inp1()[[3]],
      options = list(
        maxOptions = 7,
        create = TRUE,
        persist = TRUE,
        render = I(optCrt)
      )
    )
  })
  output$sc1a3sub1.ui <- renderUI({
    sub = strsplit(sc1conf[UI == input$sc1a3sub1]$fID, "\\|")[[1]]
    checkboxGroupInput(
      "sc1a3sub2",
      "Select which cells to show",
      inline = TRUE,
      choices = sub,
      selected = sub
    )
  })
  observeEvent(input$sc1a3sub1non, {
    sub = strsplit(sc1conf[UI == input$sc1a3sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(
      session,
      inputId = "sc1a3sub2",
      label = "Select which cells to show",
      choices = sub,
      selected = NULL,
      inline = TRUE
    )
  })
  observeEvent(input$sc1a3sub1all, {
    sub = strsplit(sc1conf[UI == input$sc1a3sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(
      session,
      inputId = "sc1a3sub2",
      label = "Select which cells to show",
      choices = sub,
      selected = sub,
      inline = TRUE
    )
  })
  sc1a3oup1 <- reactive({
    scDRcoex(
      sc1conf,
      sc1meta,
      sc1dimr,
      input$sc1a3dr,
      input$sc1a3inp1,
      input$sc1a3inp2,
      paste0("sc1assay_", input$sc1a3ass1, ".h5"),
      sc1gene[[input$sc1a3ass1]],
      input$sc1a3sub1,
      input$sc1a3sub2,
      input$sc1a3min1,
      input$sc1a3max1,
      input$sc1a3min2,
      input$sc1a3max2,
      input$sc1a3siz,
      input$sc1a3col1,
      input$sc1a3ord1,
      sList[input$sc1a3fsz],
      input$sc1a3asp,
      input$sc1a3txt
    )
  })
  output$sc1a3oup1 <- renderPlot({
    sc1a3oup1()
  })
  output$sc1a3oup1.ui <- renderUI({
    plotOutput("sc1a3oup1", height = pList2[input$sc1a3psz])
  })
  output$sc1a3oup1.dl <- downloadHandler(
    filename = function() {
      paste0(
        "sc1",
        input$sc1a3dr,
        "_",
        input$sc1a3inp1,
        "_",
        input$sc1a3inp2,
        ".",
        input$sc1a3oup1.f
      )
    },
    content = function(file) {
      ggsav(
        file,
        height = input$sc1a3oup1.h,
        width = input$sc1a3oup1.w,
        plot = sc1a3oup1()
      )
    }
  )

  output$sc1a3oup2 <- renderPlot({
    scDRcoexLeg(
      input$sc1a3inp1,
      input$sc1a3inp2,
      input$sc1a3col1,
      sList[input$sc1a3fsz]
    )
  })
  output$sc1a3oup2.ui <- renderUI({
    plotOutput("sc1a3oup2", height = "300px")
  })
  output$sc1a3oup2.dl <- downloadHandler(
    filename = function() {
      paste0(
        "sc1",
        input$sc1a3dr,
        "_",
        input$sc1a3inp1,
        "_",
        input$sc1a3inp2,
        "_leg.",
        input$sc1a3oup2.f
      )
    },
    content = function(file) {
      ggsav(
        file,
        height = 3,
        width = 4,
        plot = scDRcoexLeg(
          input$sc1a3inp1,
          input$sc1a3inp2,
          input$sc1a3col1,
          sList[input$sc1a3fsz]
        )
      )
    }
  )
  output$sc1a3.dt <- renderDataTable({
    ggData = scDRcoexNum(
      sc1conf,
      sc1meta,
      input$sc1a3inp1,
      input$sc1a3inp2,
      paste0("sc1assay_", input$sc1a3ass1, ".h5"),
      sc1gene[[input$sc1a3ass1]],
      input$sc1a3sub1,
      input$sc1a3sub2
    )
    datatable(
      ggData,
      rownames = FALSE,
      extensions = "Buttons",
      options = list(
        pageLength = -1,
        dom = "tB",
        buttons = c("copy", "csv", "excel")
      )
    ) %>%
      formatRound(columns = c("percent"), digits = 2)
  })

  ### Functions for tab B1
  getGsc1b1inp1 <- reactive({
    req(gsub("^Assay: ", "", input$sc1b1ass1))
    if (gsub("^Assay: ", "", input$sc1b1ass1) == "Cell Information") {
      res <- sc1conf[is.na(fID)]$UI
      resDef <- sc1conf[is.na(fID)]$UI[1]
      resLen <- length(res)
    } else {
      res <- names(sc1gene[[gsub("^Assay: ", "", input$sc1b1ass1)]])
      resDef <- sc1def$gene1[[gsub("^Assay: ", "", input$sc1b1ass1)]]
      resLen <- 7
    }
    return(list(res, resDef, resLen))
  })
  observeEvent(gsub("^Assay: ", "", input$sc1b1ass1), {
    updateSelectizeInput(
      session,
      "sc1b1inp2",
      choices = getGsc1b1inp1()[[1]],
      server = TRUE,
      selected = getGsc1b1inp1()[[2]],
      options = list(
        maxOptions = getGsc1b1inp1()[[3]],
        create = TRUE,
        persist = TRUE,
        render = I(optCrt)
      )
    )
  })
  output$sc1b1sub1.ui <- renderUI({
    sub = strsplit(sc1conf[UI == input$sc1b1sub1]$fID, "\\|")[[1]]
    checkboxGroupInput(
      "sc1b1sub2",
      "Select which cells to show",
      inline = TRUE,
      choices = sub,
      selected = sub
    )
  })
  observeEvent(input$sc1b1sub1non, {
    sub = strsplit(sc1conf[UI == input$sc1b1sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(
      session,
      inputId = "sc1b1sub2",
      label = "Select which cells to show",
      choices = sub,
      selected = NULL,
      inline = TRUE
    )
  })
  observeEvent(input$sc1b1sub1all, {
    sub = strsplit(sc1conf[UI == input$sc1b1sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(
      session,
      inputId = "sc1b1sub2",
      label = "Select which cells to show",
      choices = sub,
      selected = sub,
      inline = TRUE
    )
  })

  sc1b1oup <- reactive({
    scVioBox(
      sc1conf,
      sc1meta,
      input$sc1b1inp1,
      input$sc1b1inp2,
      "sc1assay_",
      sc1gene,
      input$sc1b1ass1,
      input$sc1b1sub1,
      input$sc1b1sub2,
      input$sc1b1typ,
      input$sc1b1pts,
      input$sc1b1stg,
      input$sc1b1stp1,
      input$sc1b1stp2,
      input$sc1b1siz,
      sList[input$sc1b1fsz],
      input$sc1b1noi
    )
  })
  output$sc1b1oup <- renderPlot({
    sc1b1oup()
  })
  output$sc1b1oup.ui <- renderUI({
    plotOutput("sc1b1oup", height = pList2[input$sc1b1psz])
  })
  output$sc1b1oup.dl <- downloadHandler(
    filename = function() {
      paste0(
        "sc1",
        input$sc1b1typ,
        "_",
        input$sc1b1inp1,
        "_",
        input$sc1b1inp2,
        ".",
        input$sc1b1oup.f
      )
    },
    content = function(file) {
      ggsave(
        file,
        height = input$sc1b1oup.h,
        width = input$sc1b1oup.w,
        plot = sc1b1oup()
      )
    }
  )

  ### Functions for tab B2
  observeEvent(input$sc1b2inp1, {
    sub = strsplit(sc1conf[UI == input$sc1b2inp2]$fID, "\\|")[[1]]
    updateSelectizeInput(
      session,
      "sc1b2ord1",
      choices = c("Original order", sub),
      server = TRUE,
      selected = "Original order",
      options = list(
        create = TRUE,
        persist = TRUE,
        render = I(optCrt)
      )
    )
  })
  output$sc1b2sub1.ui <- renderUI({
    sub = strsplit(sc1conf[UI == input$sc1b2sub1]$fID, "\\|")[[1]]
    checkboxGroupInput(
      "sc1b2sub2",
      "Select which cells to show",
      inline = TRUE,
      choices = sub,
      selected = sub
    )
  })
  observeEvent(input$sc1b2sub1non, {
    sub = strsplit(sc1conf[UI == input$sc1b2sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(
      session,
      inputId = "sc1b2sub2",
      label = "Select which cells to show",
      choices = sub,
      selected = NULL,
      inline = TRUE
    )
  })
  observeEvent(input$sc1b2sub1all, {
    sub = strsplit(sc1conf[UI == input$sc1b2sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(
      session,
      inputId = "sc1b2sub2",
      label = "Select which cells to show",
      choices = sub,
      selected = sub,
      inline = TRUE
    )
  })

  sc1b2oup <- reactive({
    scProp(
      sc1conf,
      sc1meta,
      input$sc1b2inp1,
      input$sc1b2inp2,
      input$sc1b2sub1,
      input$sc1b2sub2,
      input$sc1b2ord1,
      input$sc1b2ord2,
      input$sc1b2typ,
      input$sc1b2flp,
      sList[input$sc1b2fsz]
    )
  })
  output$sc1b2oup <- renderPlot({
    sc1b2oup()
  })
  output$sc1b2oup.ui <- renderUI({
    plotOutput("sc1b2oup", height = pList2[input$sc1b2psz])
  })
  output$sc1b2oup.dl <- downloadHandler(
    filename = function() {
      paste0(
        "sc1",
        input$sc1b2typ,
        "_",
        input$sc1b2inp1,
        "_",
        input$sc1b2inp2,
        ".",
        input$sc1b2oup.f
      )
    },
    content = function(file) {
      ggsave(
        file,
        height = input$sc1b2oup.h,
        width = input$sc1b2oup.w,
        plot = sc1b2oup()
      )
    }
  )

  ### Functions for tab B3
  getGsc1b3inp1 <- reactive({
    req(input$sc1b3ass1)
    resDef <- paste0(sc1def$genes[[input$sc1b3ass1]], collapse = ", ")
    return(resDef)
  })
  observeEvent(input$sc1b3ass1, {
    updateTextAreaInput(session, "sc1b3inp", value = getGsc1b3inp1())
  })
  output$sc1b3sub1.ui <- renderUI({
    sub = strsplit(sc1conf[UI == input$sc1b3sub1]$fID, "\\|")[[1]]
    checkboxGroupInput(
      "sc1b3sub2",
      "Select which cells to show",
      inline = TRUE,
      choices = sub,
      selected = sub
    )
  })
  observeEvent(input$sc1b3sub1non, {
    sub = strsplit(sc1conf[UI == input$sc1b3sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(
      session,
      inputId = "sc1b3sub2",
      label = "Select which cells to show",
      choices = sub,
      selected = NULL,
      inline = TRUE
    )
  })
  observeEvent(input$sc1b3sub1all, {
    sub = strsplit(sc1conf[UI == input$sc1b3sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(
      session,
      inputId = "sc1b3sub2",
      label = "Select which cells to show",
      choices = sub,
      selected = sub,
      inline = TRUE
    )
  })
  output$sc1b3oupTxt <- renderUI({
    geneList = scGeneList(input$sc1b3inp, sc1gene[[input$sc1b3ass1]])
    if (nrow(geneList) > 50) {
      HTML("More than 50 input genes! Please reduce the gene list!")
    } else {
      oup = paste0(
        nrow(geneList[present == TRUE]),
        " genes OK and will be plotted"
      )
      if (nrow(geneList[present == FALSE]) > 0) {
        oup = paste0(
          oup,
          "<br/>",
          nrow(geneList[present == FALSE]),
          " genes not found (",
          paste0(geneList[present == FALSE]$gene, collapse = ", "),
          ")"
        )
      }
      HTML(oup)
    }
  })

  sc1b3oup <- reactive({
    scBubbHeat(
      sc1conf,
      sc1meta,
      input$sc1b3inp,
      input$sc1b3grp,
      input$sc1b3plt,
      paste0("sc1assay_", input$sc1b3ass1, ".h5"),
      sc1gene[[input$sc1b3ass1]],
      input$sc1b3sub1,
      input$sc1b3sub2,
      input$sc1b3scl,
      input$sc1b3row,
      input$sc1b3col,
      cList[[input$sc1b3cols]],
      sList[input$sc1b3fsz],
      input$sc1b3exp,
      input$sc1b3max
    )
  })
  output$sc1b3oup <- renderPlot({
    sc1b3oup()
  })
  output$sc1b3oup.ui <- renderUI({
    plotOutput("sc1b3oup", height = pList3[input$sc1b3psz])
  })
  output$sc1b3oup.dl <- downloadHandler(
    filename = function() {
      paste0("sc1", input$sc1b3plt, "_", input$sc1b3grp, ".", input$sc1b3oup.f)
    },
    content = function(file) {
      ggsave(
        file,
        height = input$sc1b3oup.h,
        width = input$sc1b3oup.w,
        plot = sc1b3oup()
      )
    }
  )
  output$deg_plot_lvl1 <- renderPlot({
    req(input$gene_input)
    gene <- input$gene_input

    df_lvl1 <- subset(
      mast_de,
      hgnc_symbol == gene & celltype_level == "level1"
    )
    if (nrow(df_lvl1) == 0) {
      return(NULL)
    }

    ggplot(
      df_lvl1,
      aes(
        x = reorder(celltypes_abbrev, avg_log2FC),
        y = avg_log2FC,
        fill = -log10(p_val_adj)
      )
    ) +
      geom_col() +
      coord_flip() +
      labs(
        title = paste("Level 1 Celltypes: DE for", gene),
        x = "",
        y = "Avg.logFC",
        fill = "-log10(adj.P.Val)"
      ) +
      theme_minimal()
  })

  output$deg_plot_lvl2 <- renderPlot({
    req(input$gene_input)
    gene <- input$gene_input

    df_lvl2 <- subset(
      mast_de,
      hgnc_symbol == gene & celltype_level == "level2"
    )
    if (nrow(df_lvl2) == 0) {
      return(NULL)
    }

    ggplot(
      df_lvl2,
      aes(
        x = reorder(celltypes_abbrev, avg_log2FC),
        y = avg_log2FC,
        fill = -log10(p_val_adj)
      )
    ) +
      geom_col() +
      coord_flip() +
      labs(
        title = paste("Level 2 Celltypes: DE for", gene),
        x = "",
        y = "Avg.logFC",
        fill = "-log10(adj.P.Val)"
      ) +
      theme_minimal()
  })
  output$tradeseq_pseudotime_plot <- renderPlot({
    req(input$pseudotime_gene_input)
    gene <- input$pseudotime_gene_input
    pseudotime_genes <- subset(all_pseudotime_sig_genes, hgnc_symbol == gene)

    plot_obj <- plotSmoothers(
      sce_subset,
      counts,
      pseudotime_genes$ensembl_gene_id,
      curvesCols = pal_npg("nrc")(4) #,
      #size = 0.2,
      #lwd = 0.8
    ) +
      scale_colour_npg(
        labels = c(
          "EndoMT-SMC - AD",
          "EndoMT-SMC - Control",
          "EndoMT-PC - AD",
          "EndoMT-PC - Control"
        )
      ) +
      ggtitle(pseudotime_genes$hgnc_symbol) +
      labs(colour = "Lineage") +
      theme(
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)
      )

    if (!is.null(plot_obj)) {
      plot_obj
    } else {
      plot.new()
      text(0.5, 0.5, "No plot available for selected gene")
    }
  })
})
