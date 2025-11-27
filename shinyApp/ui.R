library(shiny)
library(shinyhelper)
library(data.table)
library(Matrix)
library(DT)
library(magrittr)
sc1conf = readRDS("./sc1conf.rds")
sc1def = readRDS("./sc1def.rds")

### Start server code
shinyUI(fluidPage(
  ### HTML formatting of error messages

  tags$head(tags$style(HTML(
    ".shiny-output-error-validation {color: red; font-weight: bold;}"
  ))),
  list(tags$style(HTML(
    ".navbar-default .navbar-nav { font-weight: bold; font-size: 16px; }"
  ))),

  ### Page title
  titlePanel("AD Vascular-Parenchymal Paired-seq"),
  navbarPage(
    NULL,
    ### Tab1.a1: Zoom-enable Dimred
    tabPanel(
      HTML("Zoom-enable Dimred"),
      h4(
        "Zoom-enable reduced dimensions overlaid with cell info or assay expression"
      ),
      "In this tab, users can visualise either cell information and gene expression ",
      "in a zoom-enabled low-dimensional representions plots.",
      br(),
      br(),
      fluidRow(
        column(
          3,
          fluidRow(
            column(
              12,
              selectInput(
                "sc1a1dr",
                "Reduction:",
                choices = sc1def$dimrd,
                selected = sc1def$dimrd[1]
              )
            )
          )
        ), # End of column (6 space)
        column(
          3,
          actionButton("sc1a1togL", "Toggle to subset cells"),
          conditionalPanel(
            condition = "input.sc1a1togL % 2 == 1",
            selectInput(
              "sc1a1sub1",
              "Cell information to subset:",
              choices = sc1conf[grp == TRUE]$UI,
              selected = sc1def$grp1
            ),
            uiOutput("sc1a1sub1.ui"),
            actionButton(
              "sc1a1sub1all",
              "Select all groups",
              class = "btn btn-primary"
            ),
            actionButton(
              "sc1a1sub1non",
              "Deselect all groups",
              class = "btn btn-primary"
            )
          )
        ), # End of column (6 space)
        column(
          6,
          actionButton("sc1a1tog0", "Toggle graphics controls"),
          conditionalPanel(
            condition = "input.sc1a1tog0 % 2 == 1",
            fluidRow(
              column(
                6,
                sliderInput(
                  "sc1a1siz",
                  "Point size:",
                  min = 0,
                  max = 4,
                  value = 1.25,
                  step = 0.25
                ),
                radioButtons(
                  "sc1a1psz",
                  "Plot size:",
                  choices = c("Small", "Medium", "Large"),
                  selected = "Medium",
                  inline = TRUE
                ),
                radioButtons(
                  "sc1a1fsz",
                  "Font size:",
                  choices = c("Small", "Medium", "Large"),
                  selected = "Medium",
                  inline = TRUE
                )
              ),
              column(
                6,
                radioButtons(
                  "sc1a1asp",
                  "Aspect ratio:",
                  choices = c("Square", "Fixed", "Free"),
                  selected = "Square",
                  inline = TRUE
                ),
                checkboxInput("sc1a1txt", "Show axis text", value = FALSE)
              )
            )
          )
        ) # End of column (6 space)
      ), # End of fluidRow (4 space)
      fluidRow(
        column(
          3,
          style = "border-right: 2px solid black",
          h4("Information to plot"),
          selectInput(
            "sc1a1ass1",
            "Data type to colour plot:",
            choices = c("Cell Information", paste0("Assay: ", sc1def$assay)),
            selected = "Cell Information"
          ),
          selectInput(
            "sc1a1inp1",
            "Cell Info / Feature Name:",
            choices = NULL
          ) %>%
            helper(
              type = "inline",
              size = "m",
              fade = TRUE,
              title = "Cell Info / Gene to colour cells by",
              content = c(
                "Select cell info / feature to colour cells",
                "- Categorical covariates have a fixed colour palette",
                paste0(
                  "- Continuous covariates / gene expression are coloured ",
                  "in a Blue-Yellow-Red colour scheme, which can be ",
                  "changed in the plot controls"
                )
              )
            ),
          strong("Draw box on plot below to zoom: "),
          plotOutput(
            "sc1a1oup1.br",
            height = "400px",
            brush = brushOpts(id = "sc1a1inp1.br", resetOnNew = TRUE)
          ),
          actionButton("sc1a1tog1", "Toggle plot controls"),
          conditionalPanel(
            condition = "input.sc1a1tog1 % 2 == 1",
            radioButtons(
              "sc1a1col1",
              "Colour (Continuous data):",
              choices = c(
                "White-Red",
                "Blue-Yellow-Red",
                "Yellow-Green-Purple"
              ),
              selected = "Blue-Yellow-Red"
            ),
            numericInput(
              "sc1a1min1",
              "Min cutoff (q##):",
              min = 0,
              max = 50,
              value = 0,
              step = 1
            ),
            numericInput(
              "sc1a1max1",
              "Max cutoff (q##):",
              min = 50,
              max = 100,
              value = 100,
              step = 1
            ),
            radioButtons(
              "sc1a1ord1",
              "Plot order:",
              choices = c("Max-1st", "Min-1st", "Original", "Random"),
              selected = "Original",
              inline = TRUE
            ),
            checkboxInput("sc1a1lab1", "Show cell info labels", value = TRUE)
          )
        ), # End of column (6 space)
        column(
          6,
          style = "border-right: 2px solid black",
          fluidRow(column(12, uiOutput("sc1a1oup1.ui"))),
          fluidRow(column(12, uiOutput("sc1a1oup3.ui"))),
          fluidRow(
            column(3, downloadButton("sc1a1oup1.dl", "Download Plot")),
            column(
              3,
              radioButtons(
                "sc1a1oup1.f",
                "Plot format:",
                choices = c("png", "pdf"),
                selected = "png",
                inline = TRUE
              )
            ),
            column(
              3,
              numericInput(
                "sc1a1oup1.h",
                "Plot height:",
                min = 4,
                max = 20,
                value = 8,
                step = 0.5
              )
            ),
            column(
              3,
              numericInput(
                "sc1a1oup1.w",
                "Plot width:",
                min = 4,
                max = 20,
                value = 8,
                step = 0.5
              )
            )
          ),
          fluidRow(
            column(3, downloadButton("sc1a1oup3.dl", "Download Legend")),
            column(
              3,
              radioButtons(
                "sc1a1oup3.f",
                "Legend format:",
                choices = c("png", "pdf"),
                selected = "png",
                inline = TRUE
              )
            ),
            column(
              3,
              numericInput(
                "sc1a1oup3.h",
                "Legend height:",
                min = 0.2,
                max = 10,
                value = 2,
                step = 0.1
              )
            ),
            column(
              3,
              numericInput(
                "sc1a1oup3.w",
                "Legend width:",
                min = 4,
                max = 20,
                value = 8,
                step = 0.5
              )
            )
          ),
          br()
        ), # End of column (6 space)
        column(
          3,
          h4("Cell numbers"),
          numericInput(
            "sc1a1splt",
            "Split continuous cell info into nBins:",
            min = 2,
            max = 10,
            value = 4,
            step = 1
          ),
          dataTableOutput("sc1a1.dt")
        ) # End of column (6 space)
      ) # End of fluidRow (4 space)
    ), # End of tab (2 space)
    ### Tab1.a2: CellInfo vs AssayExpr on dimRed
    tabPanel(
      HTML("Side-by-side DimRed"),
      h4("Cell information vs assay expression on reduced dimensions"),
      "In this tab, users can visualise both cell information and gene ",
      "expression side-by-side on low-dimensional representions.",
      br(),
      br(),
      fluidRow(
        column(
          3,
          fluidRow(
            column(
              12,
              selectInput(
                "sc1a2dr",
                "Reduction:",
                choices = sc1def$dimrd,
                selected = sc1def$dimrd[1]
              )
            )
          )
        ), # End of column (6 space)
        column(
          3,
          actionButton("sc1a2togL", "Toggle to subset cells"),
          conditionalPanel(
            condition = "input.sc1a2togL % 2 == 1",
            selectInput(
              "sc1a2sub1",
              "Cell information to subset:",
              choices = sc1conf[grp == TRUE]$UI,
              selected = sc1def$grp1
            ),
            uiOutput("sc1a2sub1.ui"),
            actionButton(
              "sc1a2sub1all",
              "Select all groups",
              class = "btn btn-primary"
            ),
            actionButton(
              "sc1a2sub1non",
              "Deselect all groups",
              class = "btn btn-primary"
            )
          )
        ), # End of column (6 space)
        column(
          6,
          actionButton("sc1a2tog0", "Toggle graphics controls"),
          conditionalPanel(
            condition = "input.sc1a2tog0 % 2 == 1",
            fluidRow(
              column(
                6,
                sliderInput(
                  "sc1a2siz",
                  "Point size:",
                  min = 0,
                  max = 4,
                  value = 1.25,
                  step = 0.25
                ),
                radioButtons(
                  "sc1a2psz",
                  "Plot size:",
                  choices = c("Small", "Medium", "Large"),
                  selected = "Medium",
                  inline = TRUE
                ),
                radioButtons(
                  "sc1a2fsz",
                  "Font size:",
                  choices = c("Small", "Medium", "Large"),
                  selected = "Medium",
                  inline = TRUE
                )
              ),
              column(
                6,
                radioButtons(
                  "sc1a2asp",
                  "Aspect ratio:",
                  choices = c("Square", "Fixed", "Free"),
                  selected = "Square",
                  inline = TRUE
                ),
                checkboxInput("sc1a2txt", "Show axis text", value = FALSE)
              )
            )
          )
        ) # End of column (6 space)
      ), # End of fluidRow (4 space)
      fluidRow(
        column(
          6,
          style = "border-right: 2px solid black",
          fluidRow(
            column(
              6,
              selectInput(
                "sc1a2ass1",
                "Data type to colour plot:",
                choices = c(
                  "Cell Information",
                  paste0("Assay: ", sc1def$assay)
                ),
                selected = "Cell Information"
              ),
              selectInput(
                "sc1a2inp1",
                "Cell Info / Feature Name:",
                choices = NULL
              ) %>%
                helper(
                  type = "inline",
                  size = "m",
                  fade = TRUE,
                  title = "Cell Info / Gene to colour cells by",
                  content = c(
                    "Select cell info / feature to colour cells",
                    "- Categorical covariates have a fixed colour palette",
                    paste0(
                      "- Continuous covariates / gene expression are coloured ",
                      "in a Blue-Yellow-Red colour scheme, which can be ",
                      "changed in the plot controls"
                    )
                  )
                )
            ),
            column(
              6,
              actionButton("sc1a2tog1", "Toggle plot controls"),
              conditionalPanel(
                condition = "input.sc1a2tog1 % 2 == 1",
                radioButtons(
                  "sc1a2col1",
                  "Colour (Continuous data):",
                  choices = c(
                    "White-Red",
                    "Blue-Yellow-Red",
                    "Yellow-Green-Purple"
                  ),
                  selected = "Blue-Yellow-Red"
                ),
                numericInput(
                  "sc1a2min1",
                  "Min cutoff (q##):",
                  min = 0,
                  max = 50,
                  value = 0,
                  step = 1
                ),
                numericInput(
                  "sc1a2max1",
                  "Max cutoff (q##):",
                  min = 50,
                  max = 100,
                  value = 100,
                  step = 1
                ),
                radioButtons(
                  "sc1a2ord1",
                  "Plot order:",
                  choices = c("Max-1st", "Min-1st", "Original", "Random"),
                  selected = "Original",
                  inline = TRUE
                ),
                checkboxInput(
                  "sc1a2lab1",
                  "Show cell info labels",
                  value = TRUE
                )
              )
            )
          ),
          fluidRow(column(12, uiOutput("sc1a2oup1.ui"))),
          fluidRow(column(12, uiOutput("sc1a2oup3.ui"))),
          fluidRow(
            column(3, downloadButton("sc1a2oup1.dl", "Download Plot")),
            column(
              3,
              radioButtons(
                "sc1a2oup1.f",
                "Plot format:",
                choices = c("png", "pdf"),
                selected = "png",
                inline = TRUE
              )
            ),
            column(
              3,
              numericInput(
                "sc1a2oup1.h",
                "Plot height:",
                min = 4,
                max = 20,
                value = 8,
                step = 0.5
              )
            ),
            column(
              3,
              numericInput(
                "sc1a2oup1.w",
                "Plot width:",
                min = 4,
                max = 20,
                value = 8,
                step = 0.5
              )
            )
          ),
          fluidRow(
            column(3, downloadButton("sc1a2oup3.dl", "Download Legend")),
            column(
              3,
              radioButtons(
                "sc1a2oup3.f",
                "Legend format:",
                choices = c("png", "pdf"),
                selected = "png",
                inline = TRUE
              )
            ),
            column(
              3,
              numericInput(
                "sc1a2oup3.h",
                "Legend height:",
                min = 0.2,
                max = 10,
                value = 2,
                step = 0.1
              )
            ),
            column(
              3,
              numericInput(
                "sc1a2oup3.w",
                "Legend width:",
                min = 4,
                max = 20,
                value = 8,
                step = 0.5
              )
            )
          ),
          br()
        ), # End of column (6 space)
        column(
          6,
          fluidRow(
            column(
              6,
              selectInput(
                "sc1a2ass2",
                "Data type to colour plot:",
                choices = c(
                  "Cell Information",
                  paste0("Assay: ", sc1def$assay)
                ),
                selected = paste0("Assay: ", sc1def$assay[1])
              ),
              selectInput(
                "sc1a2inp2",
                "Cell Info / Feature Name:",
                choices = NULL
              ) %>%
                helper(
                  type = "inline",
                  size = "m",
                  fade = TRUE,
                  title = "Cell Info / Gene to colour cells by",
                  content = c(
                    "Select cell info / feature to colour cells",
                    "- Categorical covariates have a fixed colour palette",
                    paste0(
                      "- Continuous covariates / gene expression are coloured ",
                      "in a Blue-Yellow-Red colour scheme, which can be ",
                      "changed in the plot controls"
                    )
                  )
                )
            ),
            column(
              6,
              actionButton("sc1a2tog2", "Toggle plot controls"),
              conditionalPanel(
                condition = "input.sc1a2tog2 % 2 == 1",
                radioButtons(
                  "sc1a2col2",
                  "Colour (Continuous data):",
                  choices = c(
                    "White-Red",
                    "Blue-Yellow-Red",
                    "Yellow-Green-Purple"
                  ),
                  selected = "White-Red"
                ),
                numericInput(
                  "sc1a2min2",
                  "Min cutoff (q##):",
                  min = 0,
                  max = 50,
                  value = 0,
                  step = 1
                ),
                numericInput(
                  "sc1a2max2",
                  "Max cutoff (q##):",
                  min = 50,
                  max = 100,
                  value = 100,
                  step = 1
                ),
                radioButtons(
                  "sc1a2ord2",
                  "Plot order:",
                  choices = c("Max-1st", "Min-1st", "Original", "Random"),
                  selected = "Max-1st",
                  inline = TRUE
                ),
                checkboxInput(
                  "sc1a2lab2",
                  "Show cell info labels",
                  value = TRUE
                )
              )
            )
          ),
          fluidRow(column(12, uiOutput("sc1a2oup2.ui"))),
          fluidRow(column(12, uiOutput("sc1a2oup4.ui"))),
          fluidRow(
            column(3, downloadButton("sc1a2oup2.dl", "Download Plot")),
            column(
              3,
              radioButtons(
                "sc1a2oup2.f",
                "Plot format:",
                choices = c("png", "pdf"),
                selected = "png",
                inline = TRUE
              )
            ),
            column(
              3,
              numericInput(
                "sc1a2oup2.h",
                "Plot height:",
                min = 4,
                max = 20,
                value = 8,
                step = 0.5
              )
            ),
            column(
              3,
              numericInput(
                "sc1a2oup2.w",
                "Plot width:",
                min = 4,
                max = 20,
                value = 8,
                step = 0.5
              )
            )
          ),
          fluidRow(
            column(3, downloadButton("sc1a2oup4.dl", "Download Legend")),
            column(
              3,
              radioButtons(
                "sc1a2oup4.f",
                "Legend format:",
                choices = c("png", "pdf"),
                selected = "png",
                inline = TRUE
              )
            ),
            column(
              3,
              numericInput(
                "sc1a2oup4.h",
                "Legend height:",
                min = 0.2,
                max = 10,
                value = 2,
                step = 0.1
              )
            ),
            column(
              3,
              numericInput(
                "sc1a2oup4.w",
                "Legend width:",
                min = 4,
                max = 20,
                value = 8,
                step = 0.5
              )
            )
          ),
          br()
        ) # End of column (6 space)
      ), # End of fluidRow (4 space)
      h4("Relationship between left-side and right-side cell info / feature"),
      fluidRow(
        column(
          6,
          h4("Comparative plot"),
          style = "border-right: 2px solid black",
          uiOutput("sc1a2oup5.ui"),
          fluidRow(
            column(3, downloadButton("sc1a2oup5.dl", "Download Plot")),
            column(
              3,
              radioButtons(
                "sc1a2oup5.f",
                "Plot format:",
                choices = c("png", "pdf"),
                selected = "png",
                inline = TRUE
              )
            ),
            column(
              3,
              numericInput(
                "sc1a2oup5.h",
                "Plot height:",
                min = 4,
                max = 20,
                value = 8,
                step = 0.5
              )
            ),
            column(
              3,
              numericInput(
                "sc1a2oup5.w",
                "Plot width:",
                min = 4,
                max = 20,
                value = 8,
                step = 0.5
              )
            )
          )
        ),
        column(
          6,
          h4("Cell numbers / statistics"),
          numericInput("sc1a2cut", "Cutoff for Expression:", value = 0),
          dataTableOutput("sc1a2.dt")
        ) # End of column (6 space)
      ) # End of fluidRow (4 space)
    ), # End of tab (2 space)
    ### Tab1.a3: Gene coexpression plot
    tabPanel(
      HTML("Gene coexpression"),
      h4("Coexpression of two genes on reduced dimensions"),
      "In this tab, users can visualise the coexpression of two genes ",
      "on low-dimensional representions.",
      br(),
      br(),
      fluidRow(
        column(
          3,
          fluidRow(
            column(
              12,
              selectInput(
                "sc1a3dr",
                "Reduction:",
                choices = sc1def$dimrd,
                selected = sc1def$dimrd[1]
              )
            )
          )
        ), # End of column (6 space)
        column(
          3,
          actionButton("sc1a3togL", "Toggle to subset cells"),
          conditionalPanel(
            condition = "input.sc1a3togL % 2 == 1",
            selectInput(
              "sc1a3sub1",
              "Cell information to subset:",
              choices = sc1conf[grp == TRUE]$UI,
              selected = sc1def$grp1
            ),
            uiOutput("sc1a3sub1.ui"),
            actionButton(
              "sc1a3sub1all",
              "Select all groups",
              class = "btn btn-primary"
            ),
            actionButton(
              "sc1a3sub1non",
              "Deselect all groups",
              class = "btn btn-primary"
            )
          )
        ), # End of column (6 space)
        column(
          6,
          actionButton("sc1a3tog0", "Toggle graphics controls"),
          conditionalPanel(
            condition = "input.sc1a3tog0 % 2 == 1",
            fluidRow(
              column(
                6,
                sliderInput(
                  "sc1a3siz",
                  "Point size:",
                  min = 0,
                  max = 4,
                  value = 1.25,
                  step = 0.25
                ),
                radioButtons(
                  "sc1a3psz",
                  "Plot size:",
                  choices = c("Small", "Medium", "Large"),
                  selected = "Medium",
                  inline = TRUE
                ),
                radioButtons(
                  "sc1a3fsz",
                  "Font size:",
                  choices = c("Small", "Medium", "Large"),
                  selected = "Medium",
                  inline = TRUE
                )
              ),
              column(
                6,
                radioButtons(
                  "sc1a3asp",
                  "Aspect ratio:",
                  choices = c("Square", "Fixed", "Free"),
                  selected = "Square",
                  inline = TRUE
                ),
                checkboxInput("sc1a3txt", "Show axis text", value = FALSE)
              )
            )
          )
        ) # End of column (6 space)
      ), # End of fluidRow (4 space)
      fluidRow(
        column(
          3,
          style = "border-right: 2px solid black",
          h4("Assay Expression"),
          selectInput("sc1a3inp1", "Feature 1:", choices = NULL) %>%
            helper(
              type = "inline",
              size = "m",
              fade = TRUE,
              title = "Feature expression to colour cells by",
              content = c(
                "Select gene to colour cells by gene expression",
                paste0(
                  "- Feature expression are coloured in a ",
                  "White-Red colour scheme which can be ",
                  "changed in the plot controls"
                )
              )
            ),
          selectInput("sc1a3inp2", "Feature 2:", choices = NULL) %>%
            helper(
              type = "inline",
              size = "m",
              fade = TRUE,
              title = "Feature expression to colour cells by",
              content = c(
                "Select gene to colour cells by gene expression",
                paste0(
                  "- Feature expression are coloured in a ",
                  "White-Blue colour scheme which can be ",
                  "changed in the plot controls"
                )
              )
            ),
          selectInput(
            "sc1a3ass1",
            "Assay:",
            choices = sc1def$assay,
            selected = sc1def$assay[1]
          ),
          actionButton("sc1a3tog1", "Toggle plot controls"),
          conditionalPanel(
            condition = "input.sc1a3tog1 % 2 == 1",
            radioButtons(
              "sc1a3col1",
              "Colour:",
              choices = c(
                "Red (Gene1); Blue (Gene2)",
                "Orange (Gene1); Blue (Gene2)",
                "Red (Gene1); Green (Gene2)",
                "Green (Gene1); Blue (Gene2)"
              ),
              selected = "Red (Gene1); Blue (Gene2)"
            ),
            numericInput(
              "sc1a3min1",
              "Min cutoff gene 1 (q##):",
              min = 0,
              max = 50,
              value = 0,
              step = 1
            ),
            numericInput(
              "sc1a3max1",
              "Max cutoff gene 1 (q##):",
              min = 50,
              max = 100,
              value = 100,
              step = 1
            ),
            numericInput(
              "sc1a3min2",
              "Min cutoff gene 2 (q##):",
              min = 0,
              max = 50,
              value = 0,
              step = 1
            ),
            numericInput(
              "sc1a3max2",
              "Max cutoff gene 2 (q##):",
              min = 50,
              max = 100,
              value = 100,
              step = 1
            ),
            radioButtons(
              "sc1a3ord1",
              "Plot order:",
              choices = c("Max-1st", "Min-1st", "Original", "Random"),
              selected = "Max-1st",
              inline = TRUE
            )
          )
        ), # End of column (6 space)
        column(
          6,
          style = "border-right: 2px solid black",
          uiOutput("sc1a3oup1.ui"),
          fluidRow(
            column(3, downloadButton("sc1a3oup1.dl", "Download Plot")),
            column(
              3,
              radioButtons(
                "sc1a3oup1.f",
                "Plot format:",
                choices = c("png", "pdf"),
                selected = "png",
                inline = TRUE
              )
            ),
            column(
              3,
              numericInput(
                "sc1a3oup1.h",
                "Plot height:",
                min = 4,
                max = 20,
                value = 8,
                step = 0.5
              )
            ),
            column(
              3,
              numericInput(
                "sc1a3oup1.w",
                "Plot width:",
                min = 4,
                max = 20,
                value = 8,
                step = 0.5
              )
            )
          ),
          br()
        ), # End of column (6 space)
        column(
          3,
          uiOutput("sc1a3oup2.ui"),
          fluidRow(
            column(6, downloadButton("sc1a3oup2.dl", "Download Legend")),
            column(
              6,
              radioButtons(
                "sc1a3oup2.f",
                "Legend format:",
                choices = c("png", "pdf"),
                selected = "png",
                inline = TRUE
              )
            )
          ),
          br(),
          h4("Cell numbers"),
          dataTableOutput("sc1a3.dt")
        ) # End of column (6 space)
      ) # End of fluidRow (4 space)
    ), # End of tab (2 space)
    ### Tab1.b1: violinplot / boxplot
    tabPanel(
      HTML("Violinplot / Boxplot"),
      h4("Cell information / assay expression violin plot / box plot"),
      "In this tab, users can visualise the assay expression or continuous cell information ",
      "(e.g. Number of UMIs / module score) across groups of cells (e.g. libary / clusters).",
      br(),
      br(),
      fluidRow(
        column(
          3,
          style = "border-right: 2px solid black",
          selectInput(
            "sc1b1inp1",
            "Cell information (X-axis):",
            choices = sc1conf[grp == TRUE]$UI,
            selected = sc1def$grp1
          ) %>%
            helper(
              type = "inline",
              size = "m",
              fade = TRUE,
              title = "Cell information to group cells by",
              content = c(
                "Select categorical cell information to group cells by",
                "- Single cells are grouped by this categorical covariate",
                "- Plotted as the X-axis of the violin plot / box plot"
              )
            ),
          selectInput(
            "sc1b1inp2",
            "Cell Info / Feature name (Y-axis):",
            choices = NULL
          ) %>%
            helper(
              type = "inline",
              size = "m",
              fade = TRUE,
              title = "Cell Info / Gene to plot",
              content = c(
                "Select cell info / feature to plot on Y-axis",
                "- Can be continuous cell information (e.g. nUMIs / scores)",
                "- Can also be feature expression"
              )
            ),
          selectInput(
            "sc1b1ass1",
            "Data type for Y-axis:",
            choices = c("Cell Information", paste0("Assay: ", sc1def$assay)),
            selected = "Cell Information"
          ),
          radioButtons(
            "sc1b1typ",
            "Plot type:",
            choices = c("violin", "boxplot"),
            selected = "violin",
            inline = TRUE
          ),
          checkboxInput("sc1b1pts", "Show data points", value = FALSE),
          actionButton("sc1b1togT", "Toggle to perform stats test(s)"),
          conditionalPanel(
            condition = "input.sc1b1togT % 2 == 1",
            radioButtons(
              "sc1b1stg",
              "Perform global test:",
              choices = c("none", "kruskal.test", "anova"),
              selected = "none"
            ),
            radioButtons(
              "sc1b1stp1",
              "Perform pairwise test:",
              choices = c("none", "wilcox.test", "t.test"),
              selected = "none"
            ),
            textAreaInput(
              "sc1b1stp2",
              HTML(
                "Pairwise comparisons to make <br />  
                                             (Separate each pair by new line and <br />  
                                             within each pair using , or ;):"
              ),
              height = "150px",
              value = "none"
            )
          ),
          br(),
          br(),
          actionButton("sc1b1togL", "Toggle to subset cells"),
          conditionalPanel(
            condition = "input.sc1b1togL % 2 == 1",
            selectInput(
              "sc1b1sub1",
              "Cell information to subset:",
              choices = sc1conf[grp == TRUE]$UI,
              selected = sc1def$grp1
            ),
            uiOutput("sc1b1sub1.ui"),
            actionButton(
              "sc1b1sub1all",
              "Select all groups",
              class = "btn btn-primary"
            ),
            actionButton(
              "sc1b1sub1non",
              "Deselect all groups",
              class = "btn btn-primary"
            )
          ),
          br(),
          br(),
          actionButton("sc1b1tog", "Toggle graphics controls"),
          conditionalPanel(
            condition = "input.sc1b1tog % 2 == 1",
            sliderInput(
              "sc1b1siz",
              "Data point size:",
              min = 0,
              max = 4,
              value = 1.25,
              step = 0.25
            ),
            radioButtons(
              "sc1b1psz",
              "Plot size:",
              choices = c("Small", "Medium", "Large"),
              selected = "Medium",
              inline = TRUE
            ),
            radioButtons(
              "sc1b1fsz",
              "Font size:",
              choices = c("Small", "Medium", "Large"),
              selected = "Medium",
              inline = TRUE
            ),
            checkboxInput("sc1b1noi", "Add noise to assay expr", value = TRUE)
          )
        ), # End of column (6 space)
        column(
          9,
          uiOutput("sc1b1oup.ui"),
          fluidRow(
            column(2, downloadButton("sc1b1oup.dl", "Download Plot")),
            column(
              2,
              radioButtons(
                "sc1b1oup.f",
                "Plot format:",
                choices = c("png", "pdf"),
                selected = "png",
                inline = TRUE
              )
            ),
            column(
              2,
              numericInput(
                "sc1b1oup.h",
                "Plot height:",
                min = 4,
                max = 20,
                value = 8,
                step = 0.5
              )
            ),
            column(
              2,
              numericInput(
                "sc1b1oup.w",
                "Plot width:",
                min = 4,
                max = 20,
                value = 10,
                step = 0.5
              )
            )
          ),
          br()
        ) # End of column (6 space)
      ) # End of fluidRow (4 space)
    ), # End of tab (2 space)
    ### Tab1.b2: Proportion plot
    tabPanel(
      HTML("Proportion plot"),
      h4("Proportion / cell numbers across different cell information"),
      "In this tab, users can visualise the composition of single cells based on one discrete ",
      "cell information across another discrete cell information. ",
      "Usage examples include the library or cellcycle composition across clusters.",
      br(),
      br(),
      fluidRow(
        column(
          3,
          style = "border-right: 2px solid black",
          selectInput(
            "sc1b2inp1",
            "Cell information to plot (X-axis):",
            choices = sc1conf[grp == TRUE]$UI,
            selected = sc1def$grp2
          ) %>%
            helper(
              type = "inline",
              size = "m",
              fade = TRUE,
              title = "Cell information to plot cells by",
              content = c(
                "Select categorical cell information to plot cells by",
                "- Plotted as the X-axis of the proportion plot"
              )
            ),
          selectInput(
            "sc1b2inp2",
            "Cell information to group / colour by:",
            choices = sc1conf[grp == TRUE]$UI,
            selected = sc1def$grp1
          ) %>%
            helper(
              type = "inline",
              size = "m",
              fade = TRUE,
              title = "Cell information to group / colour cells by",
              content = c(
                "Select categorical cell information to group / colour cells by",
                "- Proportion / cell numbers are shown in different colours"
              )
            ),
          radioButtons(
            "sc1b2typ",
            "Plot value:",
            choices = c("Proportion", "CellNumbers"),
            selected = "Proportion",
            inline = TRUE
          ),
          checkboxInput("sc1b2flp", "Flip X/Y", value = FALSE),
          selectInput(
            "sc1b2ord1",
            "Reorder X-axis by which group:",
            choices = NULL
          ),
          radioButtons(
            "sc1b2ord2",
            "Reorder in which order:",
            choices = c("Decreasing", "Increasing"),
            selected = "Decreasing",
            inline = TRUE
          ),
          actionButton("sc1b2togL", "Toggle to subset cells"),
          conditionalPanel(
            condition = "input.sc1b2togL % 2 == 1",
            selectInput(
              "sc1b2sub1",
              "Cell information to subset:",
              choices = sc1conf[grp == TRUE]$UI,
              selected = sc1def$grp1
            ),
            uiOutput("sc1b2sub1.ui"),
            actionButton(
              "sc1b2sub1all",
              "Select all groups",
              class = "btn btn-primary"
            ),
            actionButton(
              "sc1b2sub1non",
              "Deselect all groups",
              class = "btn btn-primary"
            )
          ),
          br(),
          br(),
          actionButton("sc1b2tog", "Toggle graphics controls"),
          conditionalPanel(
            condition = "input.sc1b2tog % 2 == 1",
            radioButtons(
              "sc1b2psz",
              "Plot size:",
              choices = c("Small", "Medium", "Large"),
              selected = "Medium",
              inline = TRUE
            ),
            radioButtons(
              "sc1b2fsz",
              "Font size:",
              choices = c("Small", "Medium", "Large"),
              selected = "Medium",
              inline = TRUE
            )
          )
        ), # End of column (6 space)
        column(
          9,
          uiOutput("sc1b2oup.ui"),
          fluidRow(
            column(2, downloadButton("sc1b2oup.dl", "Download Plot")),
            column(
              2,
              radioButtons(
                "sc1b2oup.f",
                "Plot format:",
                choices = c("png", "pdf"),
                selected = "png",
                inline = TRUE
              )
            ),
            column(
              2,
              numericInput(
                "sc1b2oup.h",
                "Plot height:",
                min = 4,
                max = 20,
                value = 8,
                step = 0.5
              )
            ),
            column(
              2,
              numericInput(
                "sc1b2oup.w",
                "Plot width:",
                min = 4,
                max = 20,
                value = 10,
                step = 0.5
              )
            )
          ),
          br()
        ) # End of column (6 space)
      ) # End of fluidRow (4 space)
    ), # End of tab (2 space)
    ### Tab1.b3: Bubbleplot / Heatmap
    tabPanel(
      HTML("Bubbleplot / Heatmap"),
      h4("Gene expression bubbleplot / heatmap"),
      "In this tab, users can visualise the gene expression patterns of ",
      "multiple genes grouped by categorical cell information (e.g. library / cluster).",
      br(),
      "The normalised expression are averaged, log-transformed and then plotted.",
      br(),
      br(),
      fluidRow(
        column(
          3,
          style = "border-right: 2px solid black",
          textAreaInput(
            "sc1b3inp",
            HTML(
              "List of gene names <br /> 
                                        (Max 50 genes, separated <br /> 
                                         by , or ; or newline):"
            ),
            height = "250px",
            value = paste0(sc1def$genes[[1]], collapse = ", ")
          ) %>%
            helper(
              type = "inline",
              size = "m",
              fade = TRUE,
              title = "List of genes to plot on bubbleplot / heatmap",
              content = c(
                "Input genes to plot",
                "- Maximum 50 genes (due to ploting space limitations)",
                "- Genes should be separated by comma, semicolon or newline"
              )
            ),
          selectInput(
            "sc1b3ass1",
            "Assay:",
            choices = sc1def$assay,
            selected = sc1def$assay[1]
          ),
          selectInput(
            "sc1b3grp",
            "Group by:",
            choices = sc1conf[grp == TRUE]$UI,
            selected = sc1conf[grp == TRUE]$UI[1]
          ) %>%
            helper(
              type = "inline",
              size = "m",
              fade = TRUE,
              title = "Cell information to group cells by",
              content = c(
                "Select categorical cell information to group cells by",
                "- Single cells are grouped by this categorical covariate",
                "- Plotted as the X-axis of the bubbleplot / heatmap"
              )
            ),
          radioButtons(
            "sc1b3plt",
            "Plot type:",
            choices = c("Bubbleplot", "Heatmap"),
            selected = "Bubbleplot",
            inline = TRUE
          ),
          checkboxInput("sc1b3scl", "Scale gene expression", value = TRUE),
          checkboxInput("sc1b3row", "Cluster rows (genes)", value = TRUE),
          checkboxInput("sc1b3col", "Cluster columns (samples)", value = FALSE),
          sliderInput(
            "sc1b3max",
            "Scale.max:",
            min = 1,
            max = 11,
            value = 3,
            step = 0.1
          ),
          checkboxInput(
            "sc1b3exp",
            "expm1 before averaging (check for log-transformed data)",
            value = TRUE
          ),
          br(),
          actionButton("sc1b3togL", "Toggle to subset cells"),
          conditionalPanel(
            condition = "input.sc1b3togL % 2 == 1",
            selectInput(
              "sc1b3sub1",
              "Cell information to subset:",
              choices = sc1conf[grp == TRUE]$UI,
              selected = sc1def$grp1
            ),
            uiOutput("sc1b3sub1.ui"),
            actionButton(
              "sc1b3sub1all",
              "Select all groups",
              class = "btn btn-primary"
            ),
            actionButton(
              "sc1b3sub1non",
              "Deselect all groups",
              class = "btn btn-primary"
            )
          ),
          br(),
          br(),
          actionButton("sc1b3tog", "Toggle graphics controls"),
          conditionalPanel(
            condition = "input.sc1b3tog % 2 == 1",
            radioButtons(
              "sc1b3cols",
              "Colour scheme:",
              choices = c(
                "White-Red",
                "Blue-Yellow-Red",
                "Yellow-Green-Purple"
              ),
              selected = "Blue-Yellow-Red"
            ),
            radioButtons(
              "sc1b3psz",
              "Plot size:",
              choices = c("Small", "Medium", "Large"),
              selected = "Medium",
              inline = TRUE
            ),
            radioButtons(
              "sc1b3fsz",
              "Font size:",
              choices = c("Small", "Medium", "Large"),
              selected = "Medium",
              inline = TRUE
            )
          )
        ), # End of column (6 space)
        column(
          9,
          h4(htmlOutput("sc1b3oupTxt")),
          uiOutput("sc1b3oup.ui"),
          fluidRow(
            column(2, downloadButton("sc1b3oup.dl", "Download Plot")),
            column(
              2,
              radioButtons(
                "sc1b3oup.f",
                "Plot format:",
                choices = c("png", "pdf"),
                selected = "png",
                inline = TRUE
              )
            ),
            column(
              2,
              numericInput(
                "sc1b3oup.h",
                "Plot height:",
                min = 4,
                max = 20,
                value = 10,
                step = 0.5
              )
            ),
            column(
              2,
              numericInput(
                "sc1b3oup.w",
                "Plot width:",
                min = 4,
                max = 20,
                value = 10,
                step = 0.5
              )
            )
          ),
          br()
        ) # End of column (6 space)
      ) # End of fluidRow (4 space)
    ), # End of tab (2 space)
    tabPanel(
      "DEG Explorer",
      fluidRow(
        column(
          width = 12,
          h3("Differential Expression Explorer"),
          p(
            "Explore significant differential gene expression between AD and controls from MAST with donor, age and sex as covaraites. Enter a gene symbol to view log fold-changes across cell types, split by cell type level. Entering a gene not present will return nothing."
          ),
          br()
        )
      ),
      fluidRow(
        column(
          width = 4,
          textInput("gene_input", "Enter Gene Symbol:", value = "APOE")
        )
      ),
      fluidRow(
        column(
          width = 6,
          h4("Level 1 Cell Types"),
          plotOutput("deg_plot_lvl1")
        ),
        column(
          width = 6,
          h4("Level 2 Cell Types"),
          plotOutput("deg_plot_lvl2")
        )
      )
    ),
    tabPanel(
      "Pseudotime Explorer",
      fluidRow(
        column(
          width = 12,
          h3("Pseudotime Expression Explorer"),
          p(
            "Explore pseudotime gene expression between AD and controls from tradeSeq. Enter a gene symbol to view expression over pseudotime, coloured by lineage and AD/control status. Note that this is a subset built from just the Endothelial, Endo-MT, SMC and Pericyte celltypes. Also note that not all genes are present as this is subset to those that are significantly different across pseudotime, inputting genes that aren't present will return an error."
          ),
          br()
        )
      ),
      fluidRow(
        column(
          width = 4,
          textInput(
            "pseudotime_gene_input",
            "Enter Gene Symbol:",
            value = "ANGPT2"
          )
        )
      ),
      fluidRow(
        column(
          width = 6,
          h4("Pseudotime Gene Expression"),
          plotOutput("tradeseq_pseudotime_plot")
        )
      )
    ),

    br(),
    p("", style = "font-size: 125%;"),
    p(
      em("This webpage was made using "),
      a(
        "ShinyCell2",
        href = "https://github.com/the-ouyang-lab/ShinyCell2",
        target = "_blank"
      )
    ),
    br(),
    br(),
    br(),
    br(),
    br()
  )
))
