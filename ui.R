ui <- dashboardPage(
  dashboardHeader(title="Thin Film Evaporator",
    dropdownMenu(type = "messages",
                 messageItem(
                   from = "Author",
                   message =HTML("Gianmarco Polotti"),
                   time = "21/07/2018"
                 ),
                 messageItem(
                   from = "Licence",
                   message =HTML("GPL2"),
                   time = "21/07/2018"
                 ),
                 messageItem(
                   from = "Version",
                   message =HTML(paste("Version: ","1.00", sep=" ")),
                   time = "21/07/2018"
                 ),
                 messageItem(
                   from = "Questions",
                   message = HTML(paste("You find help tabs on each page.","Click on them to get quick suggestions", sep="<br/>")),
                   icon = icon("question"),
                   time = "21/07/2018"
                 ),
                 messageItem(
                   from = "Disclaimer",
                   message =HTML(paste( "The Author has made every",
                                        "attempt to ensure the accuracy",
                                        "and reliability of the calculation.",
                                        "However, the results are provided",
                                        "'as is' withoutwarranty of any kind.",
                                        "The Author does not accept",
                                        "any reponsibility or liability",
                                        "for the accuracy,content,",
                                        "completeness, legality, or",
                                        "reliability ofthe calculation",
                                        "contained in this software."
                   ,sep="<br/>")),
                   icon = icon("life-ring"),
                   time = "21/07/2018"
                 )
    )
  ),
  dashboardSidebar(
    sidebarMenu(id="sidemenu",
                menuItem("Plant", tabName = "plant",icon = icon("th"),
                         menuSubItem("Draw", tabName = "draw",icon = shiny::icon("file-text")),
                         menuSubItem("Geometry", tabName = "geometry",icon = shiny::icon("file-text")),
                         menuSubItem("Process", tabName = "Process",icon = shiny::icon("file-text"))
                ),
                menuItem("Feed",tabName = "feed",icon = shiny::icon("line-chart"),
                         menuSubItem("Global", tabName = "feedglobal",icon = shiny::icon("file-text")),
                         menuSubItem("Composition", tabName = "feedcomp",icon = shiny::icon("file-text"))
                ),
                menuItem("Simulation", tabName = 'pca', icon = shiny::icon("braille"),
                         menuSubItem("Mass Flow", tabName ="MFsimul" ),
                         menuSubItem("Heat Flow", tabName ="HFsimul" ),
                         menuSubItem("Temperature", tabName ="Tsimul" ),
                         menuSubItem("Compositions", tabName = "Csimul"),
                         menuSubItem("Mass Transfer", tabName = "MTsimul"),
                         menuSubItem("Fluidynamic", tabName = "Fsimul"),
                         menuSubItem("Summary", tabName = "Summary")
                ),
                hr(),
                menuItem("Source code",icon = icon("file-code-o"),href = "https://github.com/polgia0/ThinLayerEvaporator"),
                menuItem("Manual", icon = icon("file-pdf-o"),href = "ThinFilm_manual.pdf")
    ),
    HTML('<br><br><br><br>'),
    HTML('<p><center>Further Help ? <br>Contact the developer at <font color="cyan"><br> gianmarco.polotti@gmail.com </font></center>')
  ),
  dashboardBody(
    tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")),  
    fluidRow(
      useShinyalert(),  # Set up shinyalert
      tabItems(
        tabItem(tabName = "draw",
                fluidPage(
                  tabsetPanel(type = "tabs",
                              tabPanel("Process Draw",
                                      br(),
                                      br(),
                                      fluidRow(column(12,offset=2,align="center",mainPanel(imageOutput("draw")))),
                                      br(),
                                      br(),
                                      br(),
                                      br(),
                                      br(),
                                      br(),
                                      column(12, offset=0,align="center",uiOutput("examplecombo"))
                              ),
                              tabPanel("Help", mainPanel(htmlOutput("hlp_draw")))
                  )
                )
        ),
        tabItem(tabName = "geometry",
          fluidPage(
              tabsetPanel(type = "tabs",
                          tabPanel("Geometrical Data",
                              br(),
                              br(),
                              fluidRow(
                                  column(12, offset=2,align="center",
                                      box(
                                        uiOutput("uiZ"),
                                        uiOutput("uiD"),
                                        background ="light-blue",width=6)
                                  ),
                                  column(10, offset=0,align="center",
                                         uiOutput("uiBn"),
                                         uiOutput("uiH0"),
                                         uiOutput("uiRot")
                                  )
                              )
                          ),
                              tabPanel("Help", mainPanel(htmlOutput("hlp_geometry")))
                          )
              )
        ),
        tabItem(tabName = "Process",
          fluidPage(
              tabsetPanel(type = "tabs",
                          tabPanel("Process Data",
                                   br(),
                                   br(),
                                   fluidRow(
                                     column(12, offset=2,align="center",
                                            box(
                                              uiOutput("uiP"),
                                              uiOutput("uiTH"),
                                              background ="light-blue",width=6)
                                     )
                                   )
                          ),
                          tabPanel("Help", mainPanel(htmlOutput("hlp_process")))
              )
          )
        ),
        tabItem(tabName = "feedglobal",
          fluidPage(
            tabsetPanel(type = "tabs",
                        tabPanel("Feed Data",
                                 br(),
                                 br(),
                                 fluidRow(
                                   column(12, offset=2,align="center",
                                          box(
                                            uiOutput("uiLFw"),
                                            uiOutput("uiKc"),
                                            uiOutput("uiRo"),
                                            uiOutput("uiCp"),
                                            uiOutput("uiMu"),
                                            background ="light-blue",width=6)
                                   ),
                                   column(10, offset=0,align="center",
                                          uiOutput("uiW0"),
                                          uiOutput("uiWF")
                                   )
                                 )
                        ),
                        tabPanel("Help", mainPanel(htmlOutput("hlp_feedglob")))
            )
          )
        ),
        tabItem(tabName = "feedcomp",
          fluidPage(
            titlePanel("Liquid Phase Composition"),                       
            tabsetPanel(type = "tabs",
                        tabPanel("Ingredients",
                                     br(),
                                     br(),
                                     column(12,offset=1,align="center",uiOutput("soleventcombo")),
                                     br(),
                                     column(12,offset=1,align="center",uiOutput("Nonsoleventcombo")
                                     )
                        ),
                        tabPanel("Help",mainPanel(htmlOutput("hlp_feedcomp")))
            )
          )
        ),
        tabItem(tabName="MFsimul",
          fluidPage(
            titlePanel("Mass Flow Simulation"),
            tabsetPanel(type = "tabs",
                        tabPanel("Plot",
                                 br(),
                                 br(),
                                 column(12,offset=1,
                                      box(
                                          radioButtons("MFbox",label=tags$b("Scale on X"),choices=list("linear"=1,"log10"=2),selected=2),
                                      background ="light-blue",width=3)
                                  ),
                                  br(),
                                  fluidRow(plotOutput("MFplot",height=600))
                        ),
                        tabPanel("Help",mainPanel(htmlOutput("hlp_MFsimul")))
            )
          )
        ),
        tabItem(tabName="HFsimul",
                fluidPage(
                  titlePanel("Heat Flow Simulation"),
                  tabsetPanel(type="tabs",
                              tabPanel("Plot",
                                       br(),
                                       br(),
                                       column(12,offset=1,
                                              box(
                                                radioButtons("HFbox",label=tags$b("Scale on X"),choices=list("linear"=1,"log10"=2),selected=2),
                                                background ="light-blue",width=3)
                                       ),
                                       br(),
                                       fluidRow(plotOutput("HFplot",height=600))
                              ),
                              tabPanel("Help",mainPanel(htmlOutput("hlp_HFsimul")))
                  )
                )
        ),
        tabItem(tabName="Tsimul",
                fluidPage(
                  titlePanel("Internal Temperature Simulation"),
                  tabsetPanel(type="tabs",
                              tabPanel("Plot",
                                       br(),
                                       br(),
                                       column(12,offset = 1,
                                              box(
                                                radioButtons("Tbox",label = tags$b("Scale on X"),choices=list("linear"=1,"log10"=2),selected=2),
                                                background ="light-blue",width=3)
                                       ),
                                       br(),
                                       fluidRow(plotOutput("Tplot",height=600))
                              ),
                              tabPanel("Help",mainPanel(htmlOutput("hlp_Tsimul")))
                  )
                )
        ),
        tabItem(tabName="Csimul",
                fluidPage(
                  titlePanel("Composition Simulation"),
                  tabsetPanel(type="tabs",
                              tabPanel("Plot",
                                       br(),
                                       br(),
                                       column(12,offset = 1,
                                              box(
                                                radioButtons("Cbox",label = tags$b("Scale on X"),choices = list("linear" = 1, "log10" = 2),selected = 2),
                                                background ="light-blue",width=3)
                                       ),
                                       br(),
                                       fluidRow(plotOutput("Cplot",height=600))
                              ),
                              tabPanel("Help",mainPanel(htmlOutput("hlp_Csimul")))
                  )
                )
        ),
        tabItem(tabName="MTsimul",
                fluidPage(
                  titlePanel("Mass Transfer Quantities"),
                  tabsetPanel(type="tabs",
                              tabPanel("Plot",
                                       br(),
                                       br(),
                                       column(12,offset = 1,
                                              box(
                                                radioButtons("MTbox",label = tags$b("Scale on X"),choices = list("linear" = 1, "log10" = 2),selected = 2),
                                                background ="light-blue",width=3)
                                       ),
                                       br(),
                                       fluidRow(plotOutput("MTplot",height=600))
                              ),
                              tabPanel("Help", mainPanel(htmlOutput("hlp_MTsimul")))
                  )
                )
        ),
        tabItem(tabName="Fsimul",
                fluidPage(
                  titlePanel("Fluidynamic Table"),
                  tabsetPanel(type = "tabs",
                              tabPanel("Values",
                                       br(),
                                       br(),
                                       fluidRow(mainPanel(DT::dataTableOutput("viewflu"),width=12))     
                              ),
                              tabPanel("Help", mainPanel(htmlOutput("hlp_Fsimul")))
                  )
                )
        ),
        tabItem(tabName="Summary",
                fluidPage(
                  titlePanel("Input/Output Table"),
                  tabsetPanel(type = "tabs",
                              tabPanel("Summary",
                                       br(),
                                       br(),
                                       fluidRow(mainPanel(DT::dataTableOutput("viewsum"),width=12))     
                              ),
                              tabPanel("Help", mainPanel(htmlOutput("hlp_summary")))
                  )
                )
        )

      )
    )
    
  )
  
)




