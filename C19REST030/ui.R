
library(shiny)
library(shinydashboard)
library(plotly)
library(shinyWidgets)
library(flexdashboard)

# Define UI for application
shinyUI(
  
  dashboardPage(
    dashboardHeader(title = 'C19REST 0.3.0',titleWidth=250), # title of the dashboard
    
    dashboardSidebar(width = 250,
                     tags$style("[type = 'number'] {font-size:15px;height:20px; padding: 0px;}"),
                     # skin of the sliders
                     chooseSliderSkin(
                       skin = c("Modern"), #,"Shiny", "Flat", "Modern", "Nice", "Simple", "HTML5", "Round", "Square"),
                       color = "#990000" # color of the sliders
                     ),
                      
                     p('Simulering af risiko for overbelægning 
                       af patienter med COVID-19
                       givet forskellige åbningsscenarier'
                       ,align="center",style='font-size: 15px;'), # font-family: "Georgia", Times, "Times New Roman", serif; 
                      
                        
                     selectInput("scenarie", label = "Scenarie", 
                                 choices = list("uden åbning" = 1, 
                                                "grundscenarie" = 2, 
                                                "grundscenarie halv" = 3,
                                                "grundscenarie normal" = 4,
                                                "grundscenarie + alle" = 5), 
                                 selected = 1),  
                   
                     sliderInput("ICUcap",
                                 label="Kapacitet af intensivpladser til C19",
                                 min = 100,
                                 max = 1000,
                                 value = 300),
                     
                     sliderInput("Hoscap",
                                 label="Kapacitet af hospitalspladser til C19",
                                 min = 100,
                                 max = 5000,
                                 value = 2000),

                     selectInput("nrep", label = "Antal bootstraps", 
                                 choices = list("10 - kun til funktionstest!" = 10, 
                                                "100" = 100), 
                                 selected = 1),
                     
                     sliderInput("BSsel",
                                 label="fraction of bootstraps to discard",
                                 min = 0,
                                 max = .9,
                                 value = .6,
                                 step = 0.01)

               

    ),
    
    dashboardBody(
      tags$head(tags$style(HTML('
                                                          .main-header .logo {
                                                          font-family: "Georgia", Times, "Times New Roman", serif;
                                                          font-weight: bold;
                                                          font-size: 24px;
                                                                              }'))), #custom font of 'title'
      tags$style(HTML("
                                          .tabbable > .nav > li > a {background-color: #990000 ;  color: white; border-color: #AD3333} 
                                          .tabbable > .nav > li[class=active] > a {background-color: #AD3333; color:#fff; border-color: #990000}
                                                ")), # colors of tabs
      tabsetPanel(type = 'tabs',
      #### tab model ####
      tabPanel('Model',
               h3('Model struktur'),
               h5( htmlOutput('teori_text')),
               h5( htmlOutput('flow_text')),
               div(img(src='flow.png',height='330',width='660'),
                   style='text-align: center;')
               
      ),
                  
                  #### Intensiv kapacitet ####
                  tabPanel('Intensiv kapacitet',icon = icon("bar-chart-o"),
                           tags$style(".small-box.bg-red { background-color: #990000 !important; color: #fff !important; }"), # color of value box 1 
                           tags$style(".box.bg-red { background-color: #990000 !important; color: #fff !important; }"), # color of red box                              # total width of sum of columns = 12
                              column(width = 10, plotOutput('plot3BS',height = "500px"),
                                     'Figur over risikoen for forskellige belægninggrader.
                                     Figuren er fremkommet ved at lave gentagne antal simuleringer
                                     (bootstraps) med forskellige parametre. Disse simuleringer
                                     vægtes efter hvor godt de passer til den observerede epidemi
                                     (Se -Start fasen- tab).
                                     Det grønne område indikerer at 10% af simuleringerne kom til
                                     dette belægningsantal på den givne dato. Det røde område indikerer
                                     at alle simuleringer kom over dette niveau.
                                     Den sorte linie (median) angiver det niveau 50% af simuleringer kom over.'), # column with plot
                                      
                               column(width = 2, align = "center")
                           
                           
                  ),
                  #### hospital kapacitet ####
                  tabPanel('Hospital kapacitet',icon = icon("bar-chart-o"),

                           column(width = 10, plotOutput('plot2BS',height = "500px"),
                                  'Figur over risikoen for forskellige belægninggrader.
                                     Figuren er fremkommet ved at lave gentagne antal simuleringer
                                     (bootstraps) med forskellige parametre. Disse simuleringer
                                     vægtes efter hvor godt de passer til den observerede epidemi
                                     (Se -Start fasen- tab).
                                     Det grønne område indikerer at 10% af simuleringerne kom til
                                     dette belægningsantal på den givne dato. Det røde område indikerer
                                     at alle simuleringer kom over dette niveau.
                                     Den sorte linie (median) angiver det niveau 50% af simuleringer kom over.' ),
                           column(width = 2, align = "center")
                           
                  ),
                  #### total smittede ####
                  tabPanel('Total smittede',icon = icon("bar-chart-o"),
                           plotOutput('plot1BS',height = "500px"),
                           'Figur over risikoen for forskellige epidemistørrelser 
                           (alle smittede på en given dato).
                                     Figuren er fremkommet ved at lave gentagne antal simuleringer
                                     (bootstraps) med forskellige parametre. Disse simuleringer
                                     vægtes efter hvor godt de passer til den observerede epidemi 
                                     (Se -Start fasen- tab).
                                     Det grønne område indikerer at 10% af simuleringerne kom til
                                     denne epidemistørrelse på den givne dato. Det røde område indikerer
                                     at alle simuleringer kom over dette niveau.
                                     Den sorte linie (median) angiver det niveau 50% af simuleringer kom over.'
                           
                  ),
                  
                  tabPanel('Start fasen',icon = icon("bar-chart-o"),
                           plotOutput('plotInit',height = "600px"),
                           
                           'Figur over risikoen for forskellige belægninggrader.
                                     Figuren er fremkommet ved at lave gentagne antal simuleringer
                                     (bootstraps) med forskellige parametre. Disse simuleringer
                                     vægtes efter hvor godt de passer til den observerede epidemi.
                                     Det grønne område indikerer at 10% af simuleringerne kom til
                                     dette belægningsantal på den givne dato. Det røde område indikerer
                                     at alle simuleringer kom over dette niveau.
                                     Den sorte linie (median) angiver det niveau 50% af simuleringer kom over.'
                           
                  ),
                  
                  #### bootstrap param ####
                  tabPanel('Parametre',icon = icon("table"),
                           column(2,
                                  h4("Bootstrap Param"),
                                  sliderInput("r.ER", label = h6("Skalering"), min = 10, 
                                              max = 190, value = c(50, 150), step = 1,post=" %"),
                                  sliderInput("r.InfInit1", label = h6("Antal smittede 11 marts under 60"), 
                                              min = 10000, max = 150000, value = c(40000, 70000), step = 10000),
                                  sliderInput("r.InfInit2", label = h6("Antal smittede 11 marts over 60"), 
                                              min = 1000, max = 30000, value = c(5000, 14000), step = 1000)
                           ),
                           column(2,
                                  h4("Under 60 år"),
                                  sliderInput("r.kE1", label = h6("Latenstid"), min = 1, 
                                              max = 20, value = c(4, 6), step = 1),
                                  sliderInput("r.k1", label = h6("Antal dage efter symptomer før indlæggelse"), min = 0.1, 
                                              max = 21, value = c(6, 10), step = 0.1),
                                  sliderInput("r.gI231", label = h6("Antal dage til ICU efter indlæggelse"), min = 0.1, 
                                              max = 21, value = c(0.5, 2.5), step = 0.1),
                                  sliderInput("r.pI1R1", label = h6("Andel [%] med sygdom der indlægges"), 
                                              min = 0.01, max = 2, value = c(0.05, 0.5), step = 0.01 ,post = " %"), #LAEC c(0.1, 0.5)
                                  sliderInput("r.pI2R1", label = h6("Sandsynlighed for rask efter indlæggelse"), 
                                              min = 0, max = 1, value = c(0.77, 0.97), step = 0.001)
                           ),
                           column(2,
                                  h4("......"),
                                  sliderInput("r.pI3R1", label = h6("Sandsynlighed for at overleve ICU"), 
                                              min = 0, max = 1, value = c(0.7, 0.95), step = 0.001),
                                  sliderInput("r.mu31", label = h6("Antal dage på ICU før bortgang"), min = 7, 
                                              max = 28, value = c(14, 28), step = 1),
                                  sliderInput("r.recI11", label = h6("Antal dage før rask udenfor hospital"), min = 0.1, 
                                              max = 21, value = c(4, 6), step = 0.1),
                                  sliderInput("r.recI21", label = h6("Antal dage på hospital før udskrivelse"), min = 0.1, 
                                              max = 21, value = c(5, 9), step = 0.1),
                                  sliderInput("r.recI31", label = h6("Antal dage på ICU før udskrivelse"), min = 7, 
                                              max = 28, value = c(14, 28), step = 0.1, )
                           ),
                           column(2,
                                  h4("Over 60 år"),
                                  sliderInput("r.kE2", label = h6("Latenstid"), min = 1, 
                                              max = 20, value = c(4, 6), step = 1),
                                  sliderInput("r.k2", label = h6("Antal dage efter symptomer før indlæggelse"), min = 0.1, 
                                              max = 21, value = c(5, 9), step = 0.1),
                                  sliderInput("r.gI232", label = h6("Antal dage til ICU efter indlæggelse"), min = 0.1, 
                                              max = 21, value = c(0.5, 1.5), step = 0.1),
                                  sliderInput("r.pI1R2", label = h6("Andel [%] med sygdom der indlægges"), 
                                              min = 1, max = 10, value = c(5, 6.2), step = 0.1 ,post = " %"),
                                  sliderInput("r.pI2R2", label = h6("Sandsynlighed for rask efter indlæggelse"), 
                                              min = 0, max = 1, value = c(0.7, 0.9), step = 0.001)
                           ),
                           column(2,
                                  h4("......"),
                                  sliderInput("r.pI3R2", label = h6("Sandsynlighed for at overleve ICU"), 
                                              min = 0, max = 1, value = c(0.45, 0.55), step = 0.001),
                                  sliderInput("r.mu32", label = h6("Antal dage på ICU før bortgang"), min = 7, 
                                              max = 28, value = c(14, 28), step = 1),
                                  sliderInput("r.recI12", label = h6("Antal dage før rask udenfor hospital"), min = 0.1, 
                                              max = 21, value = c(4, 6), step = 0.1),
                                  sliderInput("r.recI22", label = h6("Antal dage på hospital før udskrivelse"), min = 0.1, 
                                              max = 21, value = c(5, 15), step = 0.1),
                                  sliderInput("r.recI32", label = h6("Antal dage på ICU før udskrivelse"), min = 7, 
                                              max = 28, value = c(14, 28), step = 0.1)
                           )
      
                           
                           
                  ),
                  
                  #### beta param ####
                  tabPanel('Beta',icon = icon("table"),
                           column(2,
                                  h4("beta til 14 Apr - uden åbning"),
                                  sliderInput("r.betaWIlt", label = h6("Smitterate indenfor <60 aldersgruppe"), min = 0.01, 
                                              max = .5, value = c(.12, .22), step = 0.01),
                                  sliderInput("r.betaWIge", label = h6("Smitterate indenfor >=60 aldersgruppe"), min = 0.01, 
                                              max = .5, value = c(.14, .24), step = 0.01),
                                  sliderInput("r.betaB", label = h6("Smitterate mellem aldersgrupper"), min = 0.01, 
                                              max = .5, value = c(.02, .12), step = 0.01)
                                  
                           ),
                           column(2,
                                  h4("beta fra 14 Apr - grundscenarie"),
                                  sliderInput("r.betaWIlt1", label = h6("Smitterate indenfor <60 aldersgruppe"), min = 0.01, 
                                              max = .5, value = c(.16, .26), step = 0.01),
                                  sliderInput("r.betaWIge1", label = h6("Smitterate indenfor >=60 aldersgruppe"), min = 0.01, 
                                              max = .5, value = c(.13, .23), step = 0.01),
                                  sliderInput("r.betaB1", label = h6("Smitterate mellem aldersgrupper"), min = 0.01, 
                                              max = .5, value = c(.03, .13), step = 0.01)
                                  
                           ),
                           column(2,
                                  h4("beta fra 14 Apr - grundscenarie - halv fysisk afstand"),
                                  sliderInput("r.betaWIlt2", label = h6("Smitterate indenfor <60 aldersgruppe"), min = 0.01, 
                                              max = .5, value = c(.25, .35), step = 0.01),
                                  sliderInput("r.betaWIge2", label = h6("Smitterate indenfor >=60 aldersgruppe"), min = 0.01, 
                                              max = .5, value = c(.24, .34), step = 0.01),
                                  sliderInput("r.betaB2", label = h6("Smitterate mellem aldersgrupper"), min = 0.01, 
                                              max = .5, value = c(.06, .16), step = 0.01)
                                  
                           ),
                           column(2,
                                  h4("beta fra 14 Apr - grundscenarie - normal fysisk afstand"),
                                  sliderInput("r.betaWIlt3", label = h6("Smitterate indenfor <60 aldersgruppe"), min = 0.01, 
                                              max = .6, value = c(.34, .44), step = 0.01),
                                  sliderInput("r.betaWIge3", label = h6("Smitterate indenfor >=60 aldersgruppe"), min = 0.01, 
                                              max = .6, value = c(.33, .43), step = 0.01),
                                  sliderInput("r.betaB3", label = h6("Smitterate mellem aldersgrupper"), min = 0.01, 
                                              max = .6, value = c(.1, .2), step = 0.01)
                                  
                           ),
                           column(2,
                                  h4("beta fra 14 Apr - grundscenarie + alle"),
                                  sliderInput("r.betaWIlt4", label = h6("Smitterate indenfor <60 aldersgruppe"), min = 0.01, 
                                              max = 1.2, value = c(.19, .29), step = 0.01),
                                  sliderInput("r.betaWIge4", label = h6("Smitterate indenfor >=60 aldersgruppe"), min = 0.01, 
                                              max = 1, value = c(.19, .29), step = 0.01),
                                  sliderInput("r.betaB4", label = h6("Smitterate mellem aldersgrupper"), min = 0.01, 
                                              max = 1, value = c(.04, .14), step = 0.01)
                                  
                           )
                           
                           
                           
                  ),
                  #### Følsomhedsanalyse ####
                  tabPanel('SA',icon = icon("bar-chart-o"),
                           #plotOutput('plotSens',height = "400px")
                           tableOutput("plotSensP")
                           
                  ),
                  #### tab Parameter info ####
                  tabPanel('Par. info',
                           h3('Generelle parametre'),
                           strong('Skalering'), br(),
                           "Skalering på beta parametre.", br(),
                           strong('Smitteraten indenfor aldersgruppe'), br(),
                           "Specifik smitterate hhv. inden for grupperne '<60 år' og '>=60 år'", br(),
                           strong('Smitteraten mellem aldersgrupper'), br(),
                           "Specifik smitterate mellem de to grupper '<60 år' og '>=60 år'", br(),
                           h3('Parametre der beskriver karakteristika for de to aldersgrupper'), br(),
                           strong('Latenstid'), br(),
                           " Perioden fra en person bliver smittet til 
                                       personen bliver infektiøs.", br(),
                           strong('Antal dage efter smitte før indlæggelse'), br(),
                           "Hvor lang tid går der før indlæggelse. Når personen indlægges
                           antages at den ikke længere smitter", br(),
                           strong('Sandsynlighed for rask uden indlæggelse'), br(),
                           "Baseret på danske erfaringer indtil videre.", br(),
                           strong('Sandsynlighed for rask efter indlæggelse'), br(),
                           "Baseret på danske erfaringer indtil videre.", br(),
                           strong('Sandsynlighed for at overleve ICU'), br(),
                           "Baseret på danske erfaringer indtil videre.", br(),
                           strong('Antal dage til ICU efter indlæggelse'), br(),
                           "Baseret på danske erfaringer indtil videre.", br(),
                           h3("Referencer"),
                               "Referencer forefindes i notatet af 16. April."),
                  

                  #### tab kontakt ####
                  tabPanel('Kontakt',#h5( textOutput('Kontakt_text'))
                           strong('Hoved udvikler:'), br(),
                           'Seniorforsker Kaare Græsbøll (DTU Compute)',
                           br(),
                           strong('Udviklere:'), br(), 
                           'Lektor Lasse Engbo Christiansen (DTU), Lektor Uffe H. Thygesen (DTU Compute).', br(), 
                           'Modellen er kodet af Kaare Græsbøll. ', br(),
                           'Teknisk assistance af Carsten Kirkeby (KU Sund).', br(),
                           'Shiny Dashboard template udviklet af Elisabeth Ottesen Bangsgaard (DTU Compute)', br(),
                           strong('Yderligere input:'), br(),
                           'Parametrene er fastsat med input fra ekspertgruppen under SSI.
                           Resultaterne er diskuteret med samme ekspertgruppe.'
                           

                  )
                         
      ),
      tags$head(tags$style(HTML('
                                              .skin-blue .main-header .logo {
                                              background-color: #990000;
                                              }
                                              .skin-blue .main-header .logo:hover {
                                               background-color: #990000;
                                              }
                                              .skin-blue .main-header .navbar {
                                              background-color: #990000;
                                              } 
                                  '))) # change color of 'skin' to 'DTU-red'
      
      
    )
    
    
  )
  
)
