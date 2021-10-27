#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(FactoMineR)
library(ggplot2)


geneList <- c('ISL1','CPS1','COL19A1','LAMB4','CHGA','NECAB2','CTCFL','MAGEA10','SYT5','GDPD2','H19','MYCN','SCN2A','COL2A1','NTRK2','PLAAT5','SYT9','LGALS4','CALB2','OXTR','RPL10P6','ERVH481','LERFS','RTL9','LINC00942','KCNJ18')
geneName <- c('ISL1','CPS1','COL19A1','LAMB4','CHGA','NECAB2','CTCFL','MAGEA10','SYT5','GDPD2','H19','MYCN','SCN2A','COL2A1','NTRK2','PLAAT5','SYT9','LGALS4','CALB2','OXTR','RPL10P6','ERVH48-1','LERFS','RTL9','LINC00942','KCNJ18')
mean <- c(35.7601363309151,14566.5680653128,51.3841768355949,44.1745844827057,926.646845008521,74.9864323590232,37.5685378893246,127.130420105964,149.118891413212,225.882966137232,2780.30864701087,157.799898538429,505.450553491039,267.555120568223,6061.28015830747,120.616317380566,9.54297851946113,605.056201862847,1184.99421751597,129.569305756442,140.642598767128,73.0379545442,14.4720293783997,9.21999879062841,844.573785078086,1.80284491704349)
sd <- c(171.467748293296,35621.1799455326,186.683269511824,126.339642144537,4858.76867363721,117.374093007517,151.830224041199,385.096079125658,449.754173943387,1212.75088341809,12728.3080144216,449.880230102692,2198.17123763117,949.141847855096,17918.2332041063,258.424846518006,16.6252997674351,2487.42492737406,3003.6066244316,189.136522116684,334.890423171494,237.804481638616,32.1319632752683,22.5021158530028,2779.94507488914,5.79703643980779)
            
ui <- fluidPage(

    h1("Predictive Biomarker for Treatment of NSCLC with Carboplatin"),
    
    wellPanel(
        
        fluidRow(
                
            column(4,
                   h3("Patient is "),
                   textOutput("text")
            ,id='res'
            ),
    
            column(8,
                   plotOutput("pca")
            )
        )
        
    ),
    
    h3("Input count number of each genes: "),

    
        
        uiOutput("countinput"),
        


    
    tags$head(tags$style("#text{
                                 font-size: 30px;
                                 text-align: center;
                                 margin: auto;
                                 }
                         "
    )
    )
                

)

server <- function(input, output) {
    
    output$countinput <- renderUI({
        fluidRow(
            lapply(0:1, function(i) {
                column(2,
                       lapply((floor((length(geneList)/5))*i+1):(floor((length(geneList)/5))*(i+1)), function(j) {
                           numericInput(paste0(geneList[j],'_exp'), paste0(geneName[j],' count no.'),rnorm(1, mean = mean[j], sd = sd[j]))
                       })
                )
            }),
            
            lapply(0:3, function(i) {
                column(2,
                       lapply((floor((length(geneList)/6))*i+11):(floor((length(geneList)/6))*(i+1)+10), function(j) {
                           numericInput(paste0(geneList[j],'_exp'), paste0(geneName[j],' count no.'),rnorm(1, mean = mean[j], sd = sd[j]))
                       })
                )
            })
        )
    })

    output$text <- renderText({ 
            
            sum <- input$ISL1_exp*0.00172098340359184+input$CPS1_exp*1.56627466806865e-05+input$COL19A1_exp*0.00469705083328958+input$LAMB4_exp*-0.000412476911698773+input$CHGA_exp*-1.6854428482348e-06+input$NECAB2_exp*-0.00284721850367415+input$CTCFL_exp*0.000254983756192625+input$MAGEA10_exp*-0.000132565547821168+input$SYT5_exp*-0.000183294277391049+input$GDPD2_exp*7.47960852616267e-06+input$H19_exp*-2.05938734449301e-05+input$MYCN_exp*0.000179913696390075+input$SCN2A_exp*-2.6400560270463e-05+input$COL2A1_exp*1.49343424603225e-05+input$NTRK2_exp*-2.71503807370354e-05+input$PLAAT5_exp*0.00143521017791376+input$SYT9_exp*0.115450562703337+input$LGALS4_exp*2.17508001757279e-05+input$CALB2_exp*-6.88344489754432e-06+input$OXTR_exp*0.00185960170569203+input$RPL10P6_exp*-0.000668339619166987+input$ERVH481_exp*0.000218605059688095+input$LERFS_exp*-0.00878993835027868+input$RTL9_exp*0.0192116373364546+input$LINC00942_exp*-5.34718771652081e-05+input$KCNJ18_exp*0.0323343146942858-1.559557

            ifelse(sum<0, 'Carboplatin sensitive', 'Carboplatin resistant')
        
        })
    
    output$pca <- renderPlot({ 
        
        
        newdt <- read.table('writeforPCA_26gene_allsam.txt', header= TRUE, sep='\t', row.names = NULL)
        colnames(newdt)[22] <- 'ERVH481'
        newdt <- rbind(newdt,data.frame(ISL1=input$ISL1_exp,CPS1=input$CPS1_exp,COL19A1=input$COL19A1_exp,LAMB4=input$LAMB4_exp,CHGA=input$CHGA_exp,NECAB2=input$NECAB2_exp,CTCFL=input$CTCFL_exp,MAGEA10=input$MAGEA10_exp,SYT5=input$SYT5_exp,GDPD2=input$GDPD2_exp,H19=input$H19_exp,MYCN=input$MYCN_exp,SCN2A=input$SCN2A_exp,COL2A1=input$COL2A1_exp,NTRK2=input$NTRK2_exp,PLAAT5=input$PLAAT5_exp,SYT9=input$SYT9_exp,LGALS4=input$LGALS4_exp,CALB2=input$CALB2_exp,OXTR=input$OXTR_exp,RPL10P6=input$RPL10P6_exp,ERVH481=input$ERVH481_exp,LERFS=input$LERFS_exp,RTL9=input$RTL9_exp,LINC00942=input$LINC00942_exp,KCNJ18=input$KCNJ18_exp,treatment_best_response='prediction'))
        
        res.pca = PCA(newdt, scale.unit=TRUE, quali.sup=ncol(newdt), graph=F)
        plot.PCA(res.pca, axes=c(1, 2), choix="ind",label="none", habillage=ncol(newdt),invisible = "quali",
                 ggoptions=list(size=8))+
            theme(text = element_text(size = 16),
                  title = element_text(size = 16),
                  axis.title = element_text(size = 16),
                  axis.text = element_text(size =16),
                  plot.title = element_text(size = 20)
            )
        
    })
    

    
}

# Run the application 
shinyApp(ui = ui, server = server)
