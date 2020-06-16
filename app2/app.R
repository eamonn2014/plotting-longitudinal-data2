#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Rshiny ideas from on https://gallery.shinyapps.io/multi_regression/
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list=ls())
set.seed(2345)
library(mvtnorm) 
library(rms)
library(ggplot2)
library(shiny) 
library(nlme)
library(MASS)
library(tidyverse)
library(shinyWidgets)
library(lme4)
library(DT)
library("shinyalert")

options(max.print=1000000)
fig.width <- 1300
fig.height <- 550
fig.height2 <- 450
library(shinythemes)        # more funky looking apps
p1 <- function(x) {formatC(x, format="f", digits=1)}
p2 <- function(x) {formatC(x, format="f", digits=2)}
options(width=100)

v4=90
plot.backgroundG <- colors()

# function to create longitudinal data
is.even <- function(x){ x %% 2 == 0 }


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ui <- fluidPage(theme = shinytheme("journal"), #https://www.rdocumentation.org/packages/shinythemes/versions/1.1.2
                # paper
                useShinyalert(),  # Set up shinyalert
                setBackgroundColor(
                    color = c( "#2171B5", "#F7FBFF"), 
                    gradient = "linear",
                    direction = "bottom"
                ),
                h2("Plotting longitudinal data and estimating averages with the help of a natural log transformation"),
                
                h4(p("A laboratory measurment is taken on patients over time. 
                Within visit windows the times at which the measurement is grouped for all patients into common visits. There is no intervention after baseline, so we will not be treating baseline as a covariate.
                We simulate data that has a skewed distribution. Many laboratory measurments have skewed distibutions as analyte amounts cannot be negative.
                        We describe the average
                      trend using simple summary statistics, log transform the data and calculate summary statistics followed by 
                      exponentiating back, we also fit a generalized least squares (GLS) model to the log transformed data
                      (using an autocorrelation structure of order 1 as we simulate AR(1)) to account for the fact patients are 
                      followed over time, finally back transforming the model estimates. We also fit a linear mixed model (LMM) to the data. Note summary statistics approaches ignore the 
                      longitudinal nature of the data. It cannot be emphasized enough, always plot the data.")),
            
                h4(p("The natural log transformation is used throughout. By anti-logging predictions from a model assuming a normal distribution on the logged values this will result in estimates of the median response. The first plot selection '1 Means calculated on untransformed data',


                is a plot of the raw data presenting at each timepoint arithmetic means and 95%CIs. We can see heavily skewed data and that this plot is not very illuminating. 
                 The next plot '1a Means calculated on untransformed data, antilog presentation' presents the same arithmetic mean and 95%CIs 
                 as the previous plot. But here the data are 
                 log transformed. Note the means and CIs are calculated on the untransformed data, 
                 then logged for presentation with antilog y tick values. This is not recommended 
                 but serves to show the untransformed arithmetic means more clearly presented than the first plot. 
                 Skip over '2 Medians calculated on untransformed data' which is not very illuminating.
                 The next selection '2a Medians calculated on untransformed data, antilog presentation' is generated using the same 
                 approach as that of figure 1a. That is, calculate medians and CIs; Log transform the data; Log transform medians and CIs. 
                 Plot the data replacing y-axis tick marks with antilog values. This approach is acceptable for presenting the median.
                    Next we have the best approach for calculation and presentation of the arthmetic means. The selection  
                    '3 Log transformation, calculate statistics then back transform (exponentiate)' is self explanatory.
                    Two final approaches are based on modelling. As our data
                    is longitudinal in nature these are the best approaches. 
                    First the data is logged, the model then fitted and the estimates then exponentiated. 
                    For presentation 
                    log the data, estimates and confidence bounds, plot and replace the y axis ticks with antilogs. 
                    The 6th selection shows a comparison of the approaches and the 7th selection diagnostics from the GLS model fit.  
                    
                     ")),
                
                
                h4(p("
                Note the inputs simulate the data for the analyses and the final simulated response is exponentiated. 
                     The second tab is for interest, we run linear mixed models on the simulated response  
                     and on the exponentiated simulated response. The third tab is just for fun, the first tab is repeated 
                     however here you can play around with inputs to adjust the plots.")),
                
                shinyUI(pageWithSidebar(
                    headerPanel(" "),
                    
                    sidebarPanel( width=3 ,
                                  tags$style(type="text/css", ".span8 .well { background-color: #00FFFF; }"),
                                  
                                  div(
                                      actionButton(inputId='ab1', label="Shiny",   icon = icon("th"), 
                                                   onclick ="window.open('https://raw.githubusercontent.com/eamonn2014/plotting-longitudinal-data2/master/app2/app.R', '_blank')"),   
                                      actionButton(inputId='ab1', label="R code",   icon = icon("th"), 
                                                   onclick ="window.open('https://raw.githubusercontent.com/eamonn2014/plotting-longitudinal-data2/master/plotting-longitudinal-data.R', '_blank')"),   
                                      actionButton("resample", "Simulate a new sample"),
                                      br(), br(),
                                      tags$style(".well {background-color:#b6aebd ;}"), 
                                      
                                      div(h5(tags$span(style="color:blue", "Play around with the parameters that generate 
                                                       data using the sliders below. Hit simulate for another sample. Hit the other two buttons to see the code."))),
                                      tags$head(
                                          tags$style(HTML('#ab1{background-color:orange}'))
                                      ),
                                      
                                      tags$head(
                                          tags$style(HTML('#resample{background-color:orange}'))
                                      ),
                                      
                                      
                                      selectInput("plot", div(h5(tags$span(style="color:blue", "Select a plot/estimation approach:"))),
                                                  list("1 Means calculated on untransformed data" = "plot1",
                                                       "1a Means calculated on untransformed data, antilog presentation" = "plot1x",
                                                       "2 Medians calculated on untransformed data" = "plot1a",
                                                       "2a Medians calculated on untransformed data, antilog presentation" = "plot1ax",
                                                       "3 Log transformation, calculate statistics then back transform (exponentiate)" = "plot2", 
                                                       "4 GLS model on natural log data, exponentiated estimates with 95% CI" = "plot3" ,
                                                       "5 LMM model on natural log data, exponentiated estimates with 95% CI" = "plot3a",
                                                       "6 Estimate comparison plot" = "plot4",
                                                       "7 GLS model diagnostics" = "plot5"
                                                  )),
                                      
                                      # selectInput("logplot", div(h5(tags$span(style="color:blue", "Log Transform plot?"))),
                                      #             list("transform" = "yes",
                                      #                  "no transformation" = "no" 
                                      #             )),
                                      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                      
                                      
                                      
                                      
                                      ### maximum number of possible observations
                                      
                                      sliderInput("N",
                                                  div(h5(tags$span(style="color:blue", "Total number of subjects"))),
                                                  min=20, max=500, step=5, value=100, ticks=FALSE),
                                      
                                      sliderInput("J",
                                                  div(h5(tags$span(style="color:blue", "Maximum visit in data simulation including baseline"))),
                                                  min=3, max=10, step=1, value=6, ticks=FALSE),
                                      
                                      
                                      sliderInput("autocorrelation",    
                                                  div(h5(tags$span(style="color:blue", "Auto correlation"))),  ###
                                                  min = 0, max = 1, value = c(.9), step=.1, ticks=FALSE),
                                      
                                      sliderInput("beta0", 
                                                  div(h5(tags$span(style="color:blue", "Average intercept"))),
                                                  min = -4, max = 4, step=.0001, value = c(1.1174), ticks=FALSE) ,
                                      
                                      sliderInput("beta1", 
                                                  div(h5(tags$span(style="color:blue", "Average slope"))),
                                                  min = -4, max = 4, step=.0001, value = c(-0.2859 ), ticks=FALSE) ,
                                      
                                      sliderInput("q",   
                                                  div(h5(tags$span(style="color:blue", "True intercept SD"))),
                                                  min = .0001, max = 4, step=.0001, value = c(1.2748), ticks=FALSE) ,
                                      
                                      sliderInput("s",      
                                                  div(h5(tags$span(style="color:blue", "True slope SD"))),
                                                  min = .0001, max = 4, step=.0001, value = c(0.2276), ticks=FALSE) ,
                                      
                                      sliderInput("r", 
                                                  div(h5(tags$span(style="color:blue", "True intercept slope correlation"))),
                                                  min = -1, max = 1, value = c(-0.62), step=0.05, ticks=FALSE),
                                      
                                      sliderInput("sigma",  
                                                  div(h5(tags$span(style="color:blue",  "True error SD"  ))),
                                                  min = .0001, max = 4,  value = c(0.7995 ), step=.1, ticks=FALSE),
                                      
                                      
                                      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                      
                                      div(p( strong("References:"))),  
                                      # 
                                      tags$a(href = "https://hbiostat.org/bbr/md/serial.html#fn4", "[1] Serial Data Frank Harrell"),
                                      div(p(" ")),
 
                                  )
                                  
                    ),
                    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~tab panels
                    mainPanel(width=9 ,
                              #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                              navbarPage(       
                                  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
                                  tags$style(HTML(" 
                            .navbar-default .navbar-brand {color: cyan;}
                            .navbar-default .navbar-brand:hover {color: blue;}
                            .navbar { background-color: #b6aebd;}
                            .navbar-default .navbar-nav > li > a {color:black;}
                            .navbar-default .navbar-nav > .active > a,
                            .navbar-default .navbar-nav > .active > a:focus,
                            .navbar-default .navbar-nav > .active > a:hover {color: pink;background-color: purple;}
                            .navbar-default .navbar-nav > li > a:hover {color: black;background-color:yellow;text-decoration:underline;}
                            .navbar-default .navbar-nav > li > a[data-value='t1'] {color: red;background-color: pink;}
                            .navbar-default .navbar-nav > li > a[data-value='t2'] {color: blue;background-color: lightblue;}
                            .navbar-default .navbar-nav > li > a[data-value='t3'] {color: green;background-color: lightgreen;}
     
                   ")), 
                                  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~end of section to add colour     
                                  tabPanel("Plotting longitudinal data & estimating averages with the help of a natural log transformation", 
                                           
                                           #h4(strong("Plot the data. Always plot the data.")),
                                           plotOutput("reg.plot"),
                                           h4(strong("Figure 1. Refer to the plot selection option for a description of the presentation")),
                                           
                                           div(class="span7", verbatimTextOutput("reg.summary5")),
                                           h4(strong("Table 1 Summary statistics, untransformed data")),
                                           
                                           
                                           div(class="span7", verbatimTextOutput("reg.summary3")),
                                               h4(strong("Table 2 Log transformation, calculate summary statistics then back transform (exponentiate)")),
                                               
                                               
                                               div(class="span7", verbatimTextOutput("reg.summary4")),
                                               h4(strong("Table 3 GLS model fit to natural logged data, the GLS model estimates and 95% CI are then exponentiated")),
                                               
                                               div(class="span7", verbatimTextOutput("reg.summary6")),
                                               h4(strong("Table 4 LMM model fit to natural logged data, the LMM model estimates and 95% CI are then exponentiated")),
                                           ),
                                           #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                           tabPanel("LMMs fitted to simulated data", value=3, 
                                                    
                                                    div(class="span7", verbatimTextOutput("reg.summary8")),
                                                    h4(paste("Table 5. LMM on simulated data, numeric time variable.")), 
                                                    div(class="span7", verbatimTextOutput("reg.summary10")),
                                                    h4(paste("Table 6. LMM on simulated data.")), 
                                                    div(class="span7", verbatimTextOutput("reg.summary9")),
                                                    h4(paste("Table 7. LMM on exponentiated simulated data. 
                                                             The estimates should be similar to those in selection '1 Means calculated on untransformed data', and Table 1")), 
                                           ) ,
                                           #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                  tabPanel("Plotting longitudinal data, playing with the plots",  
                                           
                                           splitLayout(
                                             
                                             # selectInput("Design",
                                             #             div(h5(tags$span(style="color:blue", "Select background colour:"))),
                                             #             choices=c(  plot.backgroundG ), width='100%',
                                             #             selected = "white"),
                                             
                                             textInput("v1", div(h5(tags$span(style="color:blue", "Line size: 0.1 -> "))), value= ".2"),
                                             textInput("v2", div(h5(tags$span(style="color:blue", "Line boldness: 0-1"))), value= ".8"),
                                             textInput("v3", div(h5(tags$span(style="color:blue", "Colour reverse: 1,-1"))), value= "1")
                                              
                                             # selectInput("Design",
                                             #             div(h5(tags$span(style="color:blue", "Select background colour:"))),
                                             #             choices=c(  plot.backgroundG ), width='22%',
                                             #             selected = "white")
                                             
                                           ),
                                           
                                           selectInput("Design",
                                                       div(h5(tags$span(style="color:blue", "Select background colour:"))),
                                                       choices=c(  plot.backgroundG ), width='22%',
                                                       selected = "white"),

                                           
                                           plotOutput("reg.plot2"),
                                           h4(strong("Figure 1. Refer to the plot selection option for a description of the presentation")),
                                           
                                           div(class="span7", verbatimTextOutput("reg.summary5a")),
                                           h4(strong("Table 1 Summary statistics, untransformed data")),
                                           
                                           div(class="span7", verbatimTextOutput("reg.summary3a")),
                                          h4(strong("Table 2 Log transformation, calculate summary statistics then back transform (exponentiate)")),
                                               
                                          div(class="span7", verbatimTextOutput("reg.summary4a")),
                                          h4(strong("Table 3 GLS model fit to natural logged data, the GLS model estimates and 95% CI are then exponentiated")),
                                          
                                          div(class="span7", verbatimTextOutput("reg.summary6a")),
                                          h4(strong("Table 4 LMM model fit to natural logged data, the LMM model estimates and 95% CI are then exponentiated")),
                                           
                                           
                                   ) 
                                           #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                  
                                  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                              )
                              #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~end tab panels 
                              
                    )
                )
                ))

server <- shinyServer(function(input, output   ) {
    
    
    #__________________________________________________________________________
    
    
    shinyalert("Welcome! \nLet's simulate longitudinal data and estimate averages",
               #"Hey Little Richard!", 
               type = "info")
    
    
    
    # --------------------------------------------------------------------------
    # This is where a new sample is instigated 
    random.sample <- reactive({
        
        # I analysed some skewed data and these are the parameters
        
        
        # Dummy line to trigger off button-press
        foo <-      input$resample
        
        n     <-    input$N
        beta0 <-  input$beta0 
        beta1 <-  input$beta1
        sigma <-  input$sigma
        ar.val <- input$autocorrelation
        tau0 <-   input$q
        tau1 <-   input$s
        tau01 <-  input$r
        m <-      input$J
    
        
        
        return(list( n=n,  beta0=beta0, beta1=beta1, sigma=sigma, 
                     tau0=tau0, tau1=tau1, tau01=tau01, ar.val=ar.val, m=m
        )) 
        
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    })
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    #  start creating data
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    
    make.data <- reactive({

        sample <- random.sample()

        n     <-  sample$n
        beta0 <-  sample$beta0
        beta1 <-  sample$beta1
        sigma <-  sample$sigma
        ar.val <- sample$ar.val
        tau0 <-   sample$tau0
        tau1 <-   sample$tau1
        tau01 <-  sample$tau01
        m <-      sample$m


        ### set number of individuals
       # n <- 100
       # beta0 <-  1.1174
       # beta1 <- -0.2859
       # ar.val <- 0.9
       # sigma <- 0.7995
       # tau0  <-  1.2748
       # tau1  <-  0.2276
       # tau01 <- -0.62

        ### maximum number of possible observations
       # m <- 6
        p <- round(runif(n,2,m))

        ### simulate observation moments (assume everybody has 1st obs)
        obs <- unlist(sapply(p, function(x) c(1, sort(sample(2:m, x-1, replace=FALSE)))))
        #obs = 10*(obs)
        ### set up data frame
        dat <- data.frame(id=rep(1:n, times=p), obs=obs)

        ### simulate (correlated) random effects for intercepts and slopes
        mu  <- c(0,0)
        S   <- matrix(c(1, tau01, tau01, 1), nrow=2)
        tau <- c(tau0, tau1)
        S   <- diag(tau) %*% S %*% diag(tau)
        U   <- mvrnorm(n, mu=mu, Sigma=S)

        ### simulate AR(1) errors and then the actual outcomes
        dat$eij <- unlist(sapply(p, function(x) arima.sim(model=list(ar=ar.val), n=x) * sqrt(1-ar.val^.2) * sigma))
        dat$yij <- (beta0 + rep(U[,1], times=p)) + (beta1 + rep(U[,2], times=p)) *  log(dat$obs) + dat$eij

        ### note: use arima.sim(model=list(ar=ar.val), n=x) * sqrt(1-ar.val^2) * sigma
        ### construction, so that the true error SD is equal to sigma

        dat$yij <- exp(dat$yij)   # here I exponenitate to get observed data


        ### create grouped data object
        dat <- groupedData( (yij) ~ obs | id, data=dat)

        dat$eij <- NULL
        names(dat )[names(dat) == "yij"]<-c("value")
        names(dat )[names(dat) == "obs"]<-c("VISIT")
        names(dat )[names(dat) == "id"]<-c("ID")
        dat$variable <- "BIOCHEM.B"

        return(list(d=dat))

    })

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
 
   # run models
   check <- reactive({
        
      data <- make.data()
      
      dat <- data$d
   
    #  f0 <-  lmer( (value) ~       VISIT * 1 + (1 + VISIT | ID), data = dat)
      
    #  f1 <-  lmer( log(value) ~    VISIT * 1 + (1 + VISIT | ID), data = dat)
      f <- tryCatch(lmer(log(value) ~    VISIT + (1 + as.numeric(VISIT) | ID), data = dat), 
                     error=function(e) e)
      
      dat$VISIT <- factor(dat$VISIT)
      
      f0 <- tryCatch(lmer(value ~    VISIT + 0+ (1 + as.numeric(VISIT) | ID), data = dat), 
                          error=function(e) e)
      
      f1 <- tryCatch(lmer(log(value) ~    VISIT + 0+ (1 + as.numeric(VISIT) | ID), data = dat), 
                     error=function(e) e)
      
      return(list(f0=f0, f1=f1, f=f))
        
    })
   
 
    
    
    
#     make.data <- reactive({
#         
#         sample <- random.sample()
#     # lets try different approach to simulation 
#  
#     n     <-  sample$n
#     beta0 <-  sample$beta0 
#     beta1 <-  sample$beta1
#     error <-  sample$sigma
#     #ar.val <- sample$ar.val
#     q <-      sample$tau0
#     s <-      sample$tau1
#     r <-      sample$tau01
#     J <-      sample$m
#     
#   
#     # inputs
#     dat <- flat.df <- NULL
#      
#     n <- 100                 # total patients
#     J <- 8                   # maximum visit of each patient
#     intercept   =  1.1174
#     slope       = -0.2859
#     trt         =  0 
#     error       =  0.7995
#     interaction =  0
#     # random effects parameters
#     q =  1.2748  # standard deviations for the intercept 
#     r = -0.62    # random effects correlation of slope internet
#     s =  0.2276  # standard deviations for slope
#     
#     
#     
#     
#     
#     #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     # set up patients and random treatment assignment
#     
#     unit.df <- data.frame(unit = c(1:N), treat = sample(1:1,   N, replace=TRUE) )  
#     #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     unit.df <-  within(unit.df, {
#         E.alpha.given.treat <-  intercept + 0      # trt is the true treatment effect
#         E.beta.given.treat  <-  slope + interaction * treat   # interaction effect
#     })
#     #
#     
#     (unit.df)[1:40,]
#     #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     # covariance matrix for random effects
#     
#     cov.matrix <- matrix(c(q^2, r * q * s, r * q * s, s^2), nrow = 2, byrow = TRUE)
#     
#     random.effects <- rmvnorm(N, mean = c(0, 0), sigma = cov.matrix)
#     
#     # add random effects
#     unit.df$alpha <- unit.df$E.alpha.given.treat + random.effects[, 1]
#     unit.df$beta <-  unit.df$E.beta.given.treat  + random.effects[, 2]
#     (unit.df)[1:40,]
#     
#     #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     # structure of experimental design, time points
#     
#     #x.grid = seq(0, 8, by = 8/J)[0:8]
#     
#     #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     # set up unbalanced visits everyone at first visit but randomly end up to max after that
#     
#     p <- round(runif(N,2,J))            # last visit for each person
#     
#     M= sum(p)
#     
#     unit <- sort(rep(c(1:N), times=p))  # id for each person
#     
#     j <- as.vector(unlist(tapply(X=unit, INDEX=list( unit), FUN=seq_along))) # count with each person
#     
#     time <- j-1  # first visit , baseline becomes 0
#     
#     within.unit.df <- data.frame(cbind(unit, j, time))  # create a data frame
#     
#     #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~end unbalanced
#     # merge design and dataand create response
#     
#     flat.df = merge(unit.df, within.unit.df)
#     
#     flat.df <-  within(flat.df, y <-  alpha + time * beta + error * rnorm(n = M))
#     
#     flat.df$treat <- as.factor(flat.df$treat)
#     flat.df$treat <- ifelse(flat.df$treat %in% 1, "Active","Placebo" )
#     flat.df[1:60,]
#     
#     
#     dat <- flat.df
#     names(dat )[names(dat) == "y"]<-c("value")
#     names(dat )[names(dat) == "j"]<-c("VISIT")
#     names(dat )[names(dat) == "unit"]<-c("ID")
#     
#     
#     dat <- dat[, c("value","VISIT","ID")]
#     dat <- groupedData( (value) ~ VISIT | ID, data=dat)
#     dat$variable <- "BIOCHEM.B"
#     
#     pd <- position_dodge(.4)
#     
#     dat$value <-  exp(dat$value)
# 
#     plot1 <-  ggplot(dat,   aes (x = VISIT, y =    (value), group = ID, color = ID)) +
#         geom_line() + geom_point() + ylab("response") + xlab("visit") +
#        stat_summary(fun=mean,geom="line", colour="black",lwd=1,aes(group=variable ) ) +
#         # geom_smooth(method=lm, se=FALSE, fullrange=TRUE )+
#         # scale_shape_manual(values=c(3, 16))+
#         #scale_color_manual(values=c('#999999','#E69F00'))+
#         theme(legend.position="top") +
#         #xlim(0, J) +
#         scale_x_continuous(breaks=c(0:J)) +
#         theme(panel.background=element_blank(),
#               # axis.text.y=element_blank(),
#               # axis.ticks.y=element_blank(),
#               # https://stackoverflow.com/questions/46482846/ggplot2-x-axis-extreme-right-tick-label-clipped-after-insetting-legend
#               # stop axis being clipped
#               plot.title=element_text(), plot.margin = unit(c(5.5,12,5.5,5.5), "pt"),
#               legend.text=element_text(size=12),
#               legend.title=element_text(size=14),
#               legend.position="none",
#               axis.text.x  = element_text(size=10),
#               axis.text.y  = element_text(size=10),
#               axis.line.x = element_line(color="black"),
#               axis.line.y = element_line(color="black"),
#               plot.caption=element_text(hjust = 0, size = 7))
# 
#     plot1

#     return(list(d=dat))
#     
# }) 

    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    
    # --------------------------------------------------------------------------
    # --------------------------------------------------------------------------
    # Fit the specified regression model
    fit.regression <- reactive({
        
        data <- make.data()
        
        df <- data$d
        
        require(rms)
        require(dplyr)
        
        df$x <- df$VISIT
        df$y <- log(df$value)
        df$g <- df$ID
        
        
        dfs<- df %>%
            group_by(x) %>%
            summarise(mu=mean(y,na.rm=T), sd=sd(y,na.rm=T), n=length(na.omit(y)), se=sd(y,na.rm=T)/sqrt(length(na.omit(y))),
                      lo = mean(y,na.rm=T)-2*sd(y,na.rm=T)/sqrt(length(na.omit(y))),
                      hi = mean(y,na.rm=T)+2*sd(y,na.rm=T)/sqrt(length(na.omit(y)) ))
        
        df1 <- merge(df,dfs,  by = intersect("x", "x"))
        
        f <-df1[,c("x","y","ID")]
        
        
        f$g <- (f$ID)
        f$x <- as.factor(f$x)
        f<-plyr::arrange(f, g, x)
        
        g <- f$g
        g <- as.vector(g)
        f$id<-rep(seq_along( rle(g)$values ), times = rle(g)$lengths )
        
        
        dd <<- datadist(f)
        options(datadist='dd')
        ################################
        
        # cp <- list(corAR1, corCAR1, corExp, corLin, corGaus,corSpher)#, corSymm)
        # 
        # z <- vector('list', length(cp))
        # 
        # for(k in 1:length(cp)) {
        # 
        # 
        #   z[[k]] <-  gls(y ~ x ,
        #                  correlation=cp[[k]](form=~ as.numeric(x)|id),
        #                  weights=varIdent(form=~1|x),
        #                  f, na.action=na.omit)
        # 
        # 
        # }
        # 
        # cat("lowest AIC?\n")
        # 
        # anova(z[[1]],z[[2]], z[[3]], z[[4]], z[[5]],  z[[6]])#), z[[7]])
        
        
        # Conditionally fit the model
        # if (input$model == "Generalized Least Squares on log transformed response") {
        
        fit.res <-  
            tryCatch(gls(y ~ x +0 ,
                         correlation=corAR1(form=~ as.numeric(x)|id), 
                         weights=varIdent(form=~1|x),
                         f, na.action=na.omit) , 
                     error=function(e) e)
        
        
        i <-  intervals(fit.res , which='coef')
        se <-  sqrt(diag(vcov(fit.res)))
        i <-  cbind(as.data.frame(i$coef), se)
        i <-  i[,c(2,4,1,3)]  
        
        
        estx <- as.data.frame(i)
        
        est <- exp(estx)            # exponentiate back
        est$x<-as.numeric(1:input$J)
        rownames(est) <- NULL
        est <-  est[,c(5,1,2,3,4)]  
        names(est) <- c('Visit','Estimate','se', 'Lower', 'Upper')
        #  }
        
        
        
        
        return(list(fit.res1=est , fit.res=fit.res ))
    })     
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # LMM
    fit.regression2 <- reactive({
        
        data <- make.data()
        
        df <- data$d
        
        require(lme4)
        require(dplyr)
        
        df$x <- df$VISIT
        df$y <- log(df$value)
        df$g <- df$ID
        
        
        dfs<- df %>%
            group_by(x) %>%
            summarise(mu=mean(y,na.rm=T), sd=sd(y,na.rm=T), n=length(na.omit(y)), se=sd(y,na.rm=T)/sqrt(length(na.omit(y))),
                      lo = mean(y,na.rm=T)-2*sd(y,na.rm=T)/sqrt(length(na.omit(y))),
                      hi = mean(y,na.rm=T)+2*sd(y,na.rm=T)/sqrt(length(na.omit(y)) ))
        
        df1 <- merge(df,dfs,  by = intersect("x", "x"))
        
        f <-df1[,c("x","y","ID")]
        
        
        f$g <- (f$ID)
        f$x <- as.factor(f$x)
        f<-plyr::arrange(f, g, x)
        
        g <- f$g
        g <- as.vector(g)
        f$id<-rep(seq_along( rle(g)$values ), times = rle(g)$lengths )
        
        ###############################
        
        fit.res <- tryCatch(lmer(y ~    x + 0+ (1 + as.numeric(x) | id), data = f), 
                            error=function(e) e)
        
        i <-  confint(fit.res )
        ci <- i[grep("^x", rownames(i)), ]
        fc <- fixef(fit.res)
        
        se <-  fixef(fit.res)  # just filling this , not needed
        i <-  cbind(as.data.frame(fc), se, ci)
        
        estx <- as.data.frame(i)
        
        est <- exp(estx)            # exponentiate back
        est$x<-as.numeric(1:input$J)
        rownames(est) <- NULL
        est <-  est[,c(5,1,2,3,4)]  
        names(est) <- c('Visit','Estimate','se', 'Lower', 'Upper')
        
        #est$se <- NULL
        return(list(fit.res2a=est , fit.res2b=fit.res ))
    })     
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    # treatment effect estimate
    data.summ <- reactive({
        
        # Get the current regression data
        data <- make.data()
        
        df1 <- data$d
        
        require(dplyr)
        require(knitr)
        
        
        # summary stats on raw data
        df_summaryx <- df1 %>%                            # the names of the new data frame and the data frame to be summarised
            group_by(VISIT, variable) %>%                 # the grouping variable
            summarise( n_PL = length(na.omit(value)),     # calculates the sample size per group
                       mean_PL = mean(value, na.rm=TRUE), # calculates the mean of each group
                       sd_PL = sd(value, na.rm=TRUE),      # calculates the sd of each group
                       SE_PL = sd(value, na.rm=TRUE)/sqrt(length(na.omit(value))), # SE of each group
                       Ns=length(unique(ID)),              # unique subjects
                       low = mean(value, na.rm=T) - qt(.975, length(na.omit(value))-1)*(sd(value, na.rm=T)/ length(na.omit(value))^.5),
                       upp = mean(value, na.rm=T) + qt(.975, length(na.omit(value))-1)*(sd(value, na.rm=T)/ length(na.omit(value))^.5),
                       median=median(value, na.rm=T) ,
                       min = min(value, na.rm=T),  
                       max = max(value, na.rm=T) ,
                       miss = sum(is.na(value))
            ) 
        
        food2 <- data.frame(df_summaryx)
        names(food2) <- c("Visit","Variable","Data points","Mean", "SD","SE","Subjects", "lower95%CI", "Upper95%CI", "Med", "Min","Max", "Missing")  
        print(kable(food2, format="pandoc", row.names = FALSE, digits = c(3)))
        return(list( summ=food2)) 
        
    })     
    
    
    
    ## log the data , calculate stats and back transform
    data.summ2 <- reactive({
        
        # Get the current regression data
        data <- make.data()
        
        df1 <- data$d
        
        df1$value <- log(df1$value)
        
        require(dplyr)
        require(knitr)
        
        df_summaryx <- df1 %>%                            # the names of the new data frame and the data frame to be summarised
            group_by(VISIT, variable) %>%                 # the grouping variable
            summarise( n_PL = length(na.omit(value)),     # calculates the sample size per group
                       mean_PL = mean(value, na.rm=TRUE), # calculates the mean of each group
                       sd_PL = sd(value, na.rm=TRUE),      # calculates the sd of each group
                       SE_PL = sd(value, na.rm=TRUE)/sqrt(length(na.omit(value))), # SE of each group
                       Ns=length(unique(ID)),              # unique subjects
                       low = mean(value, na.rm=T) - qt(.975, length(na.omit(value))-1)*(sd(value, na.rm=T)/ length(na.omit(value))^.5),
                       upp = mean(value, na.rm=T) + qt(.975, length(na.omit(value))-1)*(sd(value, na.rm=T)/ length(na.omit(value))^.5),
                       median=median(value, na.rm=T) ,
                       min = min(value, na.rm=T),  
                       max = max(value, na.rm=T) ,
                       miss = sum(is.na(value))
            ) 
        
        # food2 <- data.frame(df_summaryx)
        # names(food2) <- c("Visit","Variable","Data points","Mean", "SD","SE","No of Subjects", "lower95%CI", "Upper95%CI", "Med", "Min","Max", "Missing") 
        # food2 <- data.frame( food2[,c(1,2,3,7,13)],  apply(food2[,c(4,5,6,8,9,10,11,12)],2, exp))
        # food2 <- food2[, c(1,2,3, 6,4,7,8,9,10,11,12,13,5)]
        # print(kable(food2, format="pandoc", row.names = FALSE, digits = c(3)))
        # return(list( summ=food2)) 
        
        food2 <- data.frame(df_summaryx)
        
        food2[c(4,5,6,8,9,10,11,12)] <- lapply(food2[c(4,5,6,8,9,10,11,12)], function(x) exp(x))
        
        names(food2) <- c("Visit","Variable","Data points","Mean", "SD","SE","Subjects", "lower95%CI", "Upper95%CI", "Med", "Min","Max", "Missing")  
        print(kable(food2, format="pandoc", row.names = FALSE, digits = c(3)))
        return(list( summ=food2)) 
        
    })     
    
    #---------------------------------------------------------------------------
    # Plot the data  
    
    output$reg.plot2 <- output$reg.plot <- renderPlot({         
        
        require(ggplot2)
        
      v1 <- as.numeric(    eval(parse(text= (input$v1)) ) )
      v2 <- as.numeric(    eval(parse(text= (input$v2)) ) )
      v3 <- as.numeric(    eval(parse(text= (input$v3)) ) )
    #  v4 <- as.numeric(    eval(parse(text= (input$v4)) ) )
      
    
      
        # Get the current regression data
        data <- make.data()
        df <- data$d
        
        # get the model fit 
        f <- fit.regression()
        est <- f$fit.res1
        
        lmm.est <- fit.regression2()$fit.res2a
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        
        if (input$plot == "plot1") {   #MEANS
            
            
            df_summary <- df %>% # the names of the new data frame and the data frame to be summarised
                group_by(VISIT, variable) %>%                # the grouping variable
                summarise(mean_PL = mean(value, na.rm=TRUE),  # calculates the mean of each group
                          sd_PL = sd(value, na.rm=TRUE),      # calculates the sd of each group
                          n_PL = length(na.omit(value)),      # calculates the sample size per group
                          SE_PL = sd(value, na.rm=TRUE)/sqrt(length(na.omit(value)))) # SE of each group
            
            df_summary1 <- merge(df, df_summary)  # merge stats to dataset
            
            df_summary1$L2SE <- df_summary1$mean_PL - 2*df_summary1$SE_PL
            df_summary1$H2SE <- df_summary1$mean_PL + 2*df_summary1$SE_PL
            
            
            pr1 <- ggplot((df_summary1), aes(x = VISIT, y =value, color = ID)) +
                geom_line( size=v1, alpha=v2) +
                #scale_color_gradient(low = "blue", high = "red")+
                #scale_color_brewer(palette = "Dark2") +
                stat_summary(geom="line",  fun=mean, colour="black", lwd=0.5) +  # , linetype="dashed"
                stat_summary(geom="point", fun=mean, colour="black") +
                geom_errorbar(data=(df_summary1), 
                              aes( ymin=L2SE, ymax=H2SE ), color = "black",
                              width=0.05, lwd = 0.05) +
                scale_y_continuous(expand = c(.1,0) ) +
                
                
                
                scale_x_continuous(breaks = c(unique(df$VISIT)),
                                   labels = 
                                       c(unique(df$VISIT))
                ) +
                
                
                
                
               # geom_segment(aes(x = 1, xend = 3, y = (1), yend=(1)), color = "blue" , size=0.05, linetype="dashed", alpha=0.02) +
            #    geom_segment(aes(x = 1, xend = 3, y = (2.5), yend=(2.5)), color = "blue" , size=0.05, linetype="dashed", alpha=0.02) +
                
                theme(
                    # get rid of panel grids
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    # Change plot and panel background
                    plot.background=element_rect(fill = input$Design),                ##change background colour here
                    panel.background = element_rect(fill = 'black'),
                    # Change legend
                    legend.position = c(0.6, 0.07),
                    legend.direction = "horizontal",
                    legend.background = element_rect(fill = "black", color = NA),
                    legend.key = element_rect(color = "gray", fill = "black"),
                    legend.title = element_text(color = "white"),
                    legend.text = element_text(color = "white")
                ) +
                
                # theme(axis.text.y   = element_text(size=10),
                #       axis.text.x   = element_text(size=10),
                #       axis.title.y  = element_text(size=14),
                #       axis.title.x  = element_text(size=14),
                #       panel.background = element_blank(),
                #       panel.grid.major = element_blank(), 
                #       panel.grid.minor = element_blank(),
                #       legend.position="none",
                #       legend.text=NULL,
                #       legend.title=NULL,
            #       axis.line = element_line(colour = "black", size=0.05),
            #       panel.border = element_rect(colour = "black", fill=NA, size=0.05)
            # ) +
            
            
            
            EnvStats::stat_n_text(size = 4, y.pos = max(df_summary1$value, na.rm=T)*1.1 , y.expand.factor=0, 
                                  angle = 0, hjust = .5, family = "mono", fontface = "plain") + #295 bold
                
                theme(panel.background=element_blank(),
                      # axis.text.y=element_blank(),
                      # axis.ticks.y=element_blank(),
                      # https://stackoverflow.com/questions/46482846/ggplot2-x-axis-extreme-right-tick-label-clipped-after-insetting-legend
                      # stop axis being clipped
                      plot.title=element_text(), plot.margin = unit(c(5.5,12,5.5,5.5), "pt"),
                      legend.text=element_text(size=12),
                      legend.title=element_text(size=14),
                      legend.position="none",
                      axis.text.x  = element_text(size=10),
                      axis.text.y  = element_text(size=10),
                      axis.line.x = element_line(color="black"),
                      axis.line.y = element_line(color="black"),
                      plot.caption=element_text(hjust = 0, size = 7)) +
                           scale_color_hue(direction = v3, h.start=v4)
            
            
            
            
            
            print(pr1 + labs(y="Response", x = "Visit") + 
                      ggtitle(paste0("Individual responses ",
                                     length(unique(df$ID))," patients & arithmetic mean with 95% CI shown in black\nNumber of patient values at each time point") )
            )
        }
            
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if (input$plot == "plot1x") {   #untransformed means, transformed data
                
             #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # step 1 calcualte means and CIs
            df_summary <- df %>% # the names of the new data frame and the data frame to be summarised
                group_by(VISIT, variable) %>%                # the grouping variable
                summarise(mean_PL = mean(value, na.rm=TRUE),  # calculates the mean of each group
                          sd_PL = sd(value, na.rm=TRUE),      # calculates the sd of each group
                          n_PL = length(na.omit(value)),      # calculates the sample size per group
                          SE_PL = sd(value, na.rm=TRUE)/sqrt(length(na.omit(value)))) # SE of each group
            
            df_summary1 <- merge(df, df_summary)  # merge stats to dataset
            
            df_summary1$L2SE <- df_summary1$mean_PL - 2*df_summary1$SE_PL
            df_summary1$H2SE <- df_summary1$mean_PL + 2*df_summary1$SE_PL
            
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # step 2 log the raw data
            df_summary1$lvalue <- log(df_summary1$value)  #loG THE DATA
            
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # step 3 log my summary stats
            df_summary1$mean_PL<- log(df_summary1$mean_PL)
            df_summary1$L2SE<- log(df_summary1$L2SE)
            df_summary1$H2SE <- log(df_summary1$H2SE)
            #df2 <- unique(df_summary1[,c("VISIT","mean_PL","L2SE", "H2SE")])
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # plot the raw data
            
            pr1 <- ggplot((df_summary1), aes(x = VISIT, y =lvalue, color = ID)) +
                geom_line(size=v1, alpha=v2) +
                geom_point( data=df_summary1, aes(x=VISIT , y=mean_PL), colour="red")           +
                geom_line( data=df_summary1, aes(x=VISIT , y=mean_PL), colour="red") +
                geom_errorbar(data=df_summary1, 
                              aes( x=VISIT , y=mean_PL,ymin=L2SE, ymax=H2SE ), color = "black",
                              width=0.05, lwd = 0.05) +
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                # create ticks, replace log values with antlogs  
                
                scale_y_continuous(
                    breaks= 
                        log(c(0.01,      0.1,      1 ,      10 ,      100) )  ,  # this is where the values go
                    labels= c(0.01,      0.1 ,     1,       10 ,      100)) +      
                
                scale_x_continuous(breaks = c(unique(df_summary1$VISIT)),
                                   labels = 
                                       c(unique(df_summary1$VISIT))
                ) +
                
               scale_color_hue(direction = v3, h.start=v4) +
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                # create visit counts log the y position 
                
                EnvStats::stat_n_text(size = 4, y.pos = log(max(df_summary1$value, na.rm=T)*1.1) , y.expand.factor=0, 
                                      angle = 0, hjust = .5, family = "mono", fontface = "plain") + #295 bold
                
                theme(panel.background=element_blank(),
                      # axis.text.y=element_blank(),
                      # axis.ticks.y=element_blank(),
                      # https://stackoverflow.com/questions/46482846/ggplot2-x-axis-extreme-right-tick-label-clipped-after-insetting-legend
                      # stop axis being clipped
                      plot.title=element_text(), plot.margin = unit(c(5.5,12,5.5,5.5), "pt"),
                      plot.background=element_rect(fill = input$Design),   
                      legend.text=element_text(size=12),
                      legend.title=element_text(size=14),
                      legend.position="none",
                      axis.text.x  = element_text(size=10),
                      axis.text.y  = element_text(size=10),
                      axis.line.x = element_line(color="black"),
                      axis.line.y = element_line(color="black"),
                      plot.caption=element_text(hjust = 0, size = 7))
            
            
            
            print(pr1 + labs(y="Response", x = "Visit") + 
                      ggtitle(paste0("Individual responses ",
                                     length(unique(df$ID))," patients & arithmetic mean with 95% CI shown in black\nNumber of patient values at each time point") )
            )
            
            
            
            
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
            
            
        } else if (input$plot == "plot1a") {    #MEDIANS
            
            
            
            df_summary <- df %>% # the names of the new data frame and the data frame to be summarised
                group_by(VISIT, variable) %>%                # the grouping variable
                summarise(mean_PL = mean(value, na.rm=TRUE),  # calculates the mean of each group
                          sd_PL = sd(value, na.rm=TRUE),      # calculates the sd of each group
                          n_PL = length(na.omit(value)),      # calculates the sample size per group
                          SE_PL = sd(value, na.rm=TRUE)/sqrt(length(na.omit(value))) ,# SE of each group
                          median = median(value, na.rm=TRUE) ,
                          medianL =  sort(value)[qbinom(c(.025), size=length(value), prob=.5)]   ,
                          medianU =  sort(value)[qbinom(c(.975), size=length(value), prob=.5)]   )
            
            df_summary1 <- merge(df, df_summary)  # merge stats to dataset
            
            
            
            
            df_summary1$L2SE <- df_summary1$mean_PL - 2*df_summary1$SE_PL
            df_summary1$H2SE <- df_summary1$mean_PL + 2*df_summary1$SE_PL
            
            df_summary1$L2SE <- df_summary1$medianL
            df_summary1$H2SE <- df_summary1$medianU
            
            pr1 <- ggplot((df_summary1), aes(x = VISIT, y =value, color = ID)) +
                geom_line( size=v1, alpha=v2) +
                #scale_color_gradient(low = "blue", high = "red")+
                #scale_color_brewer(palette = "Dark2") +
                stat_summary(geom="line",  fun=median, colour="black", lwd=0.5) +  # , linetype="dashed"
                stat_summary(geom="point", fun=median, colour="black") +
                geom_errorbar(data=(df_summary1), 
                              aes( ymin=L2SE, ymax=H2SE ), color = "black",
                              width=0.05, lwd = 0.05) +
                scale_y_continuous(expand = c(.1,0) ) +
                
                scale_x_continuous(breaks = c(unique(df$VISIT)),
                                   labels = 
                                       c(unique(df$VISIT))
                ) +
                
                
                
                
              #  geom_segment(aes(x = 1, xend = 3, y = (1), yend=(1)), color = "blue" , size=0.05, linetype="dashed", alpha=0.02) +
               # geom_segment(aes(x = 1, xend = 3, y = (2.5), yend=(2.5)), color = "blue" , size=0.05, linetype="dashed", alpha=0.02) +
                
                theme(
                    # get rid of panel grids
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    # Change plot and panel background
                    plot.background=element_rect(fill =  input$Design),
                    panel.background = element_rect(fill = 'black'),
                    # Change legend
                    legend.position = c(0.6, 0.07),
                    legend.direction = "horizontal",
                    legend.background = element_rect(fill = "black", color = NA),
                    legend.key = element_rect(color = "gray", fill = "black"),
                    legend.title = element_text(color = "white"),
                    legend.text = element_text(color = "white")
                ) +
                
                
               scale_color_hue(direction = v3, h.start=v4) +
                
                EnvStats::stat_n_text(size = 4, y.pos = max(df_summary1$value, na.rm=T)*1.1 , y.expand.factor=0, 
                                      angle = 0, hjust = .5, family = "mono", fontface = "plain") + #295 bold
                
                theme(panel.background=element_blank(),
                      # axis.text.y=element_blank(),
                      # axis.ticks.y=element_blank(),
                      # https://stackoverflow.com/questions/46482846/ggplot2-x-axis-extreme-right-tick-label-clipped-after-insetting-legend
                      # stop axis being clipped
                      plot.title=element_text(), plot.margin = unit(c(5.5,12,5.5,5.5), "pt"),
                      legend.text=element_text(size=12),
                      legend.title=element_text(size=14),
                      legend.position="none",
                      axis.text.x  = element_text(size=10),
                      axis.text.y  = element_text(size=10),
                      axis.line.x = element_line(color="black"),
                      axis.line.y = element_line(color="black"),
                      plot.caption=element_text(hjust = 0, size = 7))
            
            
            
            
            
            print(pr1 + labs(y="Response", x = "Visit") + 
                      ggtitle(paste0("Individual responses ",
                                     length(unique(df$ID))," patients & medians with 95% CI shown in black\nNumber of patient values at each time point") )
            )
            
            
        
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
        } else if (input$plot == "plot1ax") {    #MEDIANS 2
            
            
            
            df_summary <- df %>% # the names of the new data frame and the data frame to be summarised
                group_by(VISIT, variable) %>%                # the grouping variable
                summarise(mean_PL = mean(value, na.rm=TRUE),  # calculates the mean of each group
                          sd_PL = sd(value, na.rm=TRUE),      # calculates the sd of each group
                          n_PL = length(na.omit(value)),      # calculates the sample size per group
                          SE_PL = sd(value, na.rm=TRUE)/sqrt(length(na.omit(value))) ,# SE of each group
                          mean_PL = median(value, na.rm=TRUE) ,
                          medianL =  sort(value)[qbinom(c(.025), size=length(value), prob=.5)]   ,
                          medianU =  sort(value)[qbinom(c(.975), size=length(value), prob=.5)]   )
            
            df_summary1 <- merge(df, df_summary)  # merge stats to dataset
            
            
            
            
          #  df_summary1$L2SE <- df_summary1$mean_PL - 2*df_summary1$SE_PL
          #  df_summary1$H2SE <- df_summary1$mean_PL + 2*df_summary1$SE_PL
            
            df_summary1$L2SE <- df_summary1$medianL
            df_summary1$H2SE <- df_summary1$medianU
            
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # step 2 log the raw data
            df_summary1$lvalue <- log(df_summary1$value)  #loG THE DATA
            
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # step 3 log my summary stats
            df_summary1$mean_PL<- log(df_summary1$mean_PL)
            df_summary1$L2SE<- log(df_summary1$L2SE)
            df_summary1$H2SE <- log(df_summary1$H2SE)
            #df2 <- unique(df_summary1[,c("VISIT","mean_PL","L2SE", "H2SE")])
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # plot the raw data
            
            pr1 <- ggplot((df_summary1), aes(x = VISIT, y =lvalue, color = ID)) +
                geom_line( size=v1, alpha=v2) +
                geom_point( data=df_summary1, aes(x=VISIT , y=mean_PL), colour="red")           +
                geom_line( data=df_summary1, aes(x=VISIT , y=mean_PL), colour="red") +
                geom_errorbar(data=df_summary1, 
                              aes( x=VISIT , y=mean_PL,ymin=L2SE, ymax=H2SE ), color = "black",
                              width=0.05, lwd = 0.05) +
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                # create ticks, replace log values with antlogs  
                
                scale_y_continuous(
                    breaks= 
                        log(c(0.01,      0.1,      1 ,      10 ,      100) )  ,  # this is where the values go
                    labels= c(0.01,      0.1 ,     1,       10 ,      100)) +     
                
                scale_x_continuous(breaks = c(unique(df_summary1$VISIT)),
                                   labels = 
                                       c(unique(df_summary1$VISIT))
                ) +
               scale_color_hue(direction = v3, h.start=v4) +
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                # create visit counts log the y position 
                
                EnvStats::stat_n_text(size = 4, y.pos = log(max(df_summary1$value, na.rm=T)*1.1) , y.expand.factor=0, 
                                      angle = 0, hjust = .5, family = "mono", fontface = "plain") + #295 bold
                
                theme(panel.background=element_blank(),
                      # axis.text.y=element_blank(),
                      # axis.ticks.y=element_blank(),
                      # https://stackoverflow.com/questions/46482846/ggplot2-x-axis-extreme-right-tick-label-clipped-after-insetting-legend
                      # stop axis being clipped
                      plot.background=element_rect(fill = input$Design),    
                      plot.title=element_text(), plot.margin = unit(c(5.5,12,5.5,5.5), "pt"),
                      legend.text=element_text(size=12),
                      legend.title=element_text(size=14),
                      legend.position="none",
                      axis.text.x  = element_text(size=10),
                      axis.text.y  = element_text(size=10),
                      axis.line.x = element_line(color="black"),
                      axis.line.y = element_line(color="black"),
                      plot.caption=element_text(hjust = 0, size = 7))
            
            
            
            print(pr1 + labs(y="Response", x = "Visit") + 
                      ggtitle(paste0("Individual responses ",
                                     length(unique(df$ID))," patients & medians with 95% CI shown in black\nNumber of patient values at each time point") )
            )
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
            
            
        }   else if (input$plot == "plot2") {   # MEANS ON LOG THEN EXP
            
             
              
                df$lvalue <- log(df$value)  #log the values
                  

            est <- df_summary <- df %>% # the names of the new data frame and the data frame to be summarised
                group_by(VISIT, variable) %>%                # the grouping variable
                summarise(mean_PL = mean(lvalue, na.rm=TRUE),  # calculates the mean of each group
                          sd_PL = sd(lvalue, na.rm=TRUE),      # calculates the sd of each group
                          n_PL = length(na.omit(lvalue)),      # calculates the sample size per group
                          SE_PL = sd(lvalue, na.rm=TRUE)/sqrt(length(na.omit(lvalue)))) # SE of each group
            
            
            df_summary1 <- merge(df, df_summary)  # merge stats to dataset
            
            # dont want error bars to exceed plausible range
            df_summary1$L2SE <- df_summary1$mean_PL - 2*df_summary1$SE_PL
            df_summary1$H2SE <- df_summary1$mean_PL + 2*df_summary1$SE_PL
            
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
            pr1 <- ggplot((df_summary1), aes(x = VISIT, y =lvalue, color = ID)) +
                geom_line( size=v1, alpha=v2) +
                stat_summary(geom="line",  fun=mean, colour="black", lwd=0.5) +  
                stat_summary(geom="point", fun=mean, colour="black") +
                geom_errorbar(data=(df_summary1), 
                              aes( ymin=L2SE, ymax=H2SE ), color = "black",
                              width=0.05, lwd = 0.05) +
                
                scale_x_continuous(breaks = c(unique(df$VISIT)),
                                   labels = c(unique(df$VISIT))
                ) +
                
                scale_y_continuous(
                    breaks= 
                        c(    log(0.01), log(.1),  log(1) , log(10) , log(100) ) ,  # this is where the values go
                        labels= c(0.01,      0.1 ,     1,       10 ,      100)) +      
                
                # theme(axis.text.y   = element_text(size=10),
                #       axis.text.x   = element_text(size=10),
                #       axis.title.y  = element_text(size=14),
                #       axis.title.x  = element_text(size=14),
                #       panel.background = element_blank(),
                #       panel.grid.major = element_blank(), 
                #       panel.grid.minor = element_blank(),
                #       legend.position="none",
                #       legend.text=NULL,
                #       legend.title=NULL,
            #       axis.line = element_line(colour = "black", size=0.05),
            #       panel.border = element_rect(colour = "black", fill=NA, size=0.05)
            # ) +
           scale_color_hue(direction = v3, h.start=v4) +
            EnvStats::stat_n_text(size = 4, y.pos = max(df_summary1$lvalue, na.rm=T)*1.1 , y.expand.factor=0, 
                                  angle = 0, hjust = .5, family = "mono", fontface = "plain") +#295 bold
                
                theme(panel.background=element_blank(),
                      # axis.text.y=element_blank(),
                      # axis.ticks.y=element_blank(),
                      # https://stackoverflow.com/questions/46482846/ggplot2-x-axis-extreme-right-tick-label-clipped-after-insetting-legend
                      # stop axis being clipped
                      plot.background=element_rect(fill = input$Design),    
                      plot.title=element_text(), plot.margin = unit(c(5.5,12,5.5,5.5), "pt"),
                      legend.text=element_text(size=12),
                      legend.title=element_text(size=14),
                      legend.position="none",
                      axis.text.x  = element_text(size=10),
                      axis.text.y  = element_text(size=10),
                      axis.line.x = element_line(color="black"),
                      axis.line.y = element_line(color="black"),
                      plot.caption=element_text(hjust = 0, size = 7))
            
            
            print(pr1 + labs(y="Response", x = "Visit") + 
                      ggtitle(paste0("Individual responses ",length(unique(df$ID))," patients & mean response ± 2 SE shown in black\nNumber of patient values at each time point") )
            )
            
            est$L2SE <- est$mean_PL - 2*est$SE_PL
            est$H2SE <- est$mean_PL + 2*est$SE_PL
            
            
            
            

            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
            
        }  else if (input$plot == "plot3") {   # GLS ON LOG THEN EXP
            
            df$x <- df$VISIT
            df$y <- (df$value)
            df$g <- df$ID
            
            library(scales)
            
            dplot <- merge(est,  df , by.x="Visit", by.y="VISIT")
            
            dplot <- dplot[complete.cases(dplot),]
            
            
            p <- ggplot(data = dplot, aes(x=Visit , y=Estimate, group=1)) 
                
            
            p1 <-  p +geom_line(data=df, aes(x = x, y = y, color = g, group=g), size=v1, alpha=v2)  +
                theme(text = element_text(size=10),
                      axis.text.x = element_text(angle = 0, vjust = 1, hjust=.5)) +
                geom_point() + geom_line() +   guides(color=FALSE) +
                geom_errorbar(aes(ymin=Lower, ymax=Upper ), color="black", width=.05, lwd=.2) +
                
                
                scale_y_continuous(breaks=c(0.01,      0.1 ,     1,       10 ,      100), trans='log', labels = comma) +
                
                
                ylab("Response")  +
                xlab("Visit")  +
                
                
                scale_x_continuous(breaks = c(unique(df$VISIT)),
                                   c(unique(df$VISIT))) +   # labels = comma) + #
                
               scale_color_hue(direction = v3, h.start=v4) +
                EnvStats::stat_n_text(size = 4, y.pos = max(log(dplot$value), na.rm=T)*1.1 ,
                                      y.expand.factor=0,  angle = 0, hjust = .5, family = "mono", fontface = "plain") +#295 bold
                
                
                theme(panel.background=element_blank(),
                      # axis.text.y=element_blank(),
                      # axis.ticks.y=element_blank(),
                      # https://stackoverflow.com/questions/46482846/ggplot2-x-axis-extreme-right-tick-label-clipped-after-insetting-legend
                      # stop axis being clipped
                      plot.background=element_rect(fill = input$Design),    
                      plot.title=element_text(), plot.margin = unit(c(5.5,12,5.5,5.5), "pt"),
                      legend.text=element_text(size=12),
                      legend.title=element_text(size=14),
                      axis.text.x  = element_text(size=10),
                      axis.text.y  = element_text(size=10),
                      axis.line.x = element_line(color="black"),
                      axis.line.y = element_line(color="black"),
                      plot.caption=element_text(hjust = 0, size = 7))
            
            print(p1 + labs(y="Response", x = "Visit") +
                      ggtitle(paste0("Individual responses ",
                                     length(unique(dplot$ID))," patients & modelled mean response with 95% CI from GLS model shown in black\nNumber of patient values at each time point") )
            )
            
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
            
            
        }  else if (input$plot == "plot3a") {   # LMM ON LOG THEN EXP
            
            df$x <- df$VISIT
            df$y <- (df$value)
            df$g <- df$ID
            
            library(scales)
            
            dplot <- merge(lmm.est,  df , by.x="Visit", by.y="VISIT")
            
            dplot <- dplot[complete.cases(dplot),]
            
            
            p <- ggplot(data = dplot, aes(x=Visit , y=Estimate, group=1)) +
                geom_line(data=df, aes(x = x, y = y, color = g, group=g), size=v1, alpha=v2)  +
                geom_point() + geom_line() +   guides(color=FALSE) +
                geom_errorbar(aes(ymin=Lower, ymax=Upper ), color="black", width=.05, lwd=.2) 
            
            p1 <-  p +
                theme(text = element_text(size=10),
                      axis.text.x = element_text(angle = 0, vjust = 1, hjust=.5)) +
                
                scale_y_continuous(breaks=c(0.01,      0.1 ,     1,       10 ,      100), trans='log', labels = comma) +
                
                
                ylab("Response")  +
                xlab("Visit")  +
                
                
                scale_x_continuous(breaks = c(unique(df$VISIT)),
                                   c(unique(df$VISIT))) +   # labels = comma) + #
                
               scale_color_hue(direction = v3, h.start=v4) +
            
                EnvStats::stat_n_text(size = 4, y.pos = max(log(dplot$value), na.rm=T)*1.1 ,
                                      y.expand.factor=0,  angle = 0, hjust = .5, family = "mono", fontface = "plain") +#295 bold
                
                
                theme(panel.background=element_blank(),
                      # axis.text.y=element_blank(),
                      # axis.ticks.y=element_blank(),
                      # https://stackoverflow.com/questions/46482846/ggplot2-x-axis-extreme-right-tick-label-clipped-after-insetting-legend
                      # stop axis being clipped
                      plot.background=element_rect(fill = input$Design),    
                      plot.title=element_text(), plot.margin = unit(c(5.5,12,5.5,5.5), "pt"),
                      legend.text=element_text(size=12),
                      legend.title=element_text(size=14),
                      axis.text.x  = element_text(size=10),
                      axis.text.y  = element_text(size=10),
                      axis.line.x = element_line(color="black"),
                      axis.line.y = element_line(color="black"),
                      plot.caption=element_text(hjust = 0, size = 7))
            
            print(p1 + labs(y="Response", x = "Visit") +
                      ggtitle(paste0("Individual responses ",
                                     length(unique(dplot$ID))," patients & modelled mean response with 95% CI from LMM shown in black\nNumber of patient values at each time point") )
            )
            
            
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
            
        } else if (input$plot == "plot4") {  # COMPARISON PLOT
            
            #_______________________________________________________________________________________
            # nop transformation
            
            df_summary <- df %>% # the names of the new data frame and the data frame to be summarised
                group_by(VISIT, variable) %>%                # the grouping variable
                summarise(Estimate = mean(value, na.rm=TRUE),  # calculates the mean of each group
                          SE = sd(value, na.rm=TRUE)/sqrt(length(na.omit(value)))) # SE of each group
            
            df_summary$Lower <- df_summary$Estimate - 2*df_summary$SE 
            df_summary$Upper <- df_summary$Estimate + 2*df_summary$SE 
            
            notran <- df_summary
            namz <- c("Visit", "Variable","Estimate", "se", "Lower", "Upper")
            
            names(notran) <- namz
            
            notran$Variable <- NULL
            #_______________________________________________________________________________________
            # medians
            
            df_summary <- df %>% # the names of the new data frame and the data frame to be summarised
                group_by(VISIT, variable) %>%                # the grouping variable
                summarise(Estimate = median(value, na.rm=TRUE),  # calculates the mean of each group
                          #sd_PL = sd(value, na.rm=TRUE),      # calculates the sd of each group
                          #n_PL = length(na.omit(value)),      # calculates the sample size per group
                          SE_PL = sd(value, na.rm=TRUE)/sqrt(length(na.omit(value))) ,# SE of each group
                          #median = median(value, na.rm=TRUE) ,
                          medianL =  sort(value)[qbinom(c(.025), size=length(value), prob=.5)]   ,
                          medianU =  sort(value)[qbinom(c(.975), size=length(value), prob=.5)]   )
            
            meds <- df_summary
            namz <- c("Visit", "Variable","Estimate", "se", "Lower", "Upper")
            
            names(meds) <- namz
            
            meds$Variable <- NULL
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # LOG , calc mean then exp
            df$lvalue <- log(df$value)  #log the GH values!!!!!!!!!!!!!!!
            
            
            ##log
            df_summary <- df %>% # the names of the new data frame and the data frame to be summarised
                group_by(VISIT, variable) %>%                # the grouping variable
                summarise(Estimate = mean(lvalue, na.rm=TRUE),  # calculates the mean of each group
                          SE = sd(lvalue, na.rm=TRUE)/sqrt(length(na.omit(lvalue)))) # SE of each group
            
            df_summary$Lower <- df_summary$Estimate - 2*df_summary$SE 
            df_summary$Upper <- df_summary$Estimate + 2*df_summary$SE 
            df_summary[c(3:6)] <- lapply(df_summary[c(3:6)], function(x) exp(x))
            tran <- df_summary
            namz <- c("Visit", "Variable","Estimate", "se", "Lower", "Upper")
            
            names(tran) <- namz
            tran$Variable <- NULL
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            tran<- as.data.frame(tran)
            notran<- as.data.frame(notran)
            est<- as.data.frame(est)
            meds <- as.data.frame(meds)
            
            notran$Approach <- "Untransformed means"
            meds$Approach <- "Untransformed medians"
            tran$Approach <- "Log transformed means"
            est$Approach <- "GLS model estimates"
            lmm.est$Approach <- "LMM model estimates"
            dd <- rbind( notran,meds, tran, est, lmm.est)
            
            
            
            
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # plot
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
            # The errorbars overlapped, so use position_dodge to move them horizontally
            pd <- position_dodge(0.1) # move them .05 to the left and right
            
            pd1 <- ggplot(dd, aes(x= Visit, y=Estimate, colour=Approach)) + 
                geom_errorbar(aes(ymin=Lower, ymax=Upper), width=.1, position=pd) +
                geom_line(position=pd) +
                geom_point(position=pd) +
                scale_x_continuous(breaks = c(unique(df$VISIT)),
                                   labels = c(unique(df$VISIT))
                ) +
                
                theme(panel.background=element_blank(),
                      # axis.text.y=element_blank(),
                      # axis.ticks.y=element_blank(),
                      # https://stackoverflow.com/questions/46482846/ggplot2-x-axis-extreme-right-tick-label-clipped-after-insetting-legend
                      # stop axis being clipped
                      plot.background=element_rect(fill = input$Design),    
                      plot.title=element_text(), plot.margin = unit(c(5.5,12,5.5,5.5), "pt"),
                      legend.text=element_text(size=12),
                      legend.title=element_text(size=14),
                      axis.text.x  = element_text(size=10),
                      axis.text.y  = element_text(size=10),
                      axis.line.x = element_line(color="black"),
                      axis.line.y = element_line(color="black"),
                      plot.caption=element_text(hjust = 0, size = 7))
            
            print(pd1 + labs(y="Response", x = "Visit") + 
                      ggtitle(paste0("Estimates of central tendancy from the different approaches") )
            )
            
            
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
            
            
        }  else if (input$plot == "plot5") { # residual plots
            
            f <- fit.regression()
            
            fit <- f$fit.res
             
            data <- make.data()
            
            d2 <- data$d
            
            d2$resid <- r <- resid(fit)
            
            d2$fitted <- fitted(fit)
            
            yl <- ylab('Residuals')
            
            xl <- xlab("time")
            
            p1 <- ggplot(d2 , aes(x=fitted , y=resid)) + geom_point (   colour="#69b3a2") + yl
            
            p3 <- ggplot(d2 , aes(x=VISIT , y=resid )) +  geom_point ( colour="#69b3a2") + yl  + xl +
                stat_summary(fun.data ="mean_sdl", geom='smooth')
            
            p4 <- ggplot(d2 , aes(sample=resid )) + stat_qq(colour="#69b3a2") +
                geom_abline(intercept=mean(r), slope=sd(r)  ,  colour="black") +
                xlab('Residuals')   +
                ggtitle( " ")
            
     
            library(gridExtra)
            library(grid)
            df <- data.frame(Residuals = r)
            p5 <- ggplot(df, aes(x = Residuals)) +
                geom_histogram(aes(y =..density..),
                               #breaks = seq(-50, 50, by = 2),
                               colour = "black",
                               fill = "#69b3a2") +
                stat_function(fun = dnorm, args = list(mean = 0, sd = sigma(fit)  ))
            
            grid.arrange(p1,  p3, p4,p5, ncol=2,
                         
                         top = textGrob(paste0(input$Plot, " GLS model fit diagnostics"),gp=gpar(fontsize=20,font=3)))
            
            
        }
        
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   }) 
    
    
  
    #---------------------------------------------------------------------------
    # Show the summaries
    
    output$reg.summary5 <- renderPrint({
      
      summary <- data.summ()$summ
      
      #  return(list(summary))
      
    })
    
       output$reg.summary3 <- renderPrint({
        
        summary <- data.summ2()$summ
        
        
      #  return(list(summary))
        
    })
    
      output$reg.summary4 <- renderPrint({
        
        summary <- fit.regression()$fit.res1
        
        print(kable(summary, format="pandoc", row.names = FALSE, digits = c(3)))
        
        # return(list(summary))
        
    })     
    
    
     
    output$reg.summary6 <- renderPrint({
        
        summary <- fit.regression2()$fit.res2a
        
        summary$se <- NULL
        
        print(kable(summary, format="pandoc", row.names = FALSE, digits = c(3)))
        
        # return(list(summary))
        
    })     
    
    ###################################################################
    #repeats############################################################
    ###################################################################
    
    output$reg.summary5a   <- renderPrint({

      summary <- data.summ()$summ

      #  return(list(summary))
      print(kable(summary, format="pandoc", row.names = FALSE, digits = c(3)))


    })
    
    
    output$reg.summary3a   <- renderPrint({

      summary <- data.summ2()$summ

       #return(list(summary))
      print(kable(summary, format="pandoc", row.names = FALSE, digits = c(3)))

    })

    output$reg.summary4a  <- renderPrint({

      summary <- fit.regression()$fit.res1

      print(kable(summary, format="pandoc", row.names = FALSE, digits = c(3)))

      #return(list(summary))

    })

# 
 
    output$reg.summary6a  <- renderPrint({

      summary <- fit.regression2()$fit.res2a

      summary$se <- NULL

      print(kable(summary, format="pandoc", row.names = FALSE, digits = c(3)))

      # return(list(summary))

     })     
    ###################################################################
    #repeats############################################################
    ###################################################################
    
    
    
    
    
    
    
    
    
    
    output$reg.summary8 <- renderPrint({
        
        summary <- check()$f
        return(list(summary))
        
    })     
    output$reg.summary9 <- renderPrint({
        
        summary <- check()$f0
        return(list(summary))
        
    })     
    
    output$reg.summary10 <- renderPrint({
        
        summary <- check()$f1
        return(list(summary))
    })     
    
    
    
    
})



# Run the application 
shinyApp(ui = ui, server = server)