#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# r code, using the baseline as a covariate
# proving relevel 
# predict and emmeans are equivalent with corSymm not AR1symm though be careful!
# some nice plots along the way
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
    colz = c("lightblue", "blue",  "lightgreen", "darkgreen")
    
    # function to create longitudinal data
    is.even <- function(x){ x %% 2 == 0 }
    
    library(MASS)
    
    #~~~~~~~~~~~~~~
    # DATA CREATION
    #~~~~~~~~~~~~~~
    
    ### set number of individuals
    n <- 100
    beta0 <-  1.1174 # intercept
    beta1 <- -0.2859 # slope
    ar.val <- 0.9  # auto correlation
    sigma <- 0.7995   # residual
    tau0  <-  1.2748  # intercept sd
    tau1  <-  0.2276  # slope sd  
    tau01 <- -0.62 # intercept slope correlation
    
    ### maximum number of possible observations
    m <- 6
    p <- round(runif(n,2,m))
    #p <- rep(6,n) #########################################################NEW
    
    ### simulate observation moments (assume everybody has 1st obs)
    obs <- unlist(sapply(p, function(x) c(1, sort(sample(2:m, x-1, replace=FALSE)))))
    #obs<-rep(1:6,n) ########################################################NEW
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
    #xx<-unlist(sapply(p, function(x) arima.sim(model=list(ar=ar.val), n=x) * sqrt(1-ar.val^.2) * sigma))
    #dat$eij <- reshape::melt(xx)$value #############################################################NEW
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
    
    df  <- dat

#~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create wide data using 
#~~~~~~~~~~~~~~~~~~~~~~~~~~

   w1 <-   merge(df,setNames(subset(df, VISIT==1,select=c("ID", "value")),
                             c('ID', 'baseline')), by='ID')
   
   w <-  w1[!w1$VISIT %in% 1,]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CALC SUMMARY STATS ON LOG VALUES AND MERGE BACK TO DATA
# f SLIMDOWN DATA
# df1 MERGED DATA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
  require(dplyr)
  
  df$x <- df$VISIT
  df$y <- (df$value)   #not LOGGING
  df$g <- df$ID
  
  mus <- df %>%
    group_by(x) %>%
    summarise(mu=mean(y,na.rm=T),med=median(y,na.rm=T),
              sd=sd(y,na.rm=T), n=length(na.omit(y)), se=sd(y,na.rm=T)/sqrt(length(na.omit(y))),
              lo = mean(y,na.rm=T)-2*sd(y,na.rm=T)/sqrt(length(na.omit(y))),
              hi = mean(y,na.rm=T)+2*sd(y,na.rm=T)/sqrt(length(na.omit(y)) ))

 # mean( df[df$x==1 , "y"] )

  df$y <- log(df$value)   #LOGGING
  
  dfs<- df %>%
    group_by(x) %>%
    summarise(mu=mean(y,na.rm=T),med=median(y,na.rm=T), 
              sd=sd(y,na.rm=T), n=length(na.omit(y)), se=sd(y,na.rm=T)/sqrt(length(na.omit(y))),
              lo = mean(y,na.rm=T)-2*sd(y,na.rm=T)/sqrt(length(na.omit(y))),
              hi = mean(y,na.rm=T)+2*sd(y,na.rm=T)/sqrt(length(na.omit(y)) ))

  # merge to data , so we have the log data included
  df1 <- merge(df,dfs,  by = intersect("x", "x"))
  
  # slim it down, y is logged
  f <- df1[,c("x","y","ID")]
  
  # rename some variables
  f$g <- f$ID
  f$x <- as.factor(f$x)
  f <- plyr::arrange(f, g, x)
  
  g <- f$g
  g <- as.vector(g)
  f$id<-rep(seq_along( rle(g)$values ), times = rle(g)$lengths )
  
  # prep for harrell's functions
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RUN ANALYSIS ON LOG DATA WITH NO INTERCEPT AND EXPONENTIATE RESULTS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  fit.res <-  
    tryCatch(gls(y ~ x +0 ,
                 correlation=corSymm(form=~ as.numeric(x)|id), 
                 weights=varIdent(form=~1|x),
                 f, na.action=na.omit) , 
             error=function(e) e)
  
  
  i <-  intervals(fit.res , which='coef')   # see ?intervals.gls
  se <-  sqrt(diag(vcov(fit.res)))
  i <-  cbind(as.data.frame(i$coef), se)
  i <-  i[,c(2,4,1,3)]  
  
  
  estx <- as.data.frame(i)
  est <- exp(estx)            # exponentiate back
  est$x<-as.numeric(1:m)
  rownames(est) <- NULL
  est <-  est[,c(5,1,2,3,4)]  
  names(est) <- c('Visit','Estimate','se', 'Lower', 'Upper')
  
  # collect objects
  fit.res1 <- est
  fit.res <- fit.res
  foo <- f # save data
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MODEL ESTIMATES AFTER EXPONENTIATING
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  est

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LMM

# df <- d <- dat
# require(lme4)
# require(dplyr)
# 
# df$x <- df$VISIT
# df$y <- log(df$value)
# df$g <- df$ID
# 
# 
# dfs<- df %>%
#   group_by(x) %>%
#   summarise(mu=mean(y,na.rm=T), sd=sd(y,na.rm=T), n=length(na.omit(y)), se=sd(y,na.rm=T)/sqrt(length(na.omit(y))),
#             lo = mean(y,na.rm=T)-2*sd(y,na.rm=T)/sqrt(length(na.omit(y))),
#             hi = mean(y,na.rm=T)+2*sd(y,na.rm=T)/sqrt(length(na.omit(y)) ))
# 
# df1 <- merge(df,dfs,  by = intersect("x", "x"))
# 
# f <-df1[,c("x","y","ID")]
# 
# 
# f$g <- (f$ID)
# f$x <- as.factor(f$x)
# f<-plyr::arrange(f, g, x)
# 
# g <- f$g
# g <- as.vector(g)
# f$id<-rep(seq_along( rle(g)$values ), times = rle(g)$lengths )
# 
# ###############################
# 
# fit.res <- tryCatch(lmer(y ~    x + 0+ (1 + as.numeric(x) | id), data = f), 
#                     error=function(e) e)
# 
# i <-  confint(fit.res )
# ci <- i[grep("^x", rownames(i)), ]
# fc <- fixef(fit.res)
# 
# se <-  fixef(fit.res)  # just filling this , not needed
# i <-  cbind(as.data.frame(fc), se, ci)
# 
# estx <- as.data.frame(i)
# 
# est <- exp(estx)            # exponentiate back
# est$x<-as.numeric(1:m)
# rownames(est) <- NULL
# est <-  est[,c(5,1,2,3,4)]  
# names(est) <- c('Visit','Estimate','se', 'Lower', 'Upper')
# 
# fit.res2a=est 
# fit.res2b=fit.res ))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~
# CALL UP ORIGINAL DATA
#~~~~~~~~~~~~~~~~~~~~~~

  df1 <- dat  # refresh the data

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SUMMARY STATS ON ORIGINAL DATA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # quick check
  tapply(df1$value,df1$VISIT, mean)  
  
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SUMMARY STATS ON logged DATA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    df1 <- dat  # refresh the data
    # quick check
    tapply(log(df1$value),df1$VISIT, mean)
    
    df1$value <- log(df1$value)
    require(dplyr)
    
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
    
    food3 <- data.frame(df_summaryx)
    
    #food2[c(4,5,6,8,9,10,11,12)] <- lapply(food2[c(4,5,6,8,9,10,11,12)], function(x) exp(x))   # exponentiated here!
    
    names(food3) <- c("Visit","Variable","Data points","Mean", "SD","SE","Subjects", "lower95%CI", "Upper95%CI", "Med", "Min","Max", "Missing")  
    print(kable(food3, format="pandoc", row.names = FALSE, digits = c(3)))


    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # SUMMARY STATS ON logged DATA, then exponentiated
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    df1 <- dat  # refresh the data
    # quick check
    # quick check
    exp(tapply(log(df1$value),df1$VISIT, mean))
    
    df1$value <- log(df1$value)
    require(dplyr)
    
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
    
    food4 <- data.frame(df_summaryx)
    
    food4[c(4,5,6,8,9,10,11,12)] <- lapply(food4[c(4,5,6,8,9,10,11,12)], function(x) exp(x))
    
    names(food4) <- c("Visit","Variable","Data points","Mean", "SD","SE","Subjects", "lower95%CI", "Upper95%CI", "Med", "Min","Max", "Missing")  
    print(kable(food4, format="pandoc", row.names = FALSE, digits = c(3)))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CALL UP ORIGINAL DATA PLOT UNADULTERATED
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  df <-   dat  # refresh the data
  require(ggplot2)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # SUMMARY STATS ON ORIGINAL DATA
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
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
    geom_line( size=.5, alpha=0.2) +
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
      plot.background=element_rect(fill = "white"),
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
          plot.caption=element_text(hjust = 0, size = 7))
  
  
  
  
  
  print(pr1 + labs(y="Response", x = "Visit") + 
          ggtitle(paste0("Individual responses ",
                         length(unique(df$ID)),
    " patients & arithmetic mean with 95% CI shown in black\nNumber of patient values at each time point") )
  )


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# df <- d <- dat  #untransformed means, transformed data
# 
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # step 1 calcualte means and CIs
# df_summary <- df %>% # the names of the new data frame and the data frame to be summarised
#   group_by(VISIT, variable) %>%                # the grouping variable
#   summarise(mean_PL = mean(value, na.rm=TRUE),  # calculates the mean of each group
#             sd_PL = sd(value, na.rm=TRUE),      # calculates the sd of each group
#             n_PL = length(na.omit(value)),      # calculates the sample size per group
#             SE_PL = sd(value, na.rm=TRUE)/sqrt(length(na.omit(value)))) # SE of each group
# 
# df_summary1 <- merge(df, df_summary)  # merge stats to dataset
# 
# df_summary1$L2SE <- df_summary1$mean_PL - 2*df_summary1$SE_PL
# df_summary1$H2SE <- df_summary1$mean_PL + 2*df_summary1$SE_PL
# 
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # step 2 log the raw data
# df_summary1$lvalue <- log(df_summary1$value)  #loG THE DATA
# 
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # step 3 log my summary stats
# df_summary1$mean_PL<- log(df_summary1$mean_PL)
# df_summary1$L2SE<- log(df_summary1$L2SE)
# df_summary1$H2SE <- log(df_summary1$H2SE)
# #df2 <- unique(df_summary1[,c("VISIT","mean_PL","L2SE", "H2SE")])
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # plot the raw data
# 
# pr1 <- ggplot((df_summary1), aes(x = VISIT, y =lvalue, color = ID)) +
#   geom_line( size=.5, alpha=0.2) +
#   geom_point( data=df_summary1, aes(x=VISIT , y=mean_PL), colour="red")           +
#   geom_line( data=df_summary1, aes(x=VISIT , y=mean_PL), colour="red") +
#   geom_errorbar(data=df_summary1, 
#                 aes( x=VISIT , y=mean_PL,ymin=L2SE, ymax=H2SE ), color = "black",
#                 width=0.05, lwd = 0.05) +
#   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   # create ticks, replace log values with antlogs  
#   
#   scale_y_continuous(
#     breaks= 
#       log(c(0.01, .1,  1 , 10 ,100) )  ,  # this is where the values go
#     labels= c(0.01,      0.1 ,     1,       10 ,      100)) +      
#   
#   scale_x_continuous(breaks = c(unique(df_summary1$VISIT)),
#                      labels = 
#                        c(unique(df_summary1$VISIT))
#   ) +
#   
#   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   # create visit counts log the y position 
#   
#   EnvStats::stat_n_text(size = 4, y.pos = log(max(df_summary1$value, na.rm=T)*1.1) , y.expand.factor=0, 
#                         angle = 0, hjust = .5, family = "mono", fontface = "plain") + #295 bold
#   
#   theme(panel.background=element_blank(),
#         # axis.text.y=element_blank(),
#         # axis.ticks.y=element_blank(),
#         # https://stackoverflow.com/questions/46482846/ggplot2-x-axis-extreme-right-tick-label-clipped-after-insetting-legend
#         # stop axis being clipped
#         plot.title=element_text(), plot.margin = unit(c(5.5,12,5.5,5.5), "pt"),
#         legend.text=element_text(size=12),
#         legend.title=element_text(size=14),
#         legend.position="none",
#         axis.text.x  = element_text(size=10),
#         axis.text.y  = element_text(size=10),
#         axis.line.x = element_line(color="black"),
#         axis.line.y = element_line(color="black"),
#         plot.caption=element_text(hjust = 0, size = 7))
# 
# 
# 
# print(pr1 + labs(y="Response", x = "Visit") + 
#         ggtitle(paste0("Individual responses ",
#                        length(unique(df$ID))," patients & arithmetic mean with 95% CI shown in black\nNumber of patient values at each time point") )
# )

 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CALL UP ORIGINAL DATA PLOT on log scale CALL IN EXPONETIATED GLS MODEL ESTIMATES AND PLOT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  df <- dat  # refresh the data
   
  df$x <- df$VISIT
  df$y <- df$value
  df$g <- df$ID
  
  est <- fit.res1  # EXPONETIATED GLS MODEL ESTIMATES 
  
  library(scales)
  
  dplot <- merge(est,  df , by.x="Visit", by.y="VISIT")
  dplot <- dplot[complete.cases(dplot),]
  
  p <- ggplot(data = dplot, aes(x=Visit , y=Estimate, group=1)) +
    geom_point() + geom_line() +   guides(color=FALSE) +
    geom_errorbar(aes(ymin=Lower, ymax=Upper ), color="black", width=.05, lwd=.2) 
  
  p1 <-  p +geom_line(data=df, aes(x = x, y = y, color = g, group=g),size=0.5,alpha=.2)  +
    theme(text = element_text(size=10),
          axis.text.x = element_text(angle = 0, vjust = 1, hjust=.5)) +
    
    scale_y_continuous(breaks=c(0.01, 0.1 , 1, 10, 100), trans='log', labels = comma) +
    
    
    ylab("Response")  +
    xlab("Visit")  +
    
    scale_x_continuous(breaks = c(unique(df$VISIT)),
                      labels= c(unique(df$VISIT))) +   # labels = comma) + #
    
    EnvStats::stat_n_text(size = 4, y.pos = max(log(dplot$value), na.rm=T)*1.1 ,
                          y.expand.factor=0,  angle = 0, hjust = .5, family = "mono", fontface = "plain") +#295 bold
    
    theme(panel.background=element_blank(),
          # axis.text.y=element_blank(),
          # axis.ticks.y=element_blank(),
          # https://stackoverflow.com/questions/46482846/ggplot2-x-axis-extreme-right-tick-label-clipped-after-insetting-legend
          # stop axis being clipped
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  REMINDER OF SUM STATS ON RAW DATA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  df_summary 
  food2

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  CALL IN INITIAL GLS MODEL ON LOGED DATA AS REFERENCE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  summary(fit.res)  # check to earlier model so data is correct, ANALYSED ON LOG y

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  LOG Y AND CHECK MODEL IS THE SAME TO CONFIRM DATA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  df <- dat  # refresh the data
  df$x <- df$VISIT
  df$y <- df$value
  df$g <- df$ID
  df$y <- log(df$y)  
  df$x <- factor(df$x)

  summary(tryCatch(gls(y ~ x  +0 ,
               correlation=corAR1(form=~ as.numeric(x)|g), 
               weights=varIdent(form=~1|x),
               df, na.action=na.omit) , 
           error=function(e) e))
 
  fit.res # earlier model
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  CREATE LONG DATA STRAT TO ADJUST FOR BASELINE
#  THE RESPONSE IS LOGGED IN THIS DATA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 # VISIT IS REMOVED FROM THIS DATA
  df <- dat  # refresh the data
  df$value <- log(df$value)  
  w1 <-   merge(df,setNames(subset(df, VISIT==1,select=c("ID", "value")),
                            c('ID', 'baseline')), by='ID')
  
  w <- w1[!w1$VISIT %in% 1, ]

  tmp <- w
  d <- tmp
  
  d <- droplevels(d)
  
  d$y <- d$value
  d$x <- factor(d$VISIT)
  d$g <- d$ID
 
  mdata  <- d
  dd <- datadist(mdata)
  options(datadist='dd')
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  FIT MODEL WITH BASELINE AS COVARIATE, LOGED DATA HERE : mdata
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  fit.resb <-  
    tryCatch(gls(y ~ x  + baseline ,
                 correlation=corSymm(form=~ as.numeric(x)|g),  #corAR1
                 weights=varIdent(form=~1|x),
                 mdata, na.action=na.omit) , 
             error=function(e) e)
  
  
    summary(fit.resb)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #  FIT AGAIN WITH FRANK HARRELL GLS : mdata
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   (f1  <-  tryCatch(Gls(y ~ x  + baseline ,
                 correlation=corSymm(form=~ as.numeric(x)|g), 
                 weights=varIdent(form=~1|x),
                 mdata, na.action=na.omit) , 
             error=function(e) e))
  
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #  IDEA TO CHANGE DATA SO INTERCEPT IS MEANININGFUL IE FOR A INFORMATIVE VALUE , MEDIAN BASELINE
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # WARNING, the baseline variable will not be unique
    # so getting summary stats in it is not advisable
    #    
    #summary(mdata$baseline)
    #adjustment = median(mdata$baseline)  # WRONG
    
    df <- dat  # refresh the data
    df$value <- log(df$value)  
   # tapply(df$value,df$VISIT, summary)
    adjustment <- median( df[df$VISIT==1 , "value"] )
  #   adjustment = quantile(df[df$VISIT==1 , "value"], probs=.5)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #  USE PACKAGE TO GET MODEL MEANS, THIS PRODUCES UNLOGGED ESTIMATES 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
    require(emmeans)
    
    (means <- emmeans(f1, "x", at = list(baseline = adjustment), data=mdata ))        # try on frank harrell
    
   # (means1 <- emmeans(fit.res, "x", at = list(baseline = adjustment), data=mdata ))  # try on nlme, works but time consuming
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #  PREDICT ON FRANK HARRELL,WE CAN EXPONENTIATE IN THE FUNCTION !
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    (z <- Predict(f1, x , baseline=adjustment, fun=exp)) # we can exp automatically
    z$baseline <- NULL
    harrell <- z
  
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #  NOW WE HAVE harrell (baseline model predictions) and (est non baseline predictions)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #  PREDICT ON THE NLME GLS MODEL
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    new <- data.frame(x = as.factor(2:6)  , baseline=adjustment)
    # pr <- predict(fit.res , new , se.fit=TRUE) # no se outputted , so need following package
    require(AICcmodavg)
    (predictSE.gls (fit.resb, new, se.fit=T))
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
    ## no adjustment model
    dd <- datadist(foo); 
    options(datadist='dd')
    
    harrell0 <-  
      tryCatch(Gls(y ~ x  ,  # +0 does not work with Harrell's Gls
                   correlation=corSymm(form=~ as.numeric(x)|id),   # corAR1
                   weights=varIdent(form=~1|x),
                   foo, na.action=na.omit) , 
               error=function(e) e)
    
    (z0 <- Predict(harrell0,x,    fun=exp)) # we can exp automatically
    
    #x <- rms::con
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #  NOW WE HAVE harrell (baseline model predictions) and (est non baseline predictions)
    # z baseline adjustment using GLS (from f1)
    # est no baseline adjustment from gls, uses intervals.gls (from fit.res)
    # z0  no baseline adjustment from Gls.. CIS differ to est maybe z or t dist? (from harrell0)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    est; z0; z
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # use the long data for patient profiles
    # use est no baseline adjustment, analysis on log response, exp results
    # plot using log scale
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
    df <- dat # long data
  
    library(scales)
    
    dplot <- merge(est,  df , by.x="Visit", by.y="VISIT")  # could use z0 rather than est
    
    dplot <- dplot[complete.cases(dplot),]
    
    p <- ggplot(data = dplot, aes(x=Visit , y=Estimate, group=1)) +
      geom_point() + geom_line() +   guides(color=FALSE) +
      geom_errorbar(aes(ymin=Lower, ymax=Upper ), color="black", width=.05, lwd=.2) 
    
    # using df here
    p1 <- p +geom_line(data=df, aes(x = VISIT, y = value, color = ID, group=ID),size=0.5,alpha=.2)  +
      theme(text = element_text(size=10),
            axis.text.x = element_text(angle = 0, vjust = 1, hjust=.5)) +
      
      scale_y_continuous(breaks=c(0.01,  0.1 ,  1,  10 , 100), trans='log', labels = comma) +
      
    ylab("Response")  +
    xlab("Visit")  +
    
   #  geom_hline(yintercept=2.649206, linetype="dashed", color = "red") +
    
    scale_x_continuous(breaks = c(unique(df$VISIT)),
                                labels=c(unique(df$VISIT))) +   # labels = comma) + #
    
    # scale_x_continuous(breaks = as.character(as.numeric(unique(df$VISIT))),
    #                    as.character(as.numeric(unique(df$VISIT)))) +   # labels = comma) + #
    # 
    
    EnvStats::stat_n_text(size = 4, y.pos = max(log(dplot$value), na.rm=T)*1.1 ,
                          y.expand.factor=0,  angle = 0, hjust = .5, family = "mono", fontface = "plain") +#295 bold
    
    theme(panel.background=element_blank(),
          # axis.text.y=element_blank(),
          # axis.ticks.y=element_blank(),
          # https://stackoverflow.com/questions/46482846/ggplot2-x-axis-extreme-right-tick-label-clipped-after-insetting-legend
          # stop axis being clipped
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
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get all estimates together and plot
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  est1 <- est[,c(1,2,4,5)]
  names(est1) <- names(z0)   # not using this
  
  z0 <- as.data.frame(z0)
  
  z <- as.data.frame(z)
  
  emus <- dfs[,c(1,2,7,8)]
  names(emus) <- names(z0)
  emus <- cbind(emus[,1], apply(emus[,c(2,3,4)],2,exp))
  
  musx <- mus[,c(1,2,7,8)]
  names(musx) <- names(z0)
  
  musx$Approach <- "Untransformed means"
  emus$Approach <- "log data calc means & exp"
  est1$Approach <- "exp(log GLS with no baseline adj)"  # not using
  z0$Approach <- "exp(log GLS with no baseline adj)"
  z$Approach <- "exp(log GLS data with baseline adj)"
  
  
  dd <- rbind( musx, emus,  z, z0)  # visit
  
  dd$x <- (as.numeric(dd$x))
  
  # The errorbars overlapped, so use position_dodge to move them horizontally
  pd <- position_dodge(0.4) # move them .05 to the left and right
  
  pd1 <- ggplot(dd, aes(x= x, y=yhat, colour=Approach)) + 
    geom_errorbar(aes(ymin=lower, ymax=upper), width=.1, position=pd) +
    geom_line(position=pd) +
    geom_point(position=pd) +
    scale_x_continuous(breaks = c(unique(dd$x)),
                       labels = c(unique(dd$x))
    ) +
    ylab("Response")  +
    xlab("Visit")  +
    
    theme(panel.background=element_blank(),
          # axis.text.y=element_blank(),
          # axis.ticks.y=element_blank(),
          # https://stackoverflow.com/questions/46482846/ggplot2-x-axis-extreme-right-tick-label-clipped-after-insetting-legend
          # stop axis being clipped
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
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # add in raw data
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  pd1 <- ggplot(dd, aes(x= x, y=yhat, colour=Approach)) + 
    geom_errorbar(aes(ymin=lower, ymax=upper), width=.1, position=pd) +
    geom_line(position=pd) +
    geom_point(position=pd) +
    scale_x_continuous(breaks = c(unique(dd$x)),
                       labels = c(unique(dd$x))
    ) 
    
    p1 <-  pd1 + geom_line(data=df, aes(x = VISIT, y = value, color = 'black', group=ID),
                           size=0.5,alpha=.25)  +
    theme(text = element_text(size=10),
          axis.text.x = element_text(angle = 0, vjust = 1, hjust=.5)) +
    
    scale_y_continuous(breaks=c(0.01, 0.1 , 1, 10 , 100), trans='log', labels = comma) +
    
   # EnvStats::stat_n_text(size = 4, y.pos = max(log(df$value), na.rm=T)*1.1 ,
    #                        y.expand.factor=0,  angle = 0, hjust = .5, family = "mono",
     #                       fontface = "plain") + 
      
      ylab("Response")  +
      xlab("Visit")  +
      
  
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
  
  print(p1 + labs(y="Response", x = "Visit") + 
          ggtitle(paste0("Estimates of central tendancy from the different approaches") )
  )
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # use the long data for patient profiles
  # use baseline adjustment, analysis on log response, exp results
  # plot using log scale, showing only change from basleine
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  df   <- dat # long data
  
  library(scales)
  
  zz<- as.data.frame(z)

  dp <- merge(zz,  df , by.x="x", by.y="VISIT", all=TRUE)
  
  dp$VISIT <- as.numeric(as.character(dp$x))
  dp$VISIT <- ifelse(  is.na( dp$VISIT) , 1, dp$VISIT)
  dp$VISIT <- factor(dp$VISIT)
  
  dp$VISIT <- (as.numeric(dp$VISIT))
  
  p <- ggplot(data = dp, aes(x=VISIT , y=yhat, group=1)) +
    geom_point() + geom_line() +   guides(color=FALSE) +
    geom_errorbar(aes(ymin=lower, ymax=upper ), color="black", width=.05, lwd=.2) 
  
  p1 <-  p +geom_line(data=dp, aes(x = VISIT, y = value, color = ID, group=ID),size=0.5,alpha=.2)  +
    theme(text = element_text(size=10),
          axis.text.x = element_text(angle = 0, vjust = 1, hjust=.5)) +
   
    scale_y_continuous(breaks=c(0.01,  0.1 ,  1,  10 ,  100), trans='log', labels = comma) +
    
     # ylab("Response")  +
     # xlab("Visit")  +

    scale_x_continuous(breaks = unique(dp$VISIT), labels=unique(dp$VISIT)) +   # labels = comma) + #
    
    # scale_x_continuous(breaks = c(unique(df$VISIT)),
    #                    c(unique(df$VISIT))) +   # labels = comma) + #
    
    EnvStats::stat_n_text(size = 4, y.pos = max(log(dplot$value), na.rm=T)*1.1 ,
                          y.expand.factor=0,  angle = 0, hjust = .5, family = "mono", fontface = "plain") +#295 bold
    
    
    theme(panel.background=element_blank(),
          # axis.text.y=element_blank(),
          # axis.ticks.y=element_blank(),
          # https://stackoverflow.com/questions/46482846/ggplot2-x-axis-extreme-right-tick-label-clipped-after-insetting-legend
          # stop axis being clipped
          plot.title=element_text(), plot.margin = unit(c(5.5,12,5.5,5.5), "pt"),
          legend.text=element_text(size=12),
          legend.title=element_text(size=14),
          axis.text.x  = element_text(size=10),
          axis.text.y  = element_text(size=10),
          axis.line.x = element_line(color="black"),
          axis.line.y = element_line(color="black"),
          plot.caption=element_text(hjust = 0, size = 7),
          axis.title = element_text()
          ) +
  
    ylab("Response")  +
    xlab("Visit")  
  
  # print(p1 + #labs(y="Response", x = "Visit") +
  #         ggtitle(paste0("Individual responses ",
  #                        length(unique(dplot$ID)),
  # " patients & modelled mean response with 95% CI from GLS model shown in black\nNumber of patient values at each time point") )
  # )
  
  
   px <- p1 + labs(title = paste0("Individual responses ",
                            length(unique(dplot$ID)),
                            " patients & modelled mean response with 95% CI from GLS model shown in black\nNumber of patient values at each time point. Adjusting for baseline version of response."),
                #subtitle = "Plot of length by dose",
                caption = paste0("Adjusted for baseline (visit 1) , baseline= ",p2(adjustment),"")
             )
 
   print(px)
  
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# remember the intercept is when baseline is zero # so subtract the median baseline so the read out
# is the value of the response when the baseline is at median baseline 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  cp <- list(corSymm, corAR1, corCAR1, corExp, corLin, corGaus,corSpher) 
  k=1  # if you use say corAR1 the approach below wont work!!
  
# original baseline adjusted model on logged original data
fit.res <-  
  tryCatch(gls(y ~ x  + baseline ,
               correlation=cp[[k]](form=~ as.numeric(x)|g), 
               weights=varIdent(form=~1|x),
               mdata, na.action=na.omit) , 
           error=function(e) e)


summary(fit.res)
exp(fit.res$coefficients)
z # as baseline is not centered exp model intercept wont match this

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# here is a new variable
# lets adjust baseline so model intercept and predictions align
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#mdata$base <- mdata$baseline-median(mdata$baseline)
mdata$base <- mdata$baseline-adjustment
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # intercept is now better interpretation
  mdata$x <- relevel(mdata$x, ref= "2")   ##new
  
  fit.res <-  
    tryCatch(gls(y ~ x  + base ,
                 correlation=cp[[k]](form=~ as.numeric(x)|g), 
                 weights=varIdent(form=~1|x),
                 mdata, na.action=na.omit) , 
             error=function(e) e)
  
  
  summary(fit.res)
  exp(fit.res$coefficients)
  z # compare exp model intercept to prediction at visit 2 here (x var)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  mdata$x <- relevel(mdata$x, ref= "3")   ##new
  fit.res <- update(fit.res)    # run the same model
  summary(fit.res)
  exp(fit.res$coefficients)
  z # compare exp model intercept to prediction at visit 3 here (x var)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  mdata$x <- relevel(mdata$x, ref= "4")   ##new
  fit.res <- update(fit.res)    # run the same model
  summary(fit.res)
  exp(fit.res$coefficients)
  z # compare exp model intercept to prediction at visit 4 here (x var)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  mdata$x <- relevel(mdata$x, ref= "5")   ##new
  fit.res <- update(fit.res)    # run the same model
  summary(fit.res)
  exp(fit.res$coefficients)
  z # compare exp model intercept to prediction at visit 5 here (x var)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  mdata$x <- relevel(mdata$x, ref= "6")   ##new
  fit.res <- update(fit.res)    # run the same model
  summary(fit.res)
  exp(fit.res$coefficients)
  z # compare exp model intercept to prediction at visit 6 here (x var)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# check it out using Harrell's Gls
# here new variable.
  
  mdata$x <- relevel(mdata$x, ref= "4")   ##new
  dd <- datadist(mdata); 
  options(datadist='dd')
  
  (f  <-  tryCatch(Gls(y ~ x  + base ,
                      correlation=cp[[k]](form=~ as.numeric(x)|g), 
                      weights=varIdent(form=~1|x),
                      mdata, na.action=na.omit) , 
                  error=function(e) e))
  exp(f$coefficients)
  harrell  # Harrell or z
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
