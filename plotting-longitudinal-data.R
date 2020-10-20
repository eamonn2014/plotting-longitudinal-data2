# r code
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


 
    
    
    ### set number of individuals
    n <- 100
    beta0 <-  1.1174
    beta1 <- -0.2859 
    ar.val <- 0.9
    sigma <- 0.7995 
    tau0  <-  1.2748
    tau1  <-  0.2276
    tau01 <- -0.62
    
    ### maximum number of possible observations
    m <- 6
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
    
    df <- d <- dat
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  
  
   
    
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
    est$x<-as.numeric(1:m)
    rownames(est) <- NULL
    est <-  est[,c(5,1,2,3,4)]  
    names(est) <- c('Visit','Estimate','se', 'Lower', 'Upper')
    
    
    fit.res1=est
    fit.res=fit.res
    
 
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # LMM
   
    df <- d <- dat
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
    est$x<-as.numeric(1:m)
    rownames(est) <- NULL
    est <-  est[,c(5,1,2,3,4)]  
    names(est) <- c('Visit','Estimate','se', 'Lower', 'Upper')
    
    fit.res2a=est 
    fit.res2b=fit.res 
      
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
df1 <- d <- dat
    
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
   
  
    df1 <- d <- dat
    
    
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
 
     
  
 
    require(ggplot2)
    
    # Get the current regression data
    df1 <- d <- dat
    
    # get the model fit 
    
    est <- fit.res1
    
    lmm.est <-  fit.res2a
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
      
    df <- d <- dat
    
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
                             length(unique(df$ID))," patients & arithmetic mean with 95% CI shown in black\nNumber of patient values at each time point") )
      )
  
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      df <- d <- dat  #untransformed means, transformed data
      
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
        geom_line( size=.5, alpha=0.2) +
        geom_point( data=df_summary1, aes(x=VISIT , y=mean_PL), colour="red")           +
        geom_line( data=df_summary1, aes(x=VISIT , y=mean_PL), colour="red") +
        geom_errorbar(data=df_summary1, 
                      aes( x=VISIT , y=mean_PL,ymin=L2SE, ymax=H2SE ), color = "black",
                      width=0.05, lwd = 0.05) +
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # create ticks, replace log values with antlogs  
        
        scale_y_continuous(
          breaks= 
            log(c(0.01, .1,  1 , 10 ,100) )  ,  # this is where the values go
          labels= c(0.01,      0.1 ,     1,       10 ,      100)) +      
        
        scale_x_continuous(breaks = c(unique(df_summary1$VISIT)),
                           labels = 
                             c(unique(df_summary1$VISIT))
        ) +
        
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
      
      
      
      
      df <- d <- dat    #MEDIANS
      
      
      
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
        geom_line( size=.5, alpha=0.2) +
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
      
      
      df <- d <- dat    #MEDIANS 2
      
      
      
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
        geom_line( size=.5, alpha=0.2) +
        geom_point( data=df_summary1, aes(x=VISIT , y=mean_PL), colour="red")           +
        geom_line( data=df_summary1, aes(x=VISIT , y=mean_PL), colour="red") +
        geom_errorbar(data=df_summary1, 
                      aes( x=VISIT , y=mean_PL,ymin=L2SE, ymax=H2SE ), color = "black",
                      width=0.05, lwd = 0.05) +
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # create ticks, replace log values with antlogs  
        
        scale_y_continuous(
          breaks= 
            log(c(0.01, .1,  1 , 10 ,100) )  ,  # this is where the values go
          labels= c(0.01,      0.1 ,     1,       10 ,      100)) +      
        
        scale_x_continuous(breaks = c(unique(df_summary1$VISIT)),
                           labels = 
                             c(unique(df_summary1$VISIT))
        ) +
        
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
              legend.text=element_text(size=12),
              legend.title=element_text(size=14),
              legend.position="none",
              axis.text.x  = element_text(size=10),
              axis.text.y  = element_text(size=10),
              axis.line.x = element_line(color="black"),
              axis.line.y = element_line(color="black"),
              plot.caption=element_text(hjust = 0, size = 7))
      
      
      
      
      df <- d <- dat   # MEANS ON LOG THEN EXP
      
      
      
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
        geom_line( size=.5, alpha=0.2) +
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
      
      EnvStats::stat_n_text(size = 4, y.pos = max(df_summary1$lvalue, na.rm=T)*1.1 , y.expand.factor=0, 
                            angle = 0, hjust = .5, family = "mono", fontface = "plain") +#295 bold
        
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
              ggtitle(paste0("Individual responses ",length(unique(df$ID))," patients & mean response ± 2 SE shown in black\nNumber of patient values at each time point") )
      )
      
      est$L2SE <- est$mean_PL - 2*est$SE_PL
      est$H2SE <- est$mean_PL + 2*est$SE_PL
      
      
      
      
      print(pr1 + labs(y="Response", x = "Visit") + 
              ggtitle(paste0("Individual responses ",
                             length(unique(df$ID))," patients & medians with 95% CI shown in black\nNumber of patient values at each time point") )
      )
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      df <- d <- dat  # GLS ON LOG THEN EXP
      
      df$x <- df$VISIT
      df$y <- (df$value)
      df$g <- df$ID
      
      est <- fit.res1
      
      library(scales)
      
      dplot <- merge(est,  df , by.x="Visit", by.y="VISIT")
      
      dplot <- dplot[complete.cases(dplot),]
      
      
      p <- ggplot(data = dplot, aes(x=Visit , y=Estimate, group=1)) +
        geom_point() + geom_line() +   guides(color=FALSE) +
        geom_errorbar(aes(ymin=Lower, ymax=Upper ), color="black", width=.05, lwd=.2) 
      
      p1 <-  p +geom_line(data=df, aes(x = x, y = y, color = g, group=g),size=0.5,alpha=.2)  +
        theme(text = element_text(size=10),
              axis.text.x = element_text(angle = 0, vjust = 1, hjust=.5)) +
        
        scale_y_continuous(breaks=c(0.01,      0.1 ,     1,       10 ,      100), trans='log', labels = comma) +
        
        
        ylab("Response")  +
        xlab("Visit")  +
        
        
        scale_x_continuous(breaks = c(unique(df$VISIT)),
                           c(unique(df$VISIT))) +   # labels = comma) + #
        
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
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      df <- d <- dat   # LMM ON LOG THEN EXP
      
      df$x <- df$VISIT
      df$y <- (df$value)
      df$g <- df$ID
      
      library(scales)
      
      dplot <- merge(lmm.est,  df , by.x="Visit", by.y="VISIT")
      
      dplot <- dplot[complete.cases(dplot),]
      
      
      p <- ggplot(data = dplot, aes(x=Visit , y=Estimate, group=1)) +
        geom_point() + geom_line() +   guides(color=FALSE) +
        geom_errorbar(aes(ymin=Lower, ymax=Upper ), color="black", width=.05, lwd=.2) 
      
      p1 <-  p +geom_line(data=df, aes(x = x, y = y, color = g, group=g),size=0.5,alpha=.2)  +
        theme(text = element_text(size=10),
              axis.text.x = element_text(angle = 0, vjust = 1, hjust=.5)) +
        
        scale_y_continuous(breaks=c(0.01,      0.1 ,     1,       10 ,      100), trans='log', labels = comma) +
        
        
        ylab("Response")  +
        xlab("Visit")  +
        
        
        scale_x_continuous(breaks = c(unique(df$VISIT)),
                           c(unique(df$VISIT))) +   # labels = comma) + #
        
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
                             length(unique(dplot$ID))," patients & modelled mean response with 95% CI from LMM shown in black\nNumber of patient values at each time point") )
      )
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      df <- d <- dat   # COMPARISON PLOT
      
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
      
      
      
      # residual plots
       
      
      fit <-  fit.res
      
      d2 <- d <- dat 
      
     
      
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
                   
                   top = textGrob(paste0(  " GLS model fit diagnostics"),gp=gpar(fontsize=20,font=3)))
      
    