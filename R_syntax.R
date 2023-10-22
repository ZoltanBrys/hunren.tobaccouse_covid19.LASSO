# METADATA 
########################################################################################################
# Related manuscript  : Comparing the correlates of daily tobacco use 
#               and COVID-19 vaccination refusal among the Hungarian active population
# R-Code      : Zoltan Brys
# Closed      : -  
########################################################################################################


# PACKAGES AND ENVIRONMENT
########################################################################################################
#preparation
  rm(list = ls()) #deleting prev data
  if (!("stats" %in% (.packages()) )) stop("R Environment is not fully loaded!") #checking R

#packages
  library(Matrix)
  library(MASS) #for CIs
  library(fmsb) #for Nagelkerke
  library(glmnet) #for lasso.
  library(car) #for bootstrap.
  library(rcompanion) #for effect sizes
  library(margins) #for AME
  library(pROC) #for ROC-AUC 95% CI

#function for calculating Jaccard similarity measure
jaccardc <- function (x, y) 
{
  
  if (anyNA(x)) stop("NAs in x")
  if (anyNA(y)) stop("NAs in y")
  if (length(x) != length(x)) stop("Lenghts are non-equal!")
  
  c11 = sum(x == 1 & y == 1)
  c10 = sum(x == 1 & y == 0)
  c01 = sum(x == 0 & y == 1)
  
  return (c11 / (c11 + c10 + c01))
}
########################################################################################################


# LOADING AND CHECKING PREPARED INPUT DATA 
########################################################################################################
  C19 <- readRDS("C19.RDS")
#it contains the variables described in Table1 y1,y2 and x1, x2...x21.
# variable sorsz is the common key variable with the open KDK dataset

#check the input data
  if (!exists("C19"))             stop("Input data was not loaded!")
  if (class(C19) != "data.frame") stop("Input data format is not correct!")
########################################################################################################


# Characteristics of the sample (Table 1.)
########################################################################################################
  C19_sum_raw <- summary(C19)
  print(C19_sum_raw) 
########################################################################################################


# Everyday tobacco use and COVID-19 vaccination refusal (A1)
########################################################################################################
#bivariate
  y1y2sm <- glm(y2 ~ y1, family=binomial(link='logit') , data=C19) 
  summary(y1y2sm) 
  exp(coefficients(y1y2sm)[2])
  exp(confint(y1y2sm))
  NagelkerkeR2(y1y2sm)

#multivariate
  y1y2cm <- glm(y1 ~ y2 + x3 + x4, family=binomial(link='logit') , data=C19) 
  summary(y1y2cm)
  exp(coefficients(y1y2cm)[2])
  exp(confint(y1y2cm))
  NagelkerkeR2(y1y2cm)
#since 95% CI boostrap is stochastic, small difference are to be expected

#Jaccard
  jxy <- na.omit(subset(C19, select= c("y1", "y2")))
  jaccardc(jxy$y1, jxy$y2)
########################################################################################################


# Correlates of everyday tobacco use (A2)
########################################################################################################
#scaling the data
  A2C19 <- scale(na.omit(C19[ ,c(-1,-3)]))

#A2 -- first LASSO model for visualization
  glmmodA2 <- glmnet(as.matrix(A2C19[,-1]), 
                     as.matrix(A2C19[,1]), 
                     alpha = 1, 
                     family = "binomial",
                     intercept = FALSE)
  
#result objects for A2
  A2bet <- NULL #results betas
  A21se <- NULL #results lambda1se
  A2ROC <- NULL #results ROC

#1000 times repeated beta/l1se/ROC-AUC calculations of A2 with intercept set to 0
  for (c1 in 1:1000)
    {
    #finding lambda.1se
    cv_model1 <- cv.glmnet(as.matrix(A2C19[,-1]), 
                           as.matrix(A2C19[,1]), 
                           alpha = 1, 
                           family="binomial", 
                           kfold=5 ,    
                           type.measure = "auc",
                           intercept = FALSE) 

    #evaulating the model at lambda.1se
    glmmod2 <- glmnet(as.matrix(A2C19[,-1]), 
                      as.matrix(A2C19[,1]), 
                      alpha = 1, 
                      family = "binomial", 
                      lambda = cv_model1$lambda.1se,
                      intercept = FALSE)
    
    #saving results of each model
    A2bet <- rbind(A2bet, t(as.matrix(coef(glmmod2))))
    A21se <- rbind(A21se, cv_model1$lambda.1se)
    A2ROC <- rbind(A2ROC, cv_model1$cvm[which(cv_model1$lambda.1se==cv_model1$lambda)])
  }

  
  #Plot figure 2A
  png(filename = "Figure2A.png",
      width = 1280, height = 960, units = "px", 
      pointsize = 3, res=600,
      bg = "white"
  )
  
  plot(glmmodA2, xvar="lambda", label=TRUE, , ylab="", xlab="ln(λ)", 
       cex.lab=1, cex.axis=1, cex.main=1, cex.sub="1", ylim = c(-1, 1))
  abline(v = log(quantile(A21se, probs = 0.75)), 
         col=c("blue", "blue"), 
         lty=2, 
         lwd=1)
  
  dev.off()
########################################################################################################

 
# Correlates of  of COVID-19 vaccination refusal (A3)
########################################################################################################
#scaling the data
  A3C19 <- scale(na.omit(C19[ ,c(-1,-2)]))

#A3
  glmmodA3 <- glmnet(as.matrix(A3C19[,-1]), 
                     as.matrix(A3C19[,1]), 
                     alpha = 1, 
                     family = "binomial",
                     intercept = FALSE)

#results objects for A2
  A3bet <- NULL #betas
  A31se <- NULL #lambda 1 se
  A3ROC <- NULL #ROC values

#1000 times repeated Beta calcualtions
  for (c1 in 1:1000)
    {
      #finding lambda.1se
      cv_model1 <- cv.glmnet(as.matrix(A3C19[,-1]), 
                             as.matrix(A3C19[,1]), 
                             alpha = 1, 
                             family="binomial", 
                             kfold=5 ,    
                             type.measure = "auc",
                             intercept = FALSE) 
      
      #evaulating the model at lambda.1se
      glmmod2 <- glmnet(as.matrix(A3C19[,-1]), 
                        as.matrix(A3C19[,1]), 
                        alpha = 1, 
                        family = "binomial", 
                        lambda = cv_model1$lambda.1se,
                        intercept = FALSE)
      
      #saving tesults of c1 model
      A3bet <- rbind(A3bet, t(as.matrix(coef(glmmod2))))
      A31se <- rbind(A31se, cv_model1$lambda.1se)
      A3ROC <- rbind(A3ROC, cv_model1$cvm[which(cv_model1$lambda.1se==cv_model1$lambda)])
  }


#Plot figure 2B
  png(filename = "Figure2B.png",
       width = 1280, height = 960, units = "px", 
       pointsize = 3, res=600,
       bg = "white"
  )
  
  plot(glmmodA3, xvar="lambda", label=TRUE, , ylab="", xlab="ln(λ)", 
       cex.lab=1, cex.axis=1, cex.main=1, cex.sub="1", ylim = c(-1, 1))
  abline(v = log(quantile(A31se, probs = 0.75)), 
         col=c("blue", "blue"), 
         lty=2, 
         lwd=1)
  
  dev.off()
########################################################################################################  

      
# Comparison of social factors behind tobacco use and COVID-19 vaccination refusal (A4)
########################################################################################################
#ROC-AUC - table 3 model parameters
  summary(A2ROC)[c(1,3,6),]
  summary(A3ROC)[c(1,3,6),]
  wilcox.test(A2ROC, A3ROC)
  wilcoxonZ(A2ROC, A3ROC)
  abs(wilcoxonZ(A2ROC, A3ROC))/ sqrt(1000)
  
#distribution of lambda 1se - table 3 model parameters
  summary(A21se)[c(1,3,6),]
  summary(A31se)[c(1,3,6),]

#table 3 betas
  table3 <- data.frame(
    varname = character(22), 
    a2med   = numeric(22), 
    a2min   = numeric(22),
    a2max   = numeric(22),
    a3med   = numeric(22), 
    a3min   = numeric(22), 
    a3max   = numeric(22)
    )

  for (c3 in 1:22)
    {
    table3$varname[c3] <- paste0("x",c3-1)
    
    table3$a2min[c3]  <- min(A2bet[,c3])
    table3$a2med[c3]  <- median(A2bet[,c3])
    table3$a2max[c3]  <- max(A2bet[,c3])
    
    table3$a3min[c3]  <- min(A3bet[,c3])
    table3$a3med[c3]  <- median(A3bet[,c3])
    table3$a3max[c3]  <- max(A3bet[,c3])
  }


#table3
  table3r <- round(table3[,2:7], digits=1)
  rownames(table3r) <- table3$varname
  table3

#comparing perceived social norms
  wilcox.test(A2bet[,21], A3bet[,22])
  wilcoxonZ(A2bet[,21], A3bet[,22])
  abs(wilcoxonZ(A2bet[,21], A3bet[,22]))/ sqrt(1000)
  
########################################################################################################


# Sensitivity analysis
########################################################################################################
# fit  logistic regression model
  A2AMEm <- glm(y1 ~ x4 + x20 + x21, data = C19, family = binomial)
  A3AMEm <- glm(y2 ~ x4 + x20 + x21, data = C19, family = binomial)

# Calculate the Average Marginal Effects
  A2AME <- margins(A2AMEm, vce="bootstrap")
  A3AME <- margins(A3AMEm, vce="bootstrap")

#Table 4
  summary(A2AME)
  summary(A3AME)

#Nagelkerke y1 ~ x20
  y1x20 <- glm(y1 ~ x20 , family=binomial(link='logit') , data=C19) #y1 and y2 simple modell
  summary(y1x20)  
  exp(coefficients(y1x20)[2])
  exp(confint(y1x20))
  NagelkerkeR2(y1x20)

#ROC-AUC y1 ~ x20
  roc_result <- roc(response = y1x20$model$y1, predictor = predict(y1x20, type = "response"))
  auc_result <- auc(roc_result)
  ci.auc(roc_result, method = "bootstrap", of = "AUC", percent = 0.95)

#Nagelkerke y2 ~ x21
  y2x21 <- glm(y2 ~ x21, family=binomial(link='logit') , data=C19) #y1 and y2 controlled for edu
  summary(y2x21 ) 
  exp(coefficients(y2x21)[2])
  exp(confint(y2x21))
  NagelkerkeR2(y2x21)

#ROC-AUC y2 ~ x21
  roc_result <- roc(response = y2x21$model$y2, predictor = predict(y2x21, type = "response"))
  auc_result <- auc(roc_result)
  ci.auc(roc_result, method = "bootstrap", of = "AUC", percent = 0.95)
#since 95% CI boostrap is stochastic, small difference are to be expected
########################################################################################################
