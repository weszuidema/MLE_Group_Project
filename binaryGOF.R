# Four helpful goodness of fit functions
# for binary choice and ordered probit
# Chris Adolph
# chrisadolph.com
#
# (in the below, the null model assumes the model category always occurs)
#
# 1. percent correctly predicted, binary glm
#
# pcp.glm(res,             # an estimated glm object
#         y,               # the response data as a vector
#         type = "model"   # or "null" for the null model PCP
#         )                #   or "improve" for degree of improvement from null model
#
# 2. percent correctly predicted, ordered probit (requires simcf package from my website)
#
# pcp.oprobit(x,               # matrix of covariates (leave out the constant)
#             y,               # the response data as a vector
#             b,               # ordered probit coefficients, with cutpoints at the end
#             constant=1,      # location of the constant in the b vector;
#                                 if you used polr(), set to NA for no constant
#             ncat=3,          # number of possible outcomes
#             type = "model"   # or "null" for the null model PCP
#             )                #   or "improve" for degree of improvement from null model
#
# 3. concordance index (area under ROC), binary glm
#     - more useful for rare events than PCP
#
# concord.glm(res,             # an estimated glm object
#             y,               # the response data as a vector
#             type = "model"   # or "null" for the null model concordance
#             )                #   or "improve" for degree of improvement from null model
#
# 4. concordance index (area under ROC), ordered probit (requires simcf package from my website)
#     - more useful for rare events than PCP
#
# concord.oprobit(x,               # matrix of covariates (leave out the constant)
#                 y,               # the response data as a vector
#                 b,               # ordered probit coefficients, with cutpoints at the end
#                 constant=1,      # location of the constant in the b vector;
#                                     if you used polr(), set to NA for no constant
#                 ncat=3,          # number of possible outcomes
#                 type = "model"   # or "null" for the null model concordance
#                 )                #   or "improve" for degree of improvement from null model
#
pcp.glm <- function(res, y, type="model") { # other types:  null, improve
  
  pcp <- mean(round(predict(res,type="response"))==y)
  pcpNull <- max(mean(y), mean(1-y))
  pcpImprove <- (pcp-pcpNull)/(1-pcpNull)
  
  if (type=="model")
    return(pcp)
  if (type=="null")
    return(pcpNull)
  if (type=="improve")
    return(pcpImprove)
}

pcp.oprobit <- function(x, y, b, constant=1, ncat=3, type="model") { # other types:  null, improve
  
  require(simcf)
  b <- matrix(b,nrow=100,ncol=length(b),byrow=TRUE)
  simy <- oprobitsimev(x, b, constant=constant, cat=ncat)
  
  cats <- sort(unique(y))
  
  predcatN <- cats[rev(order(table(y)))][1]
  
  n <- length(y)
  pcp <- pcpNull <- predcatM <- rep(NA,n)
  for (i in 1:n) {
    predcatM[i] <- cats[rev(order(simy$pe[i,]))][1]
    pcp[i] <- predcatM[i]==y[i]
    pcpNull[i] <- predcatN==y[i]
  }
  
  pcp <- mean(pcp)
  pcpNull <- mean(pcpNull)
  pcpImprove <- (pcp-pcpNull)/(1-pcpNull)
  
  if (type=="model")
    return(pcp)
  if (type=="null")
    return(pcpNull)
  if (type=="improve")
    return(pcpImprove)
  
}


concord.glm <- function(res, y, type="model") { # other types:  null, improve
  require(verification)
  if (type=="model") {
    yhat <- predict(res,type="response")
    concord <- roc.area(y, yhat)$A
  }
  if (type=="null") {    
    yhat <- rep(max(mean(y), mean(1-y)), length(y))
    concord <- roc.area(y, yhat)$A
  }
  if (type=="improve") {
    yhat <- predict(res,type="response")
    model <- roc.area(y, yhat)$A
    yhat <- rep( max(mean(y), mean(1-y)), length(y))
    null <- roc.area(y, yhat)$A
    concord <- (model-null)/(1-null)
  }
  concord
}


concord.oprobit <- function(x, y, b, constant=1, ncat=3, type="model") { # other types:  null, improve
  require(simcf)
  require(verification)
  b <- matrix(b,nrow=100,ncol=length(b),byrow=TRUE)
  simy <- oprobitsimev(x, b, constant=constant, cat=ncat)
  cats <- sort(unique(y))
  
  if (type=="model") {
    model <- rep(NA,ncat)
    for (j in 1:ncat) {
      yhat <- simy$pe[,j]
      model[j] <- roc.area(as.numeric(y==cats[j]), yhat)$A
    }
    concord <- mean(model) 
  }
  
  if (type=="null") {
    null <- rep(NA,ncat)
    for (j in 1:ncat) {
      probs <- rep(mean(y==cats[j]), length(y))      
      null[j] <- roc.area(as.numeric(y==cats[j]), probs)$A
    }
    concord <- mean(null) 
  }
  
  if (type=="improve") {
    improve <- rep(NA,ncat)
    for (j in 1:ncat) {
      probs <- rep(mean(y==cats[j]), length(y))
      null <- roc.area(as.numeric(y==cats[j]), probs)$A
      yhat <- simy$pe[,j]
      model <- roc.area(as.numeric(y==cats[j]), yhat)$A
      improve[j] <- (model-null)/(1-null)
    }
    concord <- mean(improve) 
  }
  
  concord
}

