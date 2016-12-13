#' Code for generating analysis for poster project
#' Wesley Zuidema and Ellen Ahlness
#' MLE Poster Project

# Load library packages
library(ggplot2)
library(knitr)
library(MASS)
library(nlme)
library(boot)            # For cv.glm()
library(separationplot)  # For separation plot
library(pscl)            # Alternative PCP code
library(verification)    # For ROC area
library(tile)            # For some graphics; used by plot.binPredict()
library(RColorBrewer)    # For nice colors
source("binaryGOF.R")    # Percent correctly predicted and concordance indexes
source("binPredict.R")   # Code for making predicted vs actual plots
library(simcf)
library(xtable)
library(VGAM)
library(Amelia)
library(dplyr)
select <- dplyr::select


# Load data
load("data/wvs_data.rdata")
wvs_data <- WVS_Longitudinal_1981_2014_R_v2015_04_18
rm(WVS_Longitudinal_1981_2014_R_v2015_04_18)

# Rename variables we are using to make them easier to read
wvs_data$wave <- as.factor(wvs_data$S002)
wvs_data$country <- as.factor(wvs_data$S003)

# Have data for 20 - Andorra, 36 - Australia, 348 - Hungary, 360 - Indonesia, 
# 484 - Mexico, 840 - USA, and 858 - Uraguay
wvs_data$partypref <- wvs_data$E256

# Create binary variable to indicate support of far right party
# Leave NAs for missing so that we can impute
wvs_australia <- filter(wvs_data, country == 36)
wvs_hungary <- filter(wvs_data, country == 348)
wvs_indonesia <- filter(wvs_data, country == 360)

wvs_filter <- rbind(wvs_australia, wvs_hungary, wvs_indonesia)
wvs_filter <- wvs_filter %>%
  mutate(partypref = ifelse(E256 <= 0, NA, E256))

wvs_filter$farrightpref <- 0
wvs_filter <- wvs_filter %>%
  mutate(farrightpref = ifelse(partypref == 36003 | partypref == 36006 | partypref == 36008 | partypref == 36009 | partypref == 348001 | partypref == 348004 | partypref == 360007, 1, 0))


# Recode data to change non-responses to NA for multiple imputation

wvs_recode <- wvs_filter %>%
  mutate(immigrationpref = ifelse(E143 <= 0, NA, E143)) %>%
  mutate(age = ifelse(X003 <= 0, NA, X003)) %>%
  mutate(incomelevel = ifelse(X047 <= 0, NA, X047)) %>%
  mutate(female = ifelse(X001 <= 0, NA, X001)) %>%
  mutate(maritalstatus = ifelse(X007 <= 0, NA, X007)) %>%
  mutate(education = ifelse(X025R <= 0, NA, X025R)) %>%
  mutate(religious = ifelse(F034 <= 0, NA, F034)) %>%
  mutate(worldcitizen = ifelse(G019 <= 0, NA, G019)) %>%
  mutate(nationalpride = ifelse(G006 <= 0, NA, G006)) %>%
  mutate(infofriends = ifelse(E254 <= 0, NA, E254))

wvs_recode$female[wvs_recode$female == 1] <- 0
wvs_recode$female[wvs_recode$female == 2] <- 1

# Recode auxiliary variables for multiple imputation
wvs_recode <- wvs_recode %>%
  mutate(aux1 = ifelse(D059 <= 0, NA, D059),
         aux2 = ifelse(E007 <= 0, NA, E007),
         aux3 = ifelse(E023 <= 0, NA, E023),
         aux4 = ifelse(E046 <= 0, NA, E046),
         aux5 = ifelse(F108 <= 0, NA, F108)
         )


# Select variables for imputation model
wvs_observed <- wvs_recode %>%
  select(wave, country, farrightpref, immigrationpref, age, incomelevel, female, maritalstatus, education, religious, worldcitizen, nationalpride, aux1, aux2, aux3, aux4, aux5)

# Define all potential models here!
# All models include country fixed effects, Indonesia is baseline
model_a <-  farrightpref ~ immigrationpref + australia + hungary #Binary Model
model_b <-  farrightpref ~ immigrationpref + religious + worldcitizen + nationalpride + australia + hungary#Preferences
model_c <- farrightpref ~ age + incomelevel + female + education + married + domesticpartner + divorced + separated + widowed + australia + hungary #Characteristics
model_d <- farrightpref ~ immigrationpref + age + incomelevel + female + education + religious + worldcitizen + nationalpride + married + domesticpartner + divorced + separated + widowed + australia + hungary #Full Model

### Cross Validation
loocv <- function (obj) {
  data <- obj$data
  m <- dim(data)[1]
  form <- formula(obj)
  fam <- obj$family$family
  loo <- rep(NA, m)
  for (i in 1:m) {
    i.glm <- glm(form, data = data[-i, ], family = fam)
    loo[i] <- predict(i.glm, newdata = data[i,], family = fam, type = "response")
  }
  loo
}


##### IMPUTATION
nimp <- 5 # Use nimp=5 at minimum; 10 often a good idea
ordinal <- c("maritalstatus", "farrightpref", "immigrationpref", "age", "incomelevel", "education", "religious", "worldcitizen", "nationalpride", "aux1", "aux2", "aux3", "aux4", "aux5")
nominal <- c("country")
amelia.res <- amelia(wvs_observed, m=nimp, ords = ordinal, noms = nominal, idvars = "wave")
miData <- amelia.res$imputation

# Impuation diagnostic
# overimpute(amelia.res, var="incomelevel")


mi_data <- vector("list", nimp)
mi_results <- vector("list", nimp)
mi_loocv <- vector("list", nimp)


# Loop through analysis
for (i in 1:nimp) {
  # Break out dummy variables
  mi_data[[i]] <- miData[[i]] %>%
    mutate(married = ifelse(maritalstatus == 1, 1, 0),
           domesticpartner = ifelse(maritalstatus == 2, 1, 0),
           divorced = ifelse(maritalstatus == 3, 1, 0),
           separated = ifelse(maritalstatus == 4, 1, 0),
           widowed = ifelse(maritalstatus == 5, 1, 0),
           nevermarried = ifelse(maritalstatus == 6, 1, 0),
           australia = ifelse(country == 36, 1, 0),
           hungary = ifelse(country == 348, 1, 0))
  # Run Analysis
  mi_results[[i]] <- glm(model_d, family = binomial, data = mi_data[[i]])
  mi_loocv[[i]] <- loocv(mi_results[[i]])
}
# Draw Simulates
sims <- 10000 
simbetas <- NULL 
for (i in 1:nimp) { 
  simbetas <- rbind(simbetas, mvrnorm(sims/nimp, coef(mi_results[[i]]), vcov(mi_results[[i]])))
}

# Looped Analysis for Preferences only
mi_results_pref <- vector("list", nimp)
mi_loocv_pref <- vector("list", nimp)
for (i in 1:nimp) {
  # Break out dummy variables
  mi_data[[i]] <- miData[[i]] %>%
    mutate(married = ifelse(maritalstatus == 1, 1, 0),
           domesticpartner = ifelse(maritalstatus == 2, 1, 0),
           divorced = ifelse(maritalstatus == 3, 1, 0),
           separated = ifelse(maritalstatus == 4, 1, 0),
           widowed = ifelse(maritalstatus == 5, 1, 0),
           nevermarried = ifelse(maritalstatus == 6, 1, 0),
           australia = ifelse(country == 36, 1, 0),
           hungary = ifelse(country == 348, 1, 0))
  # Run Analysis
  mi_results_pref[[i]] <- glm(model_b, family = binomial, data = mi_data[[i]])
  mi_loocv_pref[[i]] <- loocv(mi_results_pref[[i]])
}
# Draw Simulates
sims <- 10000 
simbetas <- NULL 
for (i in 1:nimp) { 
  simbetas <- rbind(simbetas, mvrnorm(sims/nimp, coef(mi_results_pref[[i]]), vcov(mi_results_pref[[i]])))
}

# Looped Analysis for Characteriscs only
mi_results_char <- vector("list", nimp)
mi_loocv_char <- vector("list", nimp)
for (i in 1:nimp) {
  # Break out dummy variables
  mi_data[[i]] <- miData[[i]] %>%
    mutate(married = ifelse(maritalstatus == 1, 1, 0),
           domesticpartner = ifelse(maritalstatus == 2, 1, 0),
           divorced = ifelse(maritalstatus == 3, 1, 0),
           separated = ifelse(maritalstatus == 4, 1, 0),
           widowed = ifelse(maritalstatus == 5, 1, 0),
           nevermarried = ifelse(maritalstatus == 6, 1, 0),
           australia = ifelse(country == 36, 1, 0),
           hungary = ifelse(country == 348, 1, 0))
  # Run Analysis
  mi_results_char[[i]] <- glm(model_c, family = binomial, data = mi_data[[i]])
  mi_loocv_char[[i]] <- loocv(mi_results_char[[i]])
}
# Draw Simulates
sims <- 10000 
simbetas <- NULL 
for (i in 1:nimp) { 
  simbetas <- rbind(simbetas, mvrnorm(sims/nimp, coef(mi_results_char[[i]]), vcov(mi_results_char[[i]])))
}

# Looped Analysis for Immigration only
mi_results_binary <- vector("list", nimp)
mi_loocv_binary <- vector("list", nimp)
for (i in 1:nimp) {
  # Break out dummy variables
  mi_data[[i]] <- miData[[i]] %>%
    mutate(married = ifelse(maritalstatus == 1, 1, 0),
           domesticpartner = ifelse(maritalstatus == 2, 1, 0),
           divorced = ifelse(maritalstatus == 3, 1, 0),
           separated = ifelse(maritalstatus == 4, 1, 0),
           widowed = ifelse(maritalstatus == 5, 1, 0),
           nevermarried = ifelse(maritalstatus == 6, 1, 0),
           australia = ifelse(country == 36, 1, 0),
           hungary = ifelse(country == 348, 1, 0))
  # Run Analysis
  mi_results_binary[[i]] <- glm(model_a, family = binomial, data = mi_data[[i]])
  mi_loocv_binary[[i]] <- loocv(mi_results_binary[[i]])
}
# Draw Simulates
sims <- 10000 
simbetas <- NULL 
for (i in 1:nimp) { 
  simbetas <- rbind(simbetas, mvrnorm(sims/nimp, coef(mi_results_binary[[i]]), vcov(mi_results_binary[[i]])))
}


# Choose nice colors for plots
col <- brewer.pal(5, "Set1")
blue <- col[2]
orange <- col[5]


# Percent Correctly Predicted
pcp.null1 <- pcp.glm(mi_results[[1]], mi_data[[1]]$farrightpref, type="null")
pcp.null2 <- pcp.glm(mi_results[[2]], mi_data[[2]]$farrightpref, type="null")
pcp.null3 <- pcp.glm(mi_results[[3]], mi_data[[3]]$farrightpref, type="null")
pcp.null4 <- pcp.glm(mi_results[[4]], mi_data[[4]]$farrightpref, type="null")
pcp.null5 <- pcp.glm(mi_results[[5]], mi_data[[5]]$farrightpref, type="null")
pcp.1 <- pcp.glm(mi_results[[1]], mi_data[[1]]$farrightpref, type="model")
pcp.2 <- pcp.glm(mi_results[[2]], mi_data[[2]]$farrightpref, type="model")
pcp.3 <- pcp.glm(mi_results[[3]], mi_data[[3]]$farrightpref, type="model")
pcp.4 <- pcp.glm(mi_results[[4]], mi_data[[4]]$farrightpref, type="model")
pcp.5 <- pcp.glm(mi_results[[5]], mi_data[[5]]$farrightpref, type="model")
pcpi.1 <- pcp.glm(mi_results[[1]], mi_data[[1]]$farrightpref, type="improve")
pcpi.2 <- pcp.glm(mi_results[[2]], mi_data[[2]]$farrightpref, type="improve")
pcpi.3 <- pcp.glm(mi_results[[3]], mi_data[[3]]$farrightpref, type="improve")
pcpi.4 <- pcp.glm(mi_results[[4]], mi_data[[4]]$farrightpref, type="improve")
pcpi.5 <- pcp.glm(mi_results[[5]], mi_data[[5]]$farrightpref, type="improve")

separationplot(pred=mi_loocv_binary[[1]], actual=mi_results_binary[[1]]$y,
               file="figures/MI_SepPlot_Binary.pdf")

separationplot(pred=mi_loocv[[1]], actual=mi_results[[1]]$y,
               file="figures/MI_SepPlot_Full.pdf")


# Regression Table
htmlreg(l = list(mi_results[[1]], mi_results[[2]], mi_results[[3]], mi_results[[4]], mi_results[[5]]), 
        stars = numeric(0),
        custom.model.names = c("MI Data 1", "MI Data 2", "MI Data 3", "MI Data 4", "MI Data 5"),
        file = "figures/MIregtable.doc")

htmlreg(l = list(mi_results_binary[[1]], mi_results_pref[[1]], mi_results_char[[1]], mi_results[[1]]), 
        stars = numeric(0),
        file = "figures/modelcomparisonregtable.doc")



       
## Simulations
# Model parameters
pe.glm <- mi_results[[1]]$coefficients  # point estimates
vc.glm <- vcov(mi_results[[1]])         # var-cov matrix

# Simulate parameter distributions
sims <- 10000
simbetas <- mvrnorm(sims, pe.glm, vc.glm)

# Create example counterfactuals
xhyp <- cfMake(model_d, mi_data[[1]], nscen=10)

xhyp <- cfName(xhyp,"Opposes Immigration", scen=1)
xhyp <- cfChange(xhyp, "immigrationpref", x=4, xpre=1, scen=1)

xhyp <- cfName(xhyp, "National Pride", scen=2)
xhyp <- cfChange(xhyp, "nationalpride", x=1, xpre=4, scen=2)

xhyp <- cfName(xhyp, "World Citizen", scen=3)
xhyp <- cfChange(xhyp, "worldcitizen", x=1, xpre=4, scen=3)

xhyp <- cfName(xhyp, "Religious", scen=4)
xhyp <- cfChange(xhyp, "religious", x=1, xpre=3, scen=4)

xhyp <- cfName(xhyp, "Age + 1sd = 56", scen=5)
xhyp <- cfChange(xhyp, "age",
                 x=mean(na.omit(mi_data[[1]]$age))+sd(na.omit(mi_data[[1]]$age)),
                 xpre=mean(na.omit(mi_data[[1]]$age)),
                 scen=5)

xhyp <- cfName(xhyp, "Age - 1sd = 24", scen=6)
xhyp <- cfChange(xhyp, "age",
                 x=mean(na.omit(mi_data[[1]]$age))-sd(na.omit(mi_data[[1]]$age)),
                 xpre=mean(na.omit(mi_data[[1]]$age)),
                 scen=6)

xhyp <- cfName(xhyp, "Male", scen=7)
xhyp <- cfChange(xhyp, "female", x=1, xpre=2, scen=7)

xhyp <- cfName(xhyp, "Income (1st to 3rd QT)", scen=8)
xhyp <- cfChange(xhyp, "incomelevel", x=7, xpre=3, scen=8)

xhyp <- cfName(xhyp, "Widowed", scen=9)
xhyp <- cfChange(xhyp, "widowed", x=1, xpre=0, scen=9)

xhyp <- cfName(xhyp, "Divorced", scen=10)
xhyp <- cfChange(xhyp, "female", x=1, xpre=0, scen=10)

# Simulate expected probabilities for all scenarios
logitsims <- logitsimev(xhyp, simbetas, ci=0.95)
logitfds <- logitsimfd(xhyp, simbetas, ci=0.95)

# Get 3 nice colors for traces
col <- brewer.pal(3,"Dark2")

# Plot predicted probabilities for all four categories, sorted by size
sorted <- order(logitsims$pe)
scenNames <- row.names(xhyp$x)

trace1 <- ropeladder(x = logitsims$pe[sorted],
                     lower = logitsims$lower[sorted],
                     upper = logitsims$upper[sorted],
                     labels = scenNames[sorted],
                     size=0.5,
                     lex=1.5,
                     lineend="square",
                     plot=1
)

tile(trace1,
     limits = c(0,0.25),
     gridlines = list(type="xt"),
     topaxis=list(add=TRUE, at=c(0, 0.05, .1, .15, .2)),
     plottitle=list(labels="Expected Values"),
     xaxistitle=list(labels="Probability of far right support"),
     width=list(spacer=3),
     height = list(plottitle=3,xaxistitle=3.5,topaxistitle=3.5),
     output = list(file="figures/MIexpected_values")
)

# Plot First Difference
sorted <- order(logitfds$pe)
scenNames <- row.names(xhyp$x)
# Set up lineplot trace of relative risk
trace2 <- ropeladder(x = logitfds$pe[sorted],
                     lower = logitfds$lower[sorted],
                     upper = logitfds$upper[sorted],
                     labels = row.names(xhyp$x[sorted]),
                     size=0.5,
                     lex=1.5,
                     lineend="square",
                     baseline=list(at=0),
                     plot=1
)


# Plot traces using tile
tile(trace2,
     limits = c(-.1,0.1),
     gridlines = list(type="xt"),
     topaxis=list(add=TRUE, at=c(-.05, 0, .05, .1)),
     plottitle=list(labels="First Differences"),
     xaxistitle=list(labels="Probability of far right support"),
     width=list(spacer=3),
     height = list(plottitle=3,xaxistitle=3.5,topaxistitle=3.5),
     output=list(file="figures/MIfirstdifferences")
)




# Make cross-validated AVP and ROC plots; note use of newpred input in binPredict
binnedM1cv <- binPredict(mi_results_binary[[1]],  newpred=mi_loocv_binary[[1]],col=orange, label="M1: Binary Model", quantiles=TRUE, bins = 5)
binnedM4cv <- binPredict(mi_results[[1]], newpred=mi_loocv[[1]], col=orange, label="M4: Full Model", quantiles=TRUE, bins = 5)

binnedM2cv <- binPredict(mi_results_pref[[1]],  newpred=mi_loocv_pref[[1]],col=blue, label="M2: Preferences Only", quantiles=TRUE, bins = 5)
binnedM3cv <- binPredict(mi_results_char[[1]],  newpred=mi_loocv_char[[1]],col=blue, label="M3: Characteristics Only", quantiles=TRUE, bins = 5)

plot(binnedM2cv, binnedM4cv, display=c("roc"), hide=TRUE, thresholds=c(0.9, 0.8, 0.7, 0.6),
     labx=0.25, file = "figures/MI - ROC1 - Preferences vs FULL")

plot(binnedM2cv, binnedM4cv, display=c("avp"), hide=TRUE, thresholds=c(0.9, 0.8, 0.7, 0.6),
     labx=0.25, file = "figures/MI - AVP1 - Preferences vs FULL")

plot(binnedM3cv, binnedM4cv, display=c("roc"), hide=TRUE, thresholds=c(0.9, 0.8, 0.7, 0.6),
     labx=0.25, file = "figures/MI - ROC2 - Characateristics vs FULL")

plot(binnedM3cv, binnedM4cv, display=c("avp"), hide=TRUE, thresholds=c(0.9, 0.8, 0.7, 0.6),
     labx=0.25, file = "figures/MI - AVP2 - Characateristics vs FULL")

plot(binnedM2cv, binnedM4cv, display=c("roc", "avp"), hide=TRUE, thresholds=c(0.9, 0.8, 0.7, 0.6),
     labx=0.25, file = "figures/MI - Preferences vs FULL")
plot(binnedM3cv, binnedM4cv, display=c("roc", "avp"), hide=TRUE, thresholds=c(0.9, 0.8, 0.7, 0.6),
     labx=0.25, file = "figures/MI - Characateristics vs FULL")


## Testing for heterogenous age effects with vgam
model_vgam <- farrightpref ~ immigrationpref + s(age) + incomelevel + female + education + religious + worldcitizen + nationalpride + married + domesticpartner + divorced + separated + widowed + australia + hungary #Full Model

mi_vgam <- vector("list", nimp)
for (i in 1:nimp) {
  # Break out dummy variables
  mi_data[[i]] <- miData[[i]] %>%
    mutate(married = ifelse(maritalstatus == 1, 1, 0),
           domesticpartner = ifelse(maritalstatus == 2, 1, 0),
           divorced = ifelse(maritalstatus == 3, 1, 0),
           separated = ifelse(maritalstatus == 4, 1, 0),
           widowed = ifelse(maritalstatus == 5, 1, 0),
           nevermarried = ifelse(maritalstatus == 6, 1, 0),
           australia = ifelse(country == 36, 1, 0),
           hungary = ifelse(country == 348, 1, 0))
  # Run Analysis
  mi_vgam[[i]] <- vgam(model_vgam, binomialff, data = mi_data[[i]])
}

summary(mi_vgam[[1]])
summary(mi_vgam[[2]])
summary(mi_vgam[[3]])
summary(mi_vgam[[4]])
summary(mi_vgam[[5]])

pfit <- predict(fit, type = "terms", raw = TRUE, se = TRUE)
head(pfit$fitted)
head(pfit$se.fit)
pfit$df
pfit$sigma
