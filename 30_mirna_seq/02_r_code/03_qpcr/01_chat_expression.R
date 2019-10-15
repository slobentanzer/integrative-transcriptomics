#CHAT EXPRESSION IN CNTF DIFFERENTIATION####
rm(list=ls())

home= '~/GitHub/'
rootdir = paste(home, "integrative-transcriptomics", sep="")
setwd(rootdir)

library(ggplot2)

# LA-N-2 ####
expr <- readRDS("./raw_data/raw_expression_qpcr_la2.rds")

#molecular weight of CNTF: 22,931 Da
mw <- 22931
expr$molarC <- as.numeric(as.character(expr$Concentration))*1E-9/mw*1E3*1E9 #1)ng->g, 2)g->mol, 3)ml->l 4)M->nM
expr$Day <- factor(expr$Day)
expr$Concentration <- factor(expr$Concentration)

ggplot(expr, aes(Day, FoldExpression, fill = Concentration)) + 
  geom_boxplot() +
  scale_fill_brewer(palette = "Set3") +
  scale_y_continuous(breaks = seq(1:10)) + theme(legend.justification=c(0,1), legend.position=c(0,1))

ggsave("./img/cntf_time_dose_boxplot.svg")

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  #if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

expr.fold.sum <- summarySE(expr, measurevar = "FoldExpression", groupvars = c("Concentration", "Day"))

pd <- position_dodge(0.2)
ggplot(expr.fold.sum, aes(Day, FoldExpression, color = Concentration, group = Concentration)) + 
  geom_line(position = pd) + 
  geom_errorbar(aes(ymin=FoldExpression-se, ymax=FoldExpression+se), width=.1, position = pd, color = "grey30") +
  geom_point(position = pd, size=2, shape=21, fill="white") +
  scale_y_continuous(breaks = seq(1:10)) + theme(legend.justification=c(0,1), legend.position=c(0,1))

ggsave("./img/cntf_time_dose_curve.svg")

#pval
t.test(expr$RelExpression[expr$Concentration == "0" & expr$Day == "II"], 
       expr$RelExpression[expr$Concentration == "1" & expr$Day == "II"])
t.test(expr$RelExpression[expr$Concentration == "0" & expr$Day == "II"], 
       expr$RelExpression[expr$Concentration == "10" & expr$Day == "II"])#
t.test(expr$RelExpression[expr$Concentration == "0" & expr$Day == "II"], 
       expr$RelExpression[expr$Concentration == "100" & expr$Day == "II"])#

t.test(expr$RelExpression[expr$Concentration == "0" & expr$Day == "III"], 
       expr$RelExpression[expr$Concentration == "1" & expr$Day == "III"])
t.test(expr$RelExpression[expr$Concentration == "0" & expr$Day == "III"], 
       expr$RelExpression[expr$Concentration == "10" & expr$Day == "III"])#
t.test(expr$RelExpression[expr$Concentration == "0" & expr$Day == "III"], 
       expr$RelExpression[expr$Concentration == "100" & expr$Day == "III"])#

t.test(expr$RelExpression[expr$Concentration == "0" & expr$Day == "IV"], 
       expr$RelExpression[expr$Concentration == "1" & expr$Day == "IV"])
t.test(expr$RelExpression[expr$Concentration == "0" & expr$Day == "IV"], 
       expr$RelExpression[expr$Concentration == "10" & expr$Day == "IV"])
t.test(expr$RelExpression[expr$Concentration == "0" & expr$Day == "IV"], 
       expr$RelExpression[expr$Concentration == "100" & expr$Day == "IV"])#

# LA-N-5 ####
expr <- readRDS("./raw_data/raw_expression_qpcr_la5.rds")

#molecular weight of CNTF: 22,931 Da
mw <- 22931
expr$molarC <- as.numeric(as.character(gsub("ng", "", expr$Concentration)))*1E-9/mw*1E3*1E9 #1)ng->g, 2)g->mol, 3)ml->l 4)M->nM
expr$Day <- factor(expr$Day)
expr$Concentration <- factor(expr$Concentration)

ggplot(expr, aes(Day, RelExpression, fill = Concentration)) + 
  geom_boxplot() +
  scale_fill_brewer(palette = "Set3") +
  scale_y_continuous(breaks = seq(1:10)) + theme(legend.justification=c(0,1), legend.position=c(0,1))

# Supplemental Figure####
ggsave("./img/cntf_time_dose_boxplot_la5.svg")

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  #if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

expr.fold.sum <- summarySE(expr, measurevar = "RelExpression", groupvars = c("Concentration", "Day"))

pd <- position_dodge(0.2)
ggplot(expr.fold.sum, aes(Day, RelExpression, color = Concentration, group = Concentration)) + 
  geom_line(position = pd) + 
  geom_errorbar(aes(ymin=RelExpression-se, ymax=RelExpression+se), width=.1, position = pd, color = "grey30") +
  geom_point(position = pd, size=2, shape=21, fill="white") +
  scale_y_continuous(breaks = seq(1:10)) + theme(legend.justification=c(0,1), legend.position=c(0,1))

ggsave("./img/cntf_time_dose_curve_la5.svg")

#pval
t.test(expr$RelExpression[expr$Concentration == "000ng" & expr$Day == "2"], 
       expr$RelExpression[expr$Concentration == "001ng" & expr$Day == "2"])
t.test(expr$RelExpression[expr$Concentration == "000ng" & expr$Day == "2"], 
       expr$RelExpression[expr$Concentration == "010ng" & expr$Day == "2"])#
t.test(expr$RelExpression[expr$Concentration == "000ng" & expr$Day == "2"], 
       expr$RelExpression[expr$Concentration == "100ng" & expr$Day == "2"])

t.test(expr$RelExpression[expr$Concentration == "000ng" & expr$Day == "4"], 
       expr$RelExpression[expr$Concentration == "001ng" & expr$Day == "4"])
t.test(expr$RelExpression[expr$Concentration == "000ng" & expr$Day == "4"], 
       expr$RelExpression[expr$Concentration == "010ng" & expr$Day == "4"])#
t.test(expr$RelExpression[expr$Concentration == "000ng" & expr$Day == "4"], 
       expr$RelExpression[expr$Concentration == "100ng" & expr$Day == "4"])#

