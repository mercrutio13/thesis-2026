library(tidyverse)

cleandata <- function(datafile) {

  rdata <- read.csv(datafile)
  rdata |> 
    filter(Response.Type == "response", Display == "Trial") |>
    mutate(
      pid = Participant.Private.ID,
      conf = Spreadsheet..setting,
      stimulusfile = case_when(
        conf == "DPDM" ~ Spreadsheet..DPDM,
        conf == "DPNM" ~ Spreadsheet..DPNM,
        conf == "NPDM" ~ Spreadsheet..NPDM,
        conf == "NPNM" ~ Spreadsheet..NPNM
      ),
      p = str_extract(Spreadsheet..setting, ".(?=P)"),
      phonology = case_when(
        p == "N" ~ 0,
        p == "D" ~ 1,
      ),
      m = str_extract(Spreadsheet..setting, ".(?=M)"),
      morphology = case_when(
        m == "N" ~ 0,
        m == "D" ~ 1,
      ),
      sindex = as.integer(str_extract(stimulusfile, "\\d+")),
      sconst = (((sindex - 1) %% 18) %/% 3) + 1,
      speaker = str_extract(stimulusfile, ".(?=\\.mp3$)"),
      correct = as.logical(Correct)) |>
    group_by(pid) |>
    mutate(id = cur_group_id()) |>
    ungroup() |>
    select(id,
           pid,
           sindex,
           sconst,
           conf,
           phonology,
           morphology,
           speaker,
           Reaction.Time, 
           correct)
}

remove_errdata <- function(datain) {
  errdata <- datain |>
    group_by(sindex, conf) |>
    summarise(n = n(), err = (n() - sum(correct)), .groups = "drop_last") |>
    filter(conf == "NPNM" & (err > 1 | (err == 1 & n == 1)))
  anti_join(datain, 
            errdata,
            by = join_by(sindex))
}

#Load data and demographic data
data <- cleandata("data-raw.csv")

#Plot Reaction Time by participant, consider outliers
data |>
  ggplot(aes(group = id, y = Reaction.Time)) +
  geom_boxplot()

#Look for high error rates in NPNM
data |>
  group_by(sindex, conf) |>
  summarise(pcorr = sum(correct)/n(), 
            err = (n() - sum(correct)), n = n(), 
            meanrtime = mean(Reaction.Time)) |>
  filter(conf == "NPNM", err > 0) |>
  arrange(pcorr, desc(err))

#Remove sentences with high NPNM error rates
nerrdata <- remove_errdata(data)

#Correctness Analysis 
library(car)

mpcorrdata <- nerrdata |>
  group_by(id, conf) |>
  summarise(x = sum(correct),
            xp = x/n(),
            t = asin(sqrt(x/(n()+1))) + asin(sqrt((x+1)/(n()+1))),
            r = 46.47324337*t - 23,
            .groups = "drop_last") |>
  mutate(p = str_extract(conf, ".(?=P)"),
         phonology = case_when(
           p == "N" ~ 0,
           p == "D" ~ 1,
         ),
         m = str_extract(conf, ".(?=M)"),
         morphology = case_when(
           m == "N" ~ 0,
           m == "D" ~ 1,
         ),
        )

library(ggstatsplot)
mpcorrdata |>
  ggplot(aes(x = conf, y = xp, fill = conf)) +
  geom_boxplot(alpha = 0.5) +
  labs(x = "Configuration", y = "Mean percentage correct") +
  scale_y_continuous(label = function(x) {return(paste(x*100, "%"))}) +
  scale_x_discrete(limits=rev) +
  theme_ggstatsplot() +
  theme(legend.position = "none",
        axis.title = element_text(color = "darkred", size = 14)) +
  scale_fill_brewer(palette="Dark2")

  #Two-way ANOVA for phonology/morphology
summary(aov(r ~ phonology * morphology, data = mpcorrdata)) 
Anova(aov(r ~ phonology + morphology, data = mpcorrdata), type = 2) 

#Response Time Analysis
  #One-way ANOVA for correctness
summary(aov(Reaction.Time ~ correct, data = nerrdata))

nerrdata |>
  ggplot(aes(x = conf, y = Reaction.Time)) +
  geom_boxplot(mapping = aes(fill = correct), alpha =0.5)+
  labs(x = "Configuration", y = "Mean response time", fill = "Response") +
  scale_y_continuous(label = function(x) {return(paste(x, " ms"))}) +
  scale_x_discrete(limits=rev) +
  theme_ggstatsplot() +
  theme(axis.title = element_text(color = "darkred", size = 14)) +
  scale_fill_brewer(palette="Dark2", labels = c("Incorrect", "Correct"))

#Linear mixed-effect regression
corrdata <- filter(nerrdata, correct == TRUE)

library(lme4)
library(merTools)
library(lmerTest)
library(MuMIn)

corr_model <- lmer(Reaction.Time ~ 
                     phonology + morphology 
                      + (1 | id) + (1 | sconst / sindex), 
                   data = corrdata)
summary(corr_model)
fixef(corr_model)
plot(corr_model)
plotREsim(REsim(corr_model))
plotFEsim(FEsim(corr_model))
r.squaredGLMM(corr_model)

ggcoefstats(corr_model, exclude.intercept = TRUE, 
            title = "lmer Linear mixed-effects regression coefficients", 
            ggtheme = theme_ggstatsplot())

library(ggpubr)
ggqqplot(corrdata$Reaction.Time)
shapiro.test(corrdata$Reaction.Time)

#Robust linear mixed-effect regression
library(robustlmm)
library(insight)

getR2 <- function(modelin) {
  var.fix <- get_variance_fixed(modelin)
  var.ran <- get_variance_random(modelin)
  var.res <- get_variance_residual(modelin)
  R2m = var.fix/(var.fix+var.ran+var.res)
  R2c = (var.fix+var.ran)/(var.fix+var.ran+var.res)
  
  c("R2m" = R2m, "R2c" = R2c)
}

rcorr_model <- rlmer(Reaction.Time ~ 
                       phonology + morphology 
                        + (1 | id) + (1 | sconst / sindex), 
                     data = corrdata)
summary(rcorr_model)
fixef(rcorr_model)
plot(rcorr_model)
getR2(rcorr_model)

#Log transformed data
corr_log_model <- lmer(log(Reaction.Time) ~ 
                         phonology + morphology 
                          + (1 | id) + (1 | sconst / sindex), 
                       data = corrdata)
summary(corr_log_model)
fixef(corr_log_model)
plot(corr_log_model)
getR2(corr_log_model)

ggqqplot(log(corrdata$Reaction.Time))

rcorr_log_model <- rlmer(log(Reaction.Time) ~ 
                           phonology + morphology 
                            + (1 | id) + (1 | sconst / sindex), 
                         data = corrdata)
summary(rcorr_log_model)
fixef(rcorr_log_model)
plot(rcorr_log_model)
getR2(rcorr_log_model)

models <- list(
  "lmer" = corr_model,
  "rlmer" = rcorr_model,
  "lmer(log)" = corr_log_model,
  "rlmer(log)" = rcorr_log_model
)

library(modelsummary)
library(tinytable)
modelsummary(models, stars = TRUE, gof_omit = "IC|Adj|F|RMSE|Log") |>
  group_tt(j = list("Linear" = 2:3, "Log-Linear" = 4:5)) 

modelplot(models[1:2], coef_omit = "Intercept|SD", 
          title = "Coefficient estimates and 95% confidence intervals for linear models")
modelplot(models[3:4], coef_omit = "Intercept|SD", 
          title = "Coefficient estimates and 95% confidence intervals for linear models")