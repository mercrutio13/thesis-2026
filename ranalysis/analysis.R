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
      index = as.integer(str_extract(stimulusfile, "\\d+")),
      sconst = (((index - 1) %% 18) %/% 3) + 1,
      construction = case_when(
        sconst == 1 ~ "Time preposition",
        sconst == 2 ~ "Subject relative clause",
        sconst == 3 ~ "Perfect auxilliary",
        sconst == 4 ~ "Definite plural",
        sconst == 5 ~ "Indefinite plural",
        sconst == 6 ~ "Definite noun phrase",
      ),
      speaker = str_extract(stimulusfile, ".(?=\\.mp3$)"),
      correct = as.logical(Correct)) |>
    group_by(pid) |>
    mutate(id = cur_group_id()) |>
    ungroup() |>
    select(id,
           index,
           construction,
           conf,
           phonology,
           morphology,
           speaker,
           Reaction.Time, 
           correct)
}

remove_errdata <- function(datain) {
  errdata <- datain |>
    group_by(index, conf) |>
    summarise(n = n(), err = (n() - sum(correct)), .groups = "drop_last") |>
    filter(conf == "NPNM" & (err > 1 | (err == 1 & n == 1)))
  anti_join(datain, 
            errdata,
            by = join_by(index))
}

#Load data
data <- cleandata("data-raw.csv")

#Plot Reaction Time by participant, consider outliers
data |>
  ggplot(aes(group = id, y = Reaction.Time)) +
  geom_boxplot()

#Look for high error rates in NPNM
data |>
  group_by(index, conf) |>
  summarise(pcorr = sum(correct)/n(), 
            err = (n() - sum(correct)), n = n(), 
            meanrtime = mean(Reaction.Time)) |>
  filter(conf == "NPNM", err > 0) |>
  arrange(pcorr, desc(err))

#Remove sentences with high NPNM error rates
nerrdata <- remove_errdata(data)

#Correctness Analysis 
library(car)
library(ggstatsplot)

#Two-way ANOVA for phonology/morphology
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
         )
  )

summary(aov(r ~ phonology * morphology, data = mpcorrdata)) 
Anova(aov(r ~ phonology + morphology, data = mpcorrdata), type = 2) 

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

# One-way ANOVA for Norweigan speakers
mp_speaker_n <- nerrdata |>
  group_by(id, speaker) |>
  summarise(n = n(),
            x = sum(correct),
            xp = x/n(),
            t = asin(sqrt(x/(n()+1))) + asin(sqrt((x+1)/(n()+1))),
            r = 46.47324337*t - 23,
            .groups = "drop_last") |>
  filter(speaker %in% c("e", "p"))
summary(aov(r ~ speaker, data = mp_speaker_n))

mp_speaker_n |>
  ggplot(aes(x = speaker, y = xp, fill = speaker)) +
  geom_boxplot(alpha = 0.5) +
  labs(x = "Configuration", y = "Mean percentage correct") +
  scale_y_continuous(label = function(x) {return(paste(x*100, "%"))}) +
  scale_x_discrete(limits=rev) +
  theme_ggstatsplot() +
  theme(legend.position = "none",
        axis.title = element_text(color = "darkred", size = 14)) +
  scale_fill_brewer(palette="Dark2")

# One-way ANOVA for Danish speakers
mp_speaker_d <- nerrdata |>
  group_by(id, speaker) |>
  summarise(n = n(),
            x = sum(correct),
            xp = x/n(),
            t = asin(sqrt(x/(n()+1))) + asin(sqrt((x+1)/(n()+1))),
            r = 46.47324337*t - 23,
            .groups = "drop_last") |>
  filter(speaker %in% c("m", "s"))
summary(aov(r ~ speaker, data = mp_speaker_d))

mp_speaker_d |>
  ggplot(aes(x = speaker, y = xp, fill = speaker)) +
  geom_boxplot(alpha = 0.5) +
  labs(x = "Configuration", y = "Mean percentage correct") +
  scale_y_continuous(label = function(x) {return(paste(x*100, "%"))}) +
  scale_x_discrete(limits=rev) +
  theme_ggstatsplot() +
  theme(legend.position = "none",
        axis.title = element_text(color = "darkred", size = 14)) +
  scale_fill_brewer(palette="Dark2")

# One-way ANOVA for construction
mp_const <- nerrdata |>
  group_by(id, construction) |>
  summarise(n = n(),
            x = sum(correct),
            xp = x/n(),
            t = asin(sqrt(x/(n()+1))) + asin(sqrt((x+1)/(n()+1))),
            r = 46.47324337*t - 23,
            .groups = "drop_last") 

summary(aov(r ~ construction, data = mp_const))

mp_const |>
  ggplot(aes(x = construction, y = xp, fill = construction)) +
  geom_boxplot(alpha = 0.5) +
  labs(x = "Configuration", y = "Mean percentage correct") +
  scale_y_continuous(label = function(x) {return(paste(x*100, "%"))}) +
  scale_x_discrete(limits=rev) +
  theme_ggstatsplot() +
  theme(legend.position = "none",
        axis.title = element_text(color = "darkred", size = 14)) +
  scale_fill_brewer(palette="Dark2")

#Response Time Analysis
  #One-way ANOVA for correctness
summary(aov(Reaction.Time ~ correct, data = nerrdata))

nerrdata |>
  ggplot(aes(x = conf, y = Reaction.Time)) +
  geom_boxplot(mapping = aes(fill = correct), alpha =0.5)+
  labs(x = "Configuration", y = "Mean response time", fill = "Response") +
  scale_y_log10(label = function(x) {return(paste(x, "ms"))}) +
  scale_x_discrete(limits=rev) +
  theme_ggstatsplot() +
  theme(axis.title = element_text(color = "darkred", size = 14)) +
  scale_fill_brewer(palette="Dark2", labels = c("Incorrect", "Correct"))

#Removing incorrect results from dataset
corrdata <- filter(nerrdata, correct == TRUE)

#One-way ANOVA for speaker
n.data <- corrdata |>
  filter(speaker %in% c("e","p"))
summary(aov(Reaction.Time ~ speaker, data = n.data))

d.data <- corrdata |>
  filter(speaker %in% c("m","s"))
summary(aov(Reaction.Time ~ speaker, data = d.data))

n.data |>
  ggplot(aes(x = speaker, y = Reaction.Time, fill = speaker)) +
  geom_boxplot(alpha = 0.5) +
  labs(x = "Speaker", y = "Response Time") +
  scale_x_discrete(limits=rev) +
  scale_y_log10(label = function(x) {return(paste(x, "ms"))}) +
  theme_ggstatsplot() +
  theme(legend.position = "none",
        axis.title = element_text(color = "darkred", size = 14)) +
  scale_fill_brewer(palette="Dark2")

d.data |>
  ggplot(aes(x = speaker, y = Reaction.Time, fill = speaker)) +
  geom_boxplot(alpha = 0.5) +
  labs(x = "Speaker", y = "Response Time") +
  scale_x_discrete(limits=rev) +
  scale_y_log10(label = function(x) {return(paste(x, "ms"))}) +
  theme_ggstatsplot() +
  theme(legend.position = "none",
        axis.title = element_text(color = "darkred", size = 14)) +
  scale_fill_brewer(palette="Dark2")


#Construction effect analysis
Anova(aov(Reaction.Time ~ construction, data = corrdata), type=2)
corrdata |>
  ggplot(aes(x = construction, y = Reaction.Time, fill = construction)) +
  geom_boxplot(alpha =0.5) +
  labs(x = "Construction", y = "Response time", fill = "Construction") +
  scale_y_log10(label = function(x) {return(paste(x, "ms"))}) +
  theme_ggstatsplot() +
  theme(axis.title = element_text(color = "darkred", size = 14),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())

p_const <- mp_const |> 
  group_by(construction)|> 
  summarize(correctness = mean(xp)) 
mrt_const <- corrdata |> 
  group_by(construction) |> 
  summarise(mean.reaction.time = mean(Reaction.Time)) 
const_data <- inner_join(p_const, mrt_const) 

const_data|>
  arrange(correctness)
const_data|>
  arrange(desc(mean.reaction.time))

#Linear mixed-effect regression
library(lme4)
library(merTools)
library(lmerTest)
library(MuMIn)

corr_model <- lmer(Reaction.Time ~ 
                     phonology + morphology + (1 | id) + (1 | construction / index), 
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
                       phonology + morphology + (1 | id) + (1 | construction / index), 
                     data = corrdata)
summary(rcorr_model)
fixef(rcorr_model)
plot(rcorr_model)
getR2(rcorr_model)

#Log transformed data
corr_log_model <- lmer(log(Reaction.Time) ~ 
                         phonology + morphology + (1 | id) + (1 | construction / index), 
                       data = corrdata)
summary(corr_log_model)
fixef(corr_log_model)
plot(corr_log_model)
getR2(corr_log_model)

ggqqplot(log(corrdata$Reaction.Time))

rcorr_log_model <- rlmer(log(Reaction.Time) ~ 
                           phonology + morphology + (1 | id) + (1 | construction / index), 
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
          title = "Coefficient estimates and 95% confidence intervals for log-linear models")

#Equivalence test
diff_corrdata <- corrdata |>
  mutate(pmsum = phonology + morphology,
         pmdif = phonology - morphology)
corr_log_diff_model <- lmer(log(Reaction.Time) ~ 
                              pmsum + pmdif + (1 | id) + (1 | construction / index), 
                            data = diff_corrdata )
summary(corr_log_diff_model)

corr_log_comp_model <- lmer(log(Reaction.Time) ~ 
                              I(phonology + morphology) + (1 | id) + (1 | construction / index), 
                            data = corrdata )
anova(corr_log_model, corr_log_comp_model)