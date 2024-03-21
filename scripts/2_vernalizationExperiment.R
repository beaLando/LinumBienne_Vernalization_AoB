########################################################################
# LIBRARIES                                                            #
########################################################################

library(tidyverse)
library(ggplot2)
library(ggeffects)
library(patchwork)

library(car)
library(lme4)
library(sjPlot)
library(sjtable2df)
library(xlsx)

library(bmbstats)


########################################################################
# DATA                                                                 #
########################################################################

dat <- read.csv("./data/phenotype/vern_2019_end_Durham_tidy2.csv") #flowering vernalization experiment
pop.info <- read.csv("./data/climate/populations_locations_1.csv")

# Add column "Fl_doy1": Flowering initiation considering Sowing Date
dat$Fl_doy1 <- as.numeric(as.Date(dat$flDate_1, format="%d/%m/%Y")-as.Date(dat$SowDate, format="%d/%m/%Y")) 
dat$Germ_doy <- as.numeric(as.Date(dat$germDate, format="%d/%m/%Y")-as.Date(dat$SowDate, format="%d/%m/%Y"))  

str(dat)

# Check sample sizes per specie, treatment, and climate chamber
dat %>%
  group_by(spp, treatment, SowDate) %>% #different climate chambers are represented by sowing date
  dplyr::summarise(n = n()) %>%
  ungroup() %>%
  as.data.frame() %>% View()

# Create datasets for each species individually
datb <- dat[dat$spp == "bienne",] # subset L. bienne
datu <- dat[dat$spp == "usitatissimum", ] # subset L. usitatissimum
str(datb)
str(datu)

# L. bienne
## Tidy up data
### At the end, L. bienne populations 7, 11, 12, 13, 14, and Sut were filtered out because represented in one treatment only (n < 2 in at least one of vernalization or control treatment)
### For L. usitatissimum, no observation was excluded

# datb.sum0 <- datb %>%
#   filter(is.na(SowDate) == FALSE) %>%
#   droplevels() %>%
#   group_by(Pop, treatment) %>%
#   mutate(sown.N = n(),
#          sown.N.fam = length(unique(Pop_Ind))) %>%
#   ungroup() %>%
#   filter(is.na(Germ_doy) == FALSE) %>%
#   droplevels() %>%
#   group_by(Pop, treatment) %>%
#   summarise(sown.N = unique(sown.N),
#             sown.N.fam = unique(sown.N.fam),
#             alive.N = n(),
#             alive.N.fam = length(unique(Pop_Ind))) %>%
#   ungroup() %>%
#   as.data.frame() #%>% View()
  

datb <- datb %>% 
  group_by(Pop, treatment) %>%
  dplyr::mutate(n.pt = sum(is.na(Fl_doy1)==FALSE)) %>%
  ungroup() %>% #View()
  filter(n.pt > 1) %>% #filter out pop X treatment combinations for which there was only one observation (1 flowering plant)
  droplevels() %>%
  group_by(Pop) %>%
  dplyr::mutate(sum.n.pt = sum(n.pt),
                sum.n.t = length(unique(treatment))) %>%
  ungroup() %>% #View()
  filter(sum.n.t > 1 & sum.n.pt > 3) %>% #keep pops that had flowering plants in both treatments and for which there were at least 4 observations overall
  droplevels() %>%
  select(-c(n.pt, sum.n.pt)) %>%
  as.data.frame() #%>% View()


## Calculate vernalization sensitivity
datb.ps <- datb %>%
  group_by(Pop, treatment) %>% 
  dplyr::summarise(ft.mean = mean(Fl_doy1, na.rm = TRUE)) %>% #str()
  ungroup() %>% #str()
  pivot_wider(names_from = treatment, 
              values_from = ft.mean,
              id_cols = c(Pop)) %>% 
  mutate(sensitivity = (Control - Vernalization)/(mean(Control, na.rm = TRUE) - mean(Vernalization, na.rm = TRUE))) %>%
  mutate(sensitivity.abs = abs(sensitivity)) %>%
  as.data.frame() #%>% View()

## Check sample sizes
datb %>%
  distinct(Pop, treatment, Pop_Ind) %>%
  group_by(Pop, treatment) %>%
  dplyr::summarise(n = n()) %>%
  ungroup() %>%
  group_by(treatment) %>%
  summarise(n = sum(n)) %>%
  ungroup() %>%
  as.data.frame() #%>% View()

datb %>%
  group_by(Pop, treatment) %>%
  dplyr::summarise(n = n()) %>%
  ungroup() %>%
  dplyr::mutate(n.avg = median(n)) %>%
  as.data.frame() #%>% View()

datb %>%
  group_by(treatment) %>%
  dplyr::summarise(mean.fl = mean(Fl_doy1),
                   sd = sd(Fl_doy1),
                   n = n()) %>%
  ungroup() %>%
  as.data.frame() #%>% View()

datb.ps %>% View()
mean(datb.ps$sensitivity)
sd(datb.ps$sensitivity)

datb.ps %>% View()
mean(datb.ps$sensitivity.abs)
sd(datb.ps$sensitivity.abs)


# L. usitatissimum
## Clean data (same as for L. bienne, but here all cultivars should be retained in the end)

# datu.sum0 <- datu %>%
#   filter(is.na(SowDate) == FALSE) %>%
#   droplevels() %>%
#   group_by(Pop, treatment) %>%
#   mutate(sown.N = n(),
#          sown.N.fam = length(unique(Pop_Ind))) %>%
#   ungroup() %>%
#   filter(is.na(Germ_doy) == FALSE) %>%
#   droplevels() %>%
#   group_by(Pop, treatment) %>%
#   summarise(sown.N = unique(sown.N),
#             sown.N.fam = unique(sown.N.fam),
#             alive.N = n(),
#             alive.N.fam = length(unique(Pop_Ind))) %>%
#   ungroup() %>%
#   as.data.frame() #%>% View()

datu <- datu %>% 
  group_by(Pop, treatment) %>%
  dplyr::mutate(n.pt = n()) %>%
  ungroup() %>% #View()
  filter(n.pt > 1) %>%
  droplevels() %>%
  group_by(Pop) %>%
  dplyr::mutate(sum.n.pt = sum(n.pt),
                sum.n.t = length(unique(treatment))) %>%
  ungroup() %>% #View()
  filter(sum.n.t > 1 & sum.n.pt > 3) %>%
  droplevels() %>%
  select(-c(n.pt, sum.n.pt)) %>%
  as.data.frame() #%>% View()

datu <- datu %>% #adding latitude info to Lus, because it is missing here
  left_join(.,
            pop.info %>%
              filter(spp == "usitatissimum") %>%
              droplevels() %>%
              select(c(Pop, Lat, Lon, elev_wc)),
            by = "Pop") %>%
  mutate(Lat = Lat.y,
         Lon = Lon.y) %>%
  select(-c(Lat.x, Lat.y, Lon.x, Lon.y)) %>%
  as.data.frame()


## Calculate vernalization sensitivity
datu.ps <- datu %>% 
  group_by(Pop, treatment) %>% 
  dplyr::summarise(ft.mean = mean(Fl_doy1, na.rm = TRUE)) %>% #str()
  ungroup() %>% #str()
  pivot_wider(names_from = treatment, 
              values_from = ft.mean,
              id_cols = c(Pop)) %>% 
  mutate(sensitivity = (Control - Vernalization)/(mean(Control, na.rm = TRUE) - mean(Vernalization, na.rm = TRUE))) %>%
  mutate(sensitivity.abs = abs(sensitivity)) %>%
  as.data.frame() #%>% View()

## Check sample sizes, etc.
datu %>%
  group_by(Pop, treatment) %>%
  dplyr::summarise(n = n()) %>%
  ungroup() %>%
  as.data.frame() #%>% View()

datu %>%
  group_by(treatment) %>%
  dplyr::summarise(mean.fl = mean(Fl_doy1),
                   sd = sd(Fl_doy1),
                   n = n()) %>%
  ungroup() %>%
  as.data.frame() #%>% View()

datu.ps %>% View()
mean(datu.ps$sensitivity)
sd(datu.ps$sensitivity)

datu.ps %>% View()
mean(datu.ps$sensitivity.abs)
sd(datu.ps$sensitivity.abs)

# Summarize data at population level
popb <- datb %>%
  group_by(Pop, treatment) %>%
  dplyr::summarise(N = n(),
                   N.fam = length(unique(Pop_Ind)),
                   mean = mean(Fl_doy1, na.rm = TRUE),
                   sd = sd(Fl_doy1, na.rm = TRUE),
                   se = sd/sqrt(N)) %>%
  # left_join(.,
  #           datb.sum0,
  #           by = c("Pop", "treatment")) %>%
  bind_rows(.,
            dat %>% 
              filter(Pop %in% c("7", "11", "12", "13", "14", "Sut")) %>%
              filter(is.na(Fl_doy1)==FALSE) %>%
              droplevels() %>%
              group_by(Pop, treatment) %>%
              dplyr::summarise(N = n(),
                               N.fam = length(unique(Pop_Ind))) %>%
              ungroup() %>%
              as.data.frame()) %>% #add back sample sizes for excluded populations 7, 11, 12, 13, 14, and Sut
  as.data.frame() #%>% View()

popu <- datu %>%
  group_by(Pop, treatment) %>%
  dplyr::summarise(N = n(),
                   N.fam = length(unique(Pop_Ind)),
                   mean = mean(Fl_doy1, na.rm = TRUE),
                   sd = sd(Fl_doy1, na.rm = TRUE),
                   se = sd/sqrt(N)) %>%
  # left_join(.,
  #           datu.sum0,
  #           by = c("Pop", "treatment")) %>%
  as.data.frame() #%>% View()

pop <- bind_rows(popb %>%
                   mutate(species = "bienne"), 
                 popu %>%
                   mutate(species = "usitatissimum")) %>% #str()
  pivot_wider(names_from = treatment, 
              values_from = c(mean, sd, se, N, N.fam),
              id_cols = c(species, Pop)) %>% #str()
  left_join(.,
            bind_rows(datb.ps %>%
                        dplyr::mutate(species = "bienne"),
                      datu.ps %>%
                        dplyr::mutate(species = "usitatissimum")),
            by = c("Pop", "species")) %>%
  mutate_all(as.character) %>%
  as.data.frame() #%>% View()

View(pop) #NAs in sample size columns (start with N_ or N.fam_) mean sample size was 0 (no plant flowering)

# write.csv(pop, "./data/phenotype/vern_2019_end_Durham_tidy2_means.csv")


########################################################################
# ANALYSIS                                                             #
########################################################################

# MIXED EFFECT MODELS - Flowering time response to vernalization

## L. bienne
lbm <- lmer(Fl_doy1 ~ no_stem + Pop*treatment + 
             (1|Pop:Pop_Ind) + (1|treatment:SowDate:block), 
           data = datb)

### Quick Checks & Model Summaries
car::Anova(lbm, type = "II")
drop1(lbm, test = "Chisq")
summary(lbm)

par(mfrow = c(2,4))
hist(resid(lbm), breaks = 200)
plot(resid(lbm) ~ fitted(lbm))
qqnorm(resid(lbm))
plot(resid(lbm, influence.measures = "rstandard") ~ lbm@frame$`no_stem`)
boxplot(resid(lbm, influence.measures = "rstandard") ~ interaction(lbm@frame$treatment, lbm@frame$Pop))
boxplot(resid(lbm, influence.measures = "rstandard") ~ lbm@frame$SowDate)
boxplot(resid(lbm, influence.measures = "rstandard") ~ lbm@frame$block)
boxplot(resid(lbm, influence.measures = "rstandard") ~ lbm@frame$Pop_Ind)

lbm.aov <- data.frame(car::Anova(lbm, type = "II")) %>% 
  mutate(species = "bienne",
         response = "days to flowering",
         predictors = rownames(.)) %>%
  as.data.frame()

rownames(lbm.aov) <- 1:4 

lbm.summary <- mtab2df(tab_model(lbm, 
                                 p.val = "kr", 
                                 show.df = TRUE), 
                       1, output = "data.frame") %>%
  mutate(species = "bienne",
         response = "days to flowering") %>%
  as.data.frame()


# write.xlsx(lbm.aov, 
#            "./output/vern_mixedEff_summary.xlsx", 
#            sheetName = "aovTab_bienne", 
#            row.names = FALSE, 
#            append = FALSE)
# 
# write.xlsx(lbm.summary, 
#            "./output/vern_mixedEff_summary.xlsx", 
#            sheetName = "summary_bienne", 
#            row.names = FALSE, 
#            append = TRUE)

### Model Plots
preds.lbm <- ggemmeans(lbm, terms = c("Pop", "treatment")) %>%
  as.data.frame() %>% #str()
  rename(Pop = x,
         predicted.fls = predicted,
         treatment = group) %>% 
  dplyr::select(predicted.fls, conf.low, conf.high, std.error, treatment, Pop) %>%
  left_join(.,
            datb %>%
              distinct(Pop, Lat),
            by = c("Pop")) %>%
  as.data.frame() #%>% str()

p.lbm <- datb %>% 
  dplyr::select(treatment, SowDate, block, Pop, Lat, Pop_Ind, Fl_doy1) %>%
  group_by(treatment, SowDate, block, Pop, Lat, Pop_Ind) %>%
  summarise(fls = mean(Fl_doy1, na.rm = TRUE)) %>%
  ungroup() %>%
  as.data.frame() %>% #str()
  ggplot(aes(x = treatment, col = Lat, fill = Lat)) +
  geom_point(aes(y = fls, group = Pop), size = 1.5, position = position_dodge(width = 0.2), alpha = 0.3) +
  geom_errorbar(data = preds.lbm, aes(ymin = conf.low, ymax = conf.high, group = Pop), size = 1, width = 0.3, position = position_dodge(width = 0.2)) +
  geom_line(data = preds.lbm, aes(y = predicted.fls, group = Pop), position = position_dodge(width = 0.2), size = 1, alpha = 0.5) +
  geom_point(data = preds.lbm, aes(y = predicted.fls, group = Pop), shape = 21, col = "black", size = 4, position = position_dodge(width = 0.2)) +
  scale_color_gradient2(low = "navyblue", mid = "mediumpurple", high = "thistle", midpoint = 45) +
  scale_fill_gradient2(low = "navyblue", mid = "mediumpurple", high = "thistle", midpoint = 45) +
  labs(title = "L. bienne", x = "Treatment", y = "Flowering Onset", fill = "Latitude", col = "Latitude") +
  scale_x_discrete(labels = c("Control" = "No Vernalization", "Vernalization" = "Vernalization")) +
  theme_bw(base_size = 18) +
  theme(legend.position = "none",
        plot.title = element_text(face = "italic"),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        strip.background = element_rect(fill = "white", colour = "white"),
        panel.border = element_blank())

## L. usitatissimum
lum <- lmer(Fl_doy1 ~ no_stem + Pop*treatment + 
              (1|treatment:SowDate:block), 
            data = datu)

### Quick Checks & Model Summaries
car::Anova(lum, type = "II")
drop1(lum, test = "Chisq")
summary(lum)

par(mfrow = c(2,4))
hist(resid(lum), breaks = 200)
plot(resid(lum) ~ fitted(lum))
qqnorm(resid(lum))
plot(resid(lum, influence.measures = "rstandard") ~ lum@frame$`no_stem`)
boxplot(resid(lum, influence.measures = "rstandard") ~ interaction(lum@frame$treatment, lum@frame$Pop))
boxplot(resid(lum, influence.measures = "rstandard") ~ lum@frame$SowDate)
boxplot(resid(lum, influence.measures = "rstandard") ~ lum@frame$block)

lum.aov <- data.frame(car::Anova(lum, type = "II")) %>% 
  mutate(species = "usitatissimum",
         response = "days to flowering",
         predictors = rownames(.)) %>%
  as.data.frame()

rownames(lum.aov) <- 1:4 

lum.summary <- mtab2df(tab_model(lum, 
                                 p.val = "kr", 
                                 show.df = TRUE), 
                       1, output = "data.frame") %>%
  mutate(species = "usitatissimum",
         response = "days to flowering") %>%
  as.data.frame()


# write.xlsx(lum.aov, 
#            "./output/vern_mixedEff_summary.xlsx", 
#            sheetName = "aovTab_usitatissimum", 
#            row.names = FALSE, 
#            append = TRUE)
# 
# write.xlsx(lum.summary, 
#            "./output/vern_mixedEff_summary.xlsx", 
#            sheetName = "summary_usitatissimum", 
#            row.names = FALSE, 
#            append = TRUE)

### Model Plots
preds.lum <- ggemmeans(lum, terms = c("Pop", "treatment")) %>%
  as.data.frame() %>% #str()
  rename(Pop = x,
         predicted.fls = predicted,
         treatment = group) %>% 
  dplyr::select(predicted.fls, conf.low, conf.high, std.error, treatment, Pop) %>%
  left_join(.,
            datu %>%
              distinct(Pop, Lat),
            by = c("Pop")) %>%
  as.data.frame() #%>% str()

p.lum <- datu %>% 
  dplyr::select(treatment, SowDate, block, Pop, Lat, Fl_doy1) %>%
  group_by(treatment, SowDate, block, Pop, Lat) %>%
  summarise(fls = mean(Fl_doy1, na.rm = TRUE)) %>%
  ungroup() %>%
  as.data.frame() %>% #str()
  ggplot(aes(x = treatment, col = Lat, fill = Lat)) +
  geom_point(aes(y = fls, group = Pop), size = 1.5, position = position_dodge(width = 0.2), alpha = 0.3) +
  geom_errorbar(data = preds.lum, aes(ymin = conf.low, ymax = conf.high, group = Pop), size = 1, width = 0.3, position = position_dodge(width = 0.2)) +
  geom_line(data = preds.lum, aes(y = predicted.fls, group = Pop), position = position_dodge(width = 0.2), size = 1, alpha = 0.5) +
  geom_point(data = preds.lum, aes(y = predicted.fls, group = Pop), shape = 21, col = "black", size = 4, position = position_dodge(width = 0.2)) +
  scale_color_gradient2(low = "navyblue", mid = "mediumpurple", high = "thistle", midpoint = 45) +
  scale_fill_gradient2(low = "navyblue", mid = "mediumpurple", high = "thistle", midpoint = 45) +
  labs(title = "L.usitatissimum", x = "Treatment", y = "Flowering Onset", fill = "Latitude", col = "Latitude") +
  scale_x_discrete(labels = c("Control" = "No Vernalization", "Vernalization" = "Vernalization")) +
  theme_bw(base_size = 18) +
  theme(legend.position = "top",
        plot.title = element_text(face = "italic"),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        strip.background = element_rect(fill = "white", colour = "white"),
        panel.border = element_blank())

p <- p.lbm + p.lum  +
  plot_layout(guides = 'collect')

# ggsave("./output/vernExp_estMeans.pdf", p, width = 20, height = 8, units = "in")



# SENSITIVITY COMPARISONS ACROSS SPECIES

vps <- datb.ps %>%
  dplyr::mutate(spp = "bienne") %>%
  bind_rows(.,
            datu.ps %>%
              dplyr::mutate(spp = "usitatissimum")) %>%
  as.data.frame()

sens.summary <- vps %>%
  group_by(spp) %>%
  dplyr::summarise(
    mean.abs = mean(sensitivity.abs),
    sd.abs = sd(sensitivity.abs),
    n = n()
  ) %>%
  as.data.frame()

describe_data(
  vps[vps$spp == "bienne",]$sensitivity.abs,
  estimator_function = data_estimators,
  control = model_control(),
  na.rm = FALSE)

describe_data(
  vps[vps$spp == "usitatissimum",]$sensitivity.abs,
  estimator_function = data_estimators,
  control = model_control(),
  na.rm = FALSE)

diff.boot <- compare_independent_groups( #function bootstraps conf.int. for differences between two groups for different metrics
  vps[vps$spp == "bienne",]$sensitivity.abs,
  vps[vps$spp == "usitatissimum",]$sensitivity.abs,
  SESOI_lower = -0.1,
  SESOI_upper = 0.1,
  estimator_function = independent_groups_estimators,
  control = model_control(),
  na.rm = FALSE)

diff.boot <- rbind(diff.boot)[[1]]

sens.summary <- sens.summary %>%
  bind_cols(.,
            diff.boot %>%
              filter(estimator %in% c("Mean diff", "SD diff"))) %>%
  as.data.frame()

## save
# write.xlsx(sens.summary, 
#            "./output/vern_mixedEff_summary.xlsx", 
#            sheetName = "summary_sensitivity", 
#            row.names = FALSE, 
#            append = TRUE)
