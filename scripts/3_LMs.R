########################################################################
# LIBRARIES                                                            #
########################################################################

#TIDYING
library(xlsx)
library(tidyverse)

#PLOTTING
library(ggplot2)
library(sjPlot)
library(ggeffects)
library(ggfortify)
library(GGally)
library(gridExtra)
library(patchwork)
library(ggExtra)

#CORR
library(psych)

#MODEL SELECTION
library(AICcmodavg)
library(broom)


########################################################################
# DATA                                                                 #
########################################################################

# Pop info
popi <- read.csv("./data/climate/populations_locations_1.csv")
popc <- read.csv("./data/climate/climPCs.csv")

# Vernalization Experiment Means
datv <- read.csv("./data/phenotype/vern_2019_end_Durham_tidy2_means.csv")
colnames(datv)

colnames(datv) <- gsub("mean_", "StartFlowering_F1_", colnames(datv))
colnames(datv) <- gsub("sd_", "StartFlowering_sd_F1_", colnames(datv))
colnames(datv) <- gsub("se_", "StartFlowering_se_F1_", colnames(datv))

colnames(datv)

datv <- datv %>% select(-c(X, Control, Vernalization)) %>% as.data.frame()

# Data common garden -greenhouse grown F0
datp.germ <- read.csv("./data/phenotype/germination_4_2.csv", header = TRUE, stringsAsFactors = FALSE)
datp <- read.csv("./data/phenotype/gh_commonGarden_portsmouth.csv", header = TRUE, stringsAsFactors = FALSE)

datp <- datp %>%
  filter(is.na(G) | G != 2) %>% #keep only data from main greenhouse
  droplevels() %>% 
  mutate_at(c("Collected.on", "Sowdate", "GermDate", "FlFrDate"), as.character) %>%
  mutate_at(c("Collected.on", "Sowdate", "GermDate", "FlFrDate"), as.Date, format = "%Y-%m-%d") %>%
  filter(!(Pop %in% c("16", "Eden"))) %>% #this population actually corresponds to population 19, but seeds were collected when they were still immature, and just planted as a trial
  filter(FlFr != "NA") %>% #I filter what I don't know if it is flower or fruit record
  # distinct(Pop, Pop_Ind, G, T, P) %>% just checking that each family was replicated only once
  # group_by(Pop, Pop_Ind) %>%
  # summarise(n = n()) %>%
  # ungroup() %>% View()
  # group_by(Pop, Pop_Ind, FlFrDate) %>% # just to check that for each family and date of data collection, I have max 2 obs (1 for each flower and/or fruit)
  # mutate(n = n()) %>%
  # ungroup() %>%
  # filter(n >2) %>%
  # as.data.frame() %>% View() 
  distinct(FlFr, Pop, Pop_Ind, FlFrDate, .keep_all = TRUE) %>% #there were just 30 entries in total with n > 3, possibly due to mistakes during collection of data. I pick a random entry by using distinct()
  arrange(FlFrDate) %>%
  group_by(Pop, Pop_Ind, FlFr) %>% #NB: most families (Pop_Ind) were replicated only once (one pot), but there were a couple for which a few pots were present. I simply summarise data at family level.
  mutate(cumN = cumsum(no_FlFr)) %>%
  summarise(.groups = "keep", #summarise over collection date by:
            dataset = unique(dataset),
            experiment = unique(experiment),
            Country = unique(Country),
            Survival = unique(Survival),
            StartF = as.numeric(min(FlFrDate) - Sowdate[which.min(FlFrDate)]), #extracting earliest flowering time scored
            TotF_1mo = cumN[which.min(StartF >= StartF + 30)], #not sure why I calculated this
            TotF = sum(no_FlFr),
            MaxF = max(no_FlFr),
            Height = unique(Height),
            no_stems = unique(no_stem)
  ) %>% #View()
  ungroup() %>%
  droplevels() %>%
  as.data.frame()

datp.fs <- datp %>%
  select(c(Pop, Pop_Ind, FlFr, MaxF, TotF, StartF)) %>%
  filter(FlFr == "flowers") %>% #select only data related to flowering and exclude fruiting
  droplevels() %>%
  mutate(FlFr = paste0("FlFr", FlFr)) %>% #View()
  tidyr::pivot_wider(names_from = FlFr, values_from = c("MaxF", "TotF", "StartF")) %>% #View() #there are 5 families for which I have flowering obs, but survival is 0 or NA > I set survival to 1
  as.data.frame()

str(datp.fs)
colnames(datp.fs) <- c("Pop", "Pop_Ind", "Max_Flowers_F0", "Tot_Flowers_F0", "StartFlowering_F0")

datp.fs <- left_join(datp.germ %>% # fix pot survival with curated germination & survival file
                     filter(!(Pop %in% c("16", "Eden"))) %>%
                       filter(!(is.na(Pop))) %>%
                       droplevels() %>%
                       rename(Pop_Ind = Ind) %>%
                       as.data.frame(),
                     datp.fs,
                     by = c("Pop", "Pop_Ind")) %>%  #View()
  mutate(plantsXpot = if_else(is.na(plantsXpot) == TRUE, as.integer(0), plantsXpot)) %>% 
  mutate(plantsXpot = case_when(plantsXpot == as.integer(0) & StartFlowering_F0 > 0 ~ as.integer(1),
                                TRUE ~ plantsXpot)) %>% 
  mutate(flnfl = case_when(plantsXpot == as.integer(0) ~ "not_germ", #assign status to pot: no_germ (no seedling emerged at all); not_flowered (at least one seedling survived to flowering, but did not flower); flowered (at least one seedling survived to flowering and flowered)
                           plantsXpot != as.integer(0) & is.na(StartFlowering_F0) == TRUE ~ "not_flowered",
                           plantsXpot != as.integer(0) & is.na(StartFlowering_F0) == FALSE ~ "flowered")) %>%
  select(Pop, Pop_Ind, SowingDate:EmergDate, plantsXpot, flnfl, StartFlowering_F0, Max_Flowers_F0, Tot_Flowers_F0) %>%
  rename(Survival_Flowering = plantsXpot) %>%
  as.data.frame()

datp.fnf <- datp.fs %>%
  group_by(Pop) %>%
  mutate(
    fams.planted = n()
    ) %>%
  ungroup() %>%
  filter(flnfl != "not_germ") %>% #eliminate pots (=families) where there was no germination at all
  droplevels() %>%
  group_by(Pop) %>%
  mutate(
    fams.alive = n()) %>%
  ungroup() %>%
  group_by(Pop, flnfl) %>%
  summarise(
    fams.planted = unique(fams.planted),
    fams.alive = unique(fams.alive),
    fams.flnfl = n()) %>% 
  ungroup() %>% 
  spread(flnfl, fams.flnfl) %>%
  mutate(flowered = if_else(is.na(flowered) == TRUE, as.integer(0), flowered),
         not_flowered = if_else(is.na(not_flowered) == TRUE, as.integer(0), not_flowered)) %>%
  rename(fams.flowered_F0 = flowered,
         fams.notflowered_F0 = not_flowered,
         fams.planted_F0 = fams.planted,
         fams.alive_F0 = fams.alive) %>%
  droplevels() %>%
  as.data.frame()

View(datp.fnf)

datp.f.summ <- left_join(datp.fnf,
                         datp.fs %>%
                             droplevels() %>%
                             group_by(Pop) %>%
                             summarise(StartFlowering_F0_sd = as.integer(sd(StartFlowering_F0, na.rm = TRUE)),
                                       StartFlowering_F0 = as.integer(mean(StartFlowering_F0, na.rm = TRUE))) %>%
                           as.data.frame(),
                         by = "Pop") %>%
  as.data.frame()


# Final dataset to analyse with all infos
dat.f.summ <- datp.f.summ %>%
  mutate(species = "bienne") %>%
  full_join(.,
            datv,
            by = c("species", "Pop")) %>% 
  left_join(.,
            popi %>%
              mutate(Alt = if_else(is.na(Alt), elev_wc, Alt)) %>%
              select(c(Pop, Collected.on, Origin, Country, Lat, Lon, Alt, Distance.from.Coast)) %>%
              as.data.frame(),
              by = "Pop") %>%
  left_join(.,
            popc %>%
              select(-c(X)),
            by = "Pop") %>% #str()
  select(c(species, Pop, Collected.on:PC3, fams.planted_F0:StartFlowering_F0, StartFlowering_F1_Control:sensitivity.abs)) %>%
  as.data.frame() #%>% View()

# write.csv(dat.f.summ,
#           "./output/exp_summary_pop.csv")


########################################################################
# ANALYSIS                                                             #
########################################################################

# Keep L. bienne only
dat.f.summ <- droplevels(subset(dat.f.summ, species == "bienne"))


# Correlations between variables
## Flowering Onset in different experiments + Lat and climatic variables
dat.s <- dat.f.summ %>%
  select(Lat, PC1:PC3, StartFlowering_F0, StartFlowering_F1_Control, StartFlowering_F1_Vernalization, sensitivity)

ct <- corr.test(dat.s)

# write.xlsx(as.data.frame(ct$r),
#            "./output/var_corr.xlsx",
#            sheetName = "r",
#            row.names = TRUE,
#            col.names = TRUE,
#            append = FALSE)
# 
# write.xlsx(as.data.frame(ct$n),
#            "./output/var_corr.xlsx",
#            sheetName = "n",
#            row.names = TRUE,
#            col.names = TRUE,
#            append = TRUE)
# 
# write.xlsx(as.data.frame(ct$p),
#            "./output/var_corr.xlsx",
#            sheetName = "p",
#            row.names = TRUE,
#            col.names = TRUE,
#            append = TRUE)

remove(dat.s)


# F0 proportion of plants flowered
flnfl.mod <- glm(cbind(fams.flowered_F0, fams.alive_F0 - fams.flowered_F0) ~ Lat,
                 data = subset(dat.f.summ, is.na(fams.flowered_F0)==FALSE),
                 family = binomial(link = "logit"),
                 na.action = na.omit)

summary(flnfl.mod)
tidy(flnfl.mod)
glance(flnfl.mod)

car::Anova(flnfl.mod, type = "II")
hist(flnfl.mod$residuals, breaks = 5)

## Save model output
flnfl.summary <- as.data.frame(tidy(flnfl.mod)) %>%
  bind_cols(.,
            as.data.frame(glance(flnfl.mod))) %>%
  mutate(responsevar = "Flowering (Y/N)",
         model = "glm, binomial, logit") %>%
  as.data.frame()

flnfl.aov <- as.data.frame(car::Anova(flnfl.mod, type = "II")) %>%
  mutate(term = rownames(.),
         responsevar = "Flowering (Y/N)",
         model = "glm, binomial, logit") %>%
  as.data.frame()

# write.xlsx(flnfl.summary,
#            "./output/model_floweringProbability.xlsx",
#            sheetName = "summary",
#            row.names = TRUE,
#            col.names = TRUE,
#            append = FALSE)
# 
# write.xlsx(flnfl.aov,
#            "./output/model_floweringProbability.xlsx",
#            sheetName = "Anova",
#            row.names = TRUE,
#            col.names = TRUE,
#            append = TRUE)


## Make & Save plot
preds.flp <- ggemmeans(flnfl.mod, terms = c("Lat [all]")) %>%
  as.data.frame() %>% #str()
  rename(Lat = x,
         predicted.flp = predicted) %>% 
  dplyr::select(Lat, predicted.flp, conf.low, conf.high, std.error) %>%
  as.data.frame() #%>% str()

p.lbm <- subset(dat.f.summ, is.na(fams.flowered_F0)==FALSE) %>% 
  dplyr::select(Pop, Lat, fams.flowered_F0, fams.alive_F0) %>%
  mutate(fl.perc = fams.flowered_F0/fams.alive_F0) %>%
  as.data.frame() %>% #str()
  ggplot(aes(x = Lat, col = Lat)) +
  geom_point(aes(y = fl.perc), size = 2) +
  geom_smooth(data = preds.flp, aes(y = predicted.flp), size = 1, col = "grey38") +
  scale_color_gradient2(low = "navyblue", mid = "mediumpurple", high = "thistle", midpoint = 45) +
  labs(x = "Latitude", y = "Probability of Flowering", col = "Latitude") +
  theme_bw(base_size = 18) +
  theme(legend.position = "top",
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        strip.background = element_rect(fill = "white", colour = "white"),
        panel.border = element_blank())

# ggsave("./output/commonExp_floweringProbability.pdf", p.lbm, width = 8, height = 8, units = "in")


# Linear models for F0 and F1 flowering time and sensitivity to vernalization

## F0
### Model Selection
Cand.mod.fls0 <- list()
Cand.mod.fls0[[1]] <- lm(StartFlowering_F0 ~ PC1 * Lat,  data = subset(dat.f.summ, is.na(StartFlowering_F0)==FALSE))
Cand.mod.fls0[[2]] <- lm(StartFlowering_F0 ~ PC1 + Lat,  data = subset(dat.f.summ, is.na(StartFlowering_F0)==FALSE))
Cand.mod.fls0[[3]] <- lm(StartFlowering_F0 ~ PC1,  data = subset(dat.f.summ, is.na(StartFlowering_F0)==FALSE))
Cand.mod.fls0[[4]] <- lm(StartFlowering_F0 ~ Lat, data = subset(dat.f.summ, is.na(StartFlowering_F0)==FALSE))
Cand.mod.fls0[[5]] <- lm(StartFlowering_F0 ~ 1, data = subset(dat.f.summ, is.na(StartFlowering_F0)==FALSE))

mod.fls0names <- c("PC1 * Lat",
                   "PC1 + Lat",
                   "PC1",
                   "Lat",
                   "intercept only")

aic.tab.fls0 <- data.frame(aictab(cand.set = Cand.mod.fls0, modnames = mod.fls0names)) %>%
  mutate(responsevar = "flowering onset")

### Save model outputs for best model (PC1) and latitude model
#### Latitude
fls0.mod.lat <- Cand.mod.fls0[[4]]

summary(fls0.mod.lat)
tidy(fls0.mod.lat)
glance(fls0.mod.lat)

car::Anova(fls0.mod.lat, type = "II")
hist(fls0.mod.lat$residuals, breaks = 5)

#### PC1
fls0.mod.pc1 <- Cand.mod.fls0[[3]]

summary(fls0.mod.pc1)
tidy(fls0.mod.pc1)
glance(fls0.mod.pc1)

car::Anova(fls0.mod.pc1, type = "II")
hist(fls0.mod.pc1$residuals, breaks = 5)

#### Outputs
fls0.summary <- as.data.frame(tidy(fls0.mod.lat)) %>%
  bind_cols(.,
            as.data.frame(glance(fls0.mod.lat))) %>%
  mutate(responsevar = "Flowering Onset",
         predictor = "latitude",
         model = "lm, normal",
         generation = "F0") %>%
  bind_rows(.,
            as.data.frame(tidy(fls0.mod.pc1)) %>%
              bind_cols(.,
                        as.data.frame(glance(fls0.mod.pc1))) %>%
              mutate(responsevar = "Flowering Onset",
                     predictor = "pc1",
                     model = "lm, normal",
                     generation = "F0")) %>%
  as.data.frame()

fls0.aov <- as.data.frame(car::Anova(fls0.mod.lat, type = "II")) %>%
  mutate(term = rownames(.),
         responsevar = "Flowering Onset",
         predictor = "latitude",
         model = "lm, normal",
         generation = "F0") %>%
  bind_rows(.,
            as.data.frame(car::Anova(fls0.mod.pc1, type = "II")) %>%
              mutate(term = rownames(.),
                     responsevar = "Flowering Onset",
                     predictor = "pc1",
                     model = "lm, normal",
                     generation = "F0")) %>%
  as.data.frame()

preds.fls0 <- ggemmeans(fls0.mod.lat, terms = c("Lat [all]")) %>%
  as.data.frame() %>% #str()
  rename(Lat = x,
         predicted.fls0 = predicted) %>% 
  dplyr::select(Lat, predicted.fls0, conf.low, conf.high, std.error) %>%
  as.data.frame() #%>% str()

p.fls0 <- subset(dat.f.summ, is.na(StartFlowering_F0)==FALSE) %>% 
  dplyr::select(Pop, Lat, StartFlowering_F0) %>%
  as.data.frame() %>% #str()
  ggplot(aes(x = Lat, col = Lat)) +
  geom_point(aes(y = StartFlowering_F0), size = 2) +
  geom_smooth(data = preds.fls0, aes(y = predicted.fls0), size = 1, col = "grey38") +
  scale_color_gradient2(low = "navyblue", mid = "mediumpurple", high = "thistle", midpoint = 45) +
  labs(x = "Latitude", y = "Flowering Onset", col = "Latitude") +
  theme_bw(base_size = 18) +
  theme(legend.position = "top",
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        strip.background = element_rect(fill = "white", colour = "white"),
        panel.border = element_blank())

# ggsave("./output/commonExp_floweringOnset_f0.pdf", p.fls0, width = 8, height = 8, units = "in")

pf0 <- p.lbm + p.fls0 + plot_layout(guides = 'collect')
# ggsave("./output/commonExp_f0.pdf", pf0, width = 19, height = 8, units = "in")


## F1 - Flowering Onset
### Model Selection
Cand.mod.fls1 <- list()
Cand.mod.fls1[[1]] <- lm(StartFlowering_F1_Control ~ PC1 * Lat,  data = subset(dat.f.summ, is.na(StartFlowering_F1_Control)==FALSE))
Cand.mod.fls1[[2]] <- lm(StartFlowering_F1_Control ~ PC1 + Lat,  data = subset(dat.f.summ, is.na(StartFlowering_F1_Control)==FALSE))
Cand.mod.fls1[[3]] <- lm(StartFlowering_F1_Control ~ PC1,  data = subset(dat.f.summ, is.na(StartFlowering_F1_Control)==FALSE))
Cand.mod.fls1[[4]] <- lm(StartFlowering_F1_Control ~ Lat, data = subset(dat.f.summ, is.na(StartFlowering_F1_Control)==FALSE))
Cand.mod.fls1[[5]] <- lm(StartFlowering_F1_Control ~ 1, data = subset(dat.f.summ, is.na(StartFlowering_F1_Control)==FALSE))

mod.fls1names <- c("PC1 * Lat",
                   "PC1 + Lat",
                   "PC1",
                   "Lat",
                   "intercept only")

aic.tab.fls1 <- data.frame(aictab(cand.set = Cand.mod.fls1, modnames = mod.fls1names)) %>%
  mutate(responsevar = "flowering onset - control")

### Save model outputs for best model (PC1) and latitude model
#### Latitude
fls1.mod.lat <- Cand.mod.fls1[[4]]

summary(fls1.mod.lat)
tidy(fls1.mod.lat)
glance(fls1.mod.lat)

car::Anova(fls1.mod.lat, type = "II")
hist(fls1.mod.lat$residuals, breaks = 5)

#### PC1
fls1.mod.pc1 <- Cand.mod.fls1[[3]]

summary(fls1.mod.pc1)
tidy(fls1.mod.pc1)
glance(fls1.mod.pc1)

car::Anova(fls1.mod.pc1, type = "II")
hist(fls1.mod.pc1$residuals, breaks = 5)

#### Outputs
fls1.summary <- as.data.frame(tidy(fls1.mod.lat)) %>%
  bind_cols(.,
            as.data.frame(glance(fls1.mod.lat))) %>%
  mutate(responsevar = "Flowering Onset",
         predictor = "latitude",
         model = "lm, normal",
         generation = "f1") %>%
  bind_rows(.,
            as.data.frame(tidy(fls1.mod.pc1)) %>%
              bind_cols(.,
                        as.data.frame(glance(fls1.mod.pc1))) %>%
              mutate(responsevar = "Flowering Onset",
                     predictor = "pc1",
                     model = "lm, normal",
                     generation = "f1")) %>%
  as.data.frame()

fls1.aov <- as.data.frame(car::Anova(fls1.mod.lat, type = "II")) %>%
  mutate(term = rownames(.),
         responsevar = "Flowering Onset",
         predictor = "latitude",
         model = "lm, normal",
         generation = "f1") %>%
  bind_rows(.,
            as.data.frame(car::Anova(fls1.mod.pc1, type = "II")) %>%
              mutate(term = rownames(.),
                     responsevar = "Flowering Onset",
                     predictor = "pc1",
                     model = "lm, normal",
                     generation = "f1")) %>%
  as.data.frame()

preds.fls1 <- ggemmeans(fls1.mod.lat, terms = c("Lat [all]")) %>%
  as.data.frame() %>% #str()
  rename(Lat = x,
         predicted.fls1 = predicted) %>% 
  dplyr::select(Lat, predicted.fls1, conf.low, conf.high, std.error) %>%
  as.data.frame() #%>% str()

p.fls1 <- subset(dat.f.summ, is.na(StartFlowering_F1_Control)==FALSE) %>% 
  dplyr::select(Pop, Lat, StartFlowering_F1_Control) %>%
  as.data.frame() %>% #str()
  ggplot(aes(x = Lat, col = Lat)) +
  geom_point(aes(y = StartFlowering_F1_Control), size = 2) +
  geom_smooth(data = preds.fls1, aes(y = predicted.fls1), size = 1, col = "grey38") +
  scale_color_gradient2(low = "navyblue", mid = "mediumpurple", high = "thistle", midpoint = 45) +
  labs(x = "Latitude", y = "Flowering Onset", col = "Latitude") +
  theme_bw(base_size = 18) +
  theme(legend.position = "top",
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        strip.background = element_rect(fill = "white", colour = "white"),
        panel.border = element_blank())

# ggsave("./output/commonExp_floweringOnset_f1.pdf", p.fls1, width = 8, height = 8, units = "in")


## F1 - Vernalization Sensitivity
### Model Selection
Cand.mod.sns1 <- list()
Cand.mod.sns1[[1]] <- lm(sensitivity ~ PC1 * Lat,  data = subset(dat.f.summ, is.na(sensitivity)==FALSE))
Cand.mod.sns1[[2]] <- lm(sensitivity ~ PC1 + Lat,  data = subset(dat.f.summ, is.na(sensitivity)==FALSE))
Cand.mod.sns1[[3]] <- lm(sensitivity ~ PC1,  data = subset(dat.f.summ, is.na(sensitivity)==FALSE))
Cand.mod.sns1[[4]] <- lm(sensitivity ~ Lat, data = subset(dat.f.summ, is.na(sensitivity)==FALSE))
Cand.mod.sns1[[5]] <- lm(sensitivity ~ 1, data = subset(dat.f.summ, is.na(sensitivity)==FALSE))

mod.sns1names <- c("PC1 * Lat",
                   "PC1 + Lat",
                   "PC1",
                   "Lat",
                   "intercept only")

aic.tab.sns1 <- data.frame(aictab(cand.set = Cand.mod.sns1, modnames = mod.sns1names)) %>%
  mutate(responsevar = "vern. sensitivity")

### Save model outputs for best model (PC1) and latitude model
#### Latitude
sns1.mod.lat <- Cand.mod.sns1[[4]]

summary(sns1.mod.lat)
tidy(sns1.mod.lat)
glance(sns1.mod.lat)

car::Anova(sns1.mod.lat, type = "II")
hist(sns1.mod.lat$residuals, breaks = 5)

#### PC1
sns1.mod.pc1 <- Cand.mod.sns1[[3]]

summary(sns1.mod.pc1)
tidy(sns1.mod.pc1)
glance(sns1.mod.pc1)

car::Anova(sns1.mod.pc1, type = "II")
hist(sns1.mod.pc1$residuals, breaks = 5)

#### Outputs
sns1.summary <- as.data.frame(tidy(sns1.mod.lat)) %>%
  bind_cols(.,
            as.data.frame(glance(sns1.mod.lat))) %>%
  mutate(responsevar = "Vern. Sensitivity",
         predictor = "latitude",
         model = "lm, normal",
         generation = "f1") %>%
  bind_rows(.,
            as.data.frame(tidy(sns1.mod.pc1)) %>%
              bind_cols(.,
                        as.data.frame(glance(sns1.mod.pc1))) %>%
              mutate(responsevar = "Vern. Sensitivity",
                     predictor = "pc1",
                     model = "lm, normal",
                     generation = "f1")) %>%
  as.data.frame()

sns1.aov <- as.data.frame(car::Anova(sns1.mod.lat, type = "II")) %>%
  mutate(term = rownames(.),
         responsevar = "Vern. Sensitivity",
         predictor = "latitude",
         model = "lm, normal",
         generation = "f1") %>%
  bind_rows(.,
            as.data.frame(car::Anova(sns1.mod.pc1, type = "II")) %>%
              mutate(term = rownames(.),
                     responsevar = "Vern. Sensitivity",
                     predictor = "pc1",
                     model = "lm, normal",
                     generation = "f1")) %>%
  as.data.frame()

preds.sns1 <- ggemmeans(sns1.mod.lat, terms = c("Lat [all]")) %>%
  as.data.frame() %>% #str()
  rename(Lat = x,
         predicted.sns1 = predicted) %>% 
  dplyr::select(Lat, predicted.sns1, conf.low, conf.high, std.error) %>%
  as.data.frame() #%>% str()

p.sns1 <- subset(dat.f.summ, is.na(sensitivity)==FALSE) %>% 
  dplyr::select(Pop, Lat, sensitivity) %>%
  as.data.frame() %>% #str()
  ggplot(aes(x = Lat, col = Lat)) +
  geom_point(aes(y = sensitivity), size = 2) +
  geom_smooth(data = preds.sns1, aes(y = predicted.sns1), size = 1, col = "grey38") +
  scale_color_gradient2(low = "navyblue", mid = "mediumpurple", high = "thistle", midpoint = 45) +
  labs(x = "Latitude", y = "Vernalization Sensitivity", col = "Latitude") +
  theme_bw(base_size = 18) +
  theme(legend.position = "top",
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        strip.background = element_rect(fill = "white", colour = "white"),
        panel.border = element_blank())

# ggsave("./output/commonExp_vernSensitivity_f1.pdf", p.sns1, width = 8, height = 8, units = "in")

pf01 <- p.fls0 + p.sns1 + plot_layout(guides = 'collect')
# ggsave("./output/commonExp_f01.pdf", pf01, width = 19, height = 8, units = "in")

## Save all outputs
aic.tab <- bind_rows(aic.tab.fls0, aic.tab.fls1, aic.tab.sns1)
# write.xlsx(aic.tab,
#            "./output/model_floweringOnset.xlsx",
#            sheetName = "Model Selection",
#            row.names = TRUE,
#            col.names = TRUE,
#            append = FALSE)

aov.tab <- bind_rows(fls0.aov, fls1.aov, sns1.aov)
rownames(aov.tab) <- 1:nrow(aov.tab)
# write.xlsx(aov.tab,
#            "./output/model_floweringOnset.xlsx",
#            sheetName = "Anova",
#            row.names = TRUE,
#            col.names = TRUE,
#            append = TRUE)

summary.tab <- bind_rows(fls0.summary, fls1.summary, sns1.summary)
# write.xlsx(summary.tab,
#            "./output/model_floweringOnset.xlsx",
#            sheetName = "summary",
#            row.names = TRUE,
#            col.names = TRUE,
#            append = TRUE)

