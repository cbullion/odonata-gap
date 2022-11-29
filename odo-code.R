###########################################################
#
# Analysis of Odonate Community Science Data
#
# Focus on comparisons between structured and unstructured
# sampling.
#
# GLMs for Richness and Abundance
# 
# NMDS analyses for composition between methods
#
# Shannon diversity
#
# CM Bullion; 17 May 2022; updated 28 November 2022
#
##########################################################





# INSTALL AND LOAD NEEDED PACKAGES

library(rgbif)
citation("rgbif")

library(plyr)
citation("plyr")

library(ggplot2)
citation("ggplot2")

library(vegan)
citation("vegan")

library(MASS)
citation("MASS")

library(broom)
citation("broom")

library(BiodiversityR)
citation("BiodiversityR")

library(betapart)
citation("betapart")

library(gridExtra)
citation("gridExtra")


### ---------- SECTION 1. SET UP ---------- ###

# Import raw data 
### GBIF data contains observations from iNaturalist and Odonata Central. 
### Field observation data contains structured observations at each county site.
setwd("~/GitHub/odonata-gap")
GBIF_OBS <- read.csv("citsci_obs_2019.csv", header=T)     # GBIF species occurrence.
FIELD_OBS <- read.csv("field_obs_2019.csv", header=T)     # Field observation data.

# Cleaning raw data.
### Subsetting both data sets.
GBIF_OBS <- subset(GBIF_OBS, select = c("species", "County", "month"))
FIELD_OBS <- subset(FIELD_OBS, select = c("county", "month", "scientific_name"))
### Adding Method type to both data sets. 'Unstructured' = Citizen science, 'Structured' = In-person
GBIF_OBS$Method="Unstructured"
FIELD_OBS$Method="Structured"
### Reorganize columns.
GBIF_OBS <- GBIF_OBS[, c(2, 3, 4, 1)]
FIELD_OBS <- FIELD_OBS[, c(1, 2, 4, 3)]
### Standardize column names. 
names(GBIF_OBS)[names(GBIF_OBS) == "County"] <- "county"
names(FIELD_OBS)[names(FIELD_OBS) == "scientific_name"] <- "species"
### Standardize site names. 
##### Citizen Science
GBIF_OBS$county[GBIF_OBS$county == 'Cuyahoga'] <- 'A'
GBIF_OBS$county[GBIF_OBS$county == 'Guernsey'] <- 'B'
GBIF_OBS$county[GBIF_OBS$county == 'Wayne'] <- 'C'
GBIF_OBS$county[GBIF_OBS$county == 'Knott'] <- 'D'
GBIF_OBS$county[GBIF_OBS$county == 'Wise'] <- 'E'
##### Structured
FIELD_OBS$county[FIELD_OBS$county == 'Cuyahoga'] <- 'A'
FIELD_OBS$county[FIELD_OBS$county == 'Guernsey'] <- 'B'
FIELD_OBS$county[FIELD_OBS$county == 'Wayne'] <- 'C'
FIELD_OBS$county[FIELD_OBS$county == 'Knott'] <- 'D'
FIELD_OBS$county[FIELD_OBS$county == 'Wise'] <- 'E'

# Create matrix
GBIF_COMP <- subset(GBIF_OBS, month == 6 | month == 7 | month == 8, select=c(county, Method, month, species))
GBIF_COMP <- subset(GBIF_OBS, select = c("species", "county", "Method"))
FIELD_COMP <- subset(FIELD_OBS, select = c("species", "county", "Method"))
### Aggregate dataset. 
GBIF_COMP <- count(GBIF_COMP, c("species", "county","Method"))
GBIF_COMP <- reshape(GBIF_COMP, idvar = c("county", "Method"), timevar = "species", direction = "wide")
FIELD_COMP <- count(FIELD_COMP, c("species", "county", "Method"))
FIELD_COMP <- reshape(FIELD_COMP, idvar = c("county", "Method"), timevar = "species", direction = "wide")
GBIF_COMP[3] <- NULL                  # Removing NULL columns.

# Combine datasets.
TOTAL_COMP <- rbind.fill(GBIF_COMP,FIELD_COMP)
TOTAL_COMP[is.na(TOTAL_COMP)] = 0     # Replacing NAs with 0s

# Change County and Source to factors
TOTAL_COMP <- within(TOTAL_COMP, {county <- factor(county)})
TOTAL_COMP <- within(TOTAL_COMP, {Method <- factor(Method)})

# Check data structure
str(TOTAL_COMP)
head(TOTAL_COMP)
tail(TOTAL_COMP)
dim(TOTAL_COMP)

# Column sums
GBIF_SUM <- GBIF_COMP[3:84]
FIELD_SUM <-FIELD_COMP[3:29]
colSums(GBIF_SUM, na.rm = TRUE)
colSums(FIELD_SUM, na.rm = TRUE)

# Set up color and point vectors for figures.
colvec2 <- c("#2d708eFF", "#c9c93fFF", "#95ac5aFF", "#618e74FF", "#2d708eFF")
pchvec1 <- c(16, 17, 15, 18)
pchvec2 <- c(19, 15, 4, 9)
ltyvec <- c(1, 2, 3, 4)

plottitle = ""

# Functions
make.sorted.plot <- function(x){ordered <- sort(x, T)
plot(
  ordered,
  col = colvec2,
  xaxt = "n", pch = 16, cex = 2,
  ylim = c(min(ordered)*0.5, max(ordered)),
  xlim = c(0, length(x)+1),
  ylab = "Diversity measure", xlab = "Sampling Location",
  main = plottitle)
text(ordered,
     names(ordered),
     srt = -75,
     pos = 4)
}


# Plot theme
theme_set(theme_bw())





### ---------- SECTION 2. SHANNON DIVERSITY ---------- ###

# Comparing richness by method.

# Rarefactions.
levels(TOTAL_COMP$Method)
RARE_SOURCE_S <- TOTAL_COMP[which(TOTAL_COMP$Method == "Structured"),]
RARE_SOURCE_U <- TOTAL_COMP[which(TOTAL_COMP$Method == "Unstructured"),]

SP_SOURCE_S <- specaccum(RARE_SOURCE_S[3:89], method = "rarefaction", permutations = 100, gamma = "jack2")
SP_SOURCE_U <- specaccum(RARE_SOURCE_U[3:89], method = "rarefaction", permutations = 100, gamma = "jack2")

plot(SP_SOURCE_S, pch = 19, col = "#FDE725FF", xvar = c("individuals"), lty = 1, lwd = 3,
     ylab = "Species Richness", xlab = "Number of Individuals", xlim = c(0, 1000), ylim = c(0, 100))
plot(SP_SOURCE_U, add = TRUE, pch = 15, xvar = c("individuals"), lty = 2, lwd = 3, col = "#2D708EFF")
legend("topleft", legend = c("Structured","Community Sci."), lty = ltyvec, bty = "n", lwd = 3, col = colvec2)

# Rarefaction plots
plot(SP_SOURCE_S)
data <- data.frame(Sites=SP_SOURCE_S$sites, Richness=SP_SOURCE_S$richness, SD=SP_SOURCE_S$sd)
RARE_SOURCE_S <- ggplot() +
  geom_point(data, mapping=aes(x=Sites, y=Richness)) +
  geom_line(data, mapping=aes(x=Sites, y=Richness)) +
  geom_ribbon(data, mapping=aes(x=Sites, ymin=(Richness-2*SD),ymax=(Richness+2*SD)),alpha=0.2) +
  theme_bw() +
  labs(title="Species Accumulation (Structured)")
  

plot(SP_SOURCE_U)
data2 <- data.frame(Sites=SP_SOURCE_U$sites, Richness=SP_SOURCE_U$richness, SD=SP_SOURCE_U$sd)
RARE_SOURCE_U <- ggplot() +
  geom_point(data2, mapping=aes(x=Sites, y=Richness)) +
  geom_line(data2, mapping=aes(x=Sites, y=Richness)) +
  geom_ribbon(data2, mapping=aes(x=Sites, ymin=(Richness-2*SD),ymax=(Richness+2*SD)),alpha=0.2) +
  theme_bw() +
  labs(title="Species Accumulation (Unstructured)")

grid.arrange(RARE_SOURCE_S, RARE_SOURCE_U, nrow = 1)


# Species richness for each site
###  A  B  C  D  E 
### 70 32 18 11 38 
specnumber(TOTAL_COMP[3:89], groups = TOTAL_COMP$county)

# Species richness for each method
### Unstructured Structured 
### 82 27
specnumber(TOTAL_COMP[3:89], groups = TOTAL_COMP$Method)

# Shannon Diversity
### By County
#####        A        B        C        D        E 
##### 3.456312 2.733695 2.530959 1.970978 3.099765 
SHANNON_COUNTY <- subset (TOTAL_COMP, select = -c(2))
SHANNON_COUNTY <- rowsum(SHANNON_COUNTY[,c(2:88)],SHANNON_COUNTY$county,na.rm=T)
SHANNON_COUNTY_RES <- diversity(SHANNON_COUNTY)
### By Method
#####       Unstructured       Structured 
##### 3.539527 2.864405 
SHANNON_METHOD <- subset (TOTAL_COMP, select = -c(1))
SHANNON_METHOD <- rowsum(SHANNON_METHOD[,c(2:88)],SHANNON_METHOD$Method,na.rm=T)
SHANNON_METHOD_RES <- diversity(SHANNON_METHOD) # Unstructured 3.539527 | Structured 2.864405 

plottitle = "Shannon diversity by location"
make.sorted.plot(SHANNON_COUNTY_RES)
plottitle = "Shannon diversity by sampling method"
make.sorted.plot(SHANNON_METHOD_RES)




### ---------- SECTION 3. RICHNESS ---------- ###

# Aggregate each dataset.
GBIF_RICH <- aggregate(species~county+Method+month, GBIF_OBS, FUN = function(x) length(unique(x)))
FIELD_RICH <- aggregate(species~county+Method+month, FIELD_OBS, FUN = function(x) length(unique(x)))

# Combine GBIF and Field Observation data into a single dataframe. 
TOT_RICH <- rbind.fill(GBIF_RICH,FIELD_RICH)
names(TOT_RICH)[names(TOT_RICH) == "species"] <- "richness"     # Rename species column to richness.

# Subset appropriate months
SUB_RICH <- subset(TOT_RICH, month == 6 | month == 7 | month == 8, select=c(county, Method, month, richness))

# Add missing combinations
SUB_RICH[nrow(SUB_RICH) + 1,] = c("A", "Structured", 6, 0)
SUB_RICH[nrow(SUB_RICH) + 1,] = c("B", "Structured", 6, 0)
SUB_RICH[nrow(SUB_RICH) + 1,] = c("C", "Unstructured", 6, 0)
SUB_RICH[nrow(SUB_RICH) + 1,] = c("C", "Unstructured", 7, 0)
SUB_RICH[nrow(SUB_RICH) + 1,] = c("C", "Unstructured", 8, 0)
SUB_RICH[nrow(SUB_RICH) + 1,] = c("D", "Unstructured", 8, 0)
SUB_RICH[nrow(SUB_RICH) + 1,] = c("E", "Structured", 8, 0)

SUB_RICH$Method <- as.factor(SUB_RICH$Method)
SUB_RICH$richness <- as.integer(SUB_RICH$richness)
SUB_RICH$month <- as.factor(SUB_RICH$month)
SUB_RICH$county <- as.factor(SUB_RICH$county)

# GLMs
# Without interaction.   AIC: 200.75
##### Coefficients:
#####              Estimate Std. Error z value Pr(>|z|)    
##### (Intercept)  3.418168   0.605702   5.643 1.67e-08 ***
##### sourceIP    -0.007808   0.441702  -0.018  0.98590    
##### countyB     -0.991657   0.676171  -1.467  0.14249    
##### countyC     -1.682924   0.690186  -2.438  0.01475 *  
##### countyD     -2.197006   0.708543  -3.101  0.00193 ** 
##### countyE     -1.083003   0.677532  -1.598  0.10994    
##### month7      -0.499256   0.541437  -0.922  0.35648    
##### month8      -0.250449   0.537217  -0.466  0.64107  
RICH_GLM1 <- glm.nb(formula = richness ~ Method + county + month, data = SUB_RICH)
tidy(RICH_GLM1)
summary(RICH_GLM1)

# With interaction.      AIC: 186.59
##### Coefficients:
#####                   Estimate  Std. Error z value Pr(>|z|)    
##### (Intercept)       3.874e+00  4.698e-01   8.245  < 2e-16 ***
##### sourceIP         -2.151e+00  6.372e-01  -3.375 0.000738 ***
##### countyB          -1.145e+00  6.019e-01  -1.902 0.057149 .  
##### countyC          -4.143e+01  3.875e+07   0.000 0.999999    
##### countyD          -3.471e+00  7.709e-01  -4.502 6.74e-06 ***
##### countyE          -1.474e+00  6.101e-01  -2.417 0.015646 *  
##### month7           -2.959e-01  3.849e-01  -0.769 0.442064    
##### month8           -1.411e-01  3.812e-01  -0.370 0.711353    
##### sourceIP:countyB  1.081e+00  9.131e-01   1.184 0.236448    
##### sourceIP:countyC  4.203e+01  3.875e+07   0.000 0.999999    
##### sourceIP:countyD  3.279e+00  1.037e+00   3.162 0.001568 ** 
##### sourceIP:countyE  1.824e+00  9.049e-01   2.015 0.043854 *  
RICH_GLM2 <- glm.nb(formula = richness ~ Method * county + month, data = SUB_RICH)
tidy(RICH_GLM2)
glance(RICH_GLM2)
summary(RICH_GLM2)

# Generate boxplots
### Separated
RICH_BOX1 <- ggplot(SUB_RICH, aes(x=as.factor(Method), y=as.numeric(richness), fill=Method)) + 
  geom_boxplot() +
  facet_wrap(~county, nrow = 1, scale="fixed") +
  labs(title="Species richness by method per county", x="Sampling Method", y="Species Richness") +
  theme_bw() +
  theme(legend.key.size = unit(2, 'cm')) +
  guides(color=guide_legend("Sampling Method")) +
  scale_fill_manual(values=c("#fde725","#2d708e"))
RICH_BOX1

# plot(RICH_BOX1)
### Combined
RICH_BOX2 <- ggplot(SUB_RICH, aes(x=county, y=as.numeric(richness), fill=Method)) + 
  geom_boxplot() +
  labs(title="Species richness by method per county", x="Sampling Method", y="Species Richness") +
  theme_bw() +
  theme(legend.key.size = unit(2, 'cm')) +
  guides(color=guide_legend("Sampling Method")) +
  scale_fill_manual(values=c("#fde725","#2d708e"))
plot(RICH_BOX2)

# T.Test
### t = 1.6925, df = 15.639, p-value = 0.1104
RICH_TT1 <- t.test(SUB_RICH$richness ~ SUB_RICH$Method)
RICH_TT1





### ---------- SECTION 4. ABUNDANCE ---------- ###

# Aggregate each dataset.
GBIF_ABUND <- count(GBIF_OBS, c("county", "month","Method"))
FIELD_ABUND <- count(FIELD_OBS, c("county", "month", "Method"))

# Combine GBIF and Field Observation data into a single dataframe.
TOT_AB <- rbind.fill(GBIF_ABUND,FIELD_ABUND)
names(TOT_AB)[names(TOT_AB) == "freq"] <- "abundance"

# Subset appropriate months
TOT_AB <- subset(TOT_AB, month == 6 | month == 7 | month == 8, select=c(county, Method, month, abundance))

# Add missing combinations
TOT_AB[nrow(TOT_AB) + 1,] = c("A", "Structured", 6, 0)
TOT_AB[nrow(TOT_AB) + 1,] = c("B", "Structured", 6, 0)
TOT_AB[nrow(TOT_AB) + 1,] = c("C", "Unstructured", 6, 0)
TOT_AB[nrow(TOT_AB) + 1,] = c("C", "Unstructured", 7, 0)
TOT_AB[nrow(TOT_AB) + 1,] = c("C", "Unstructured", 8, 0)
TOT_AB[nrow(TOT_AB) + 1,] = c("D", "Unstructured", 8, 0)
TOT_AB[nrow(TOT_AB) + 1,] = c("E", "Structured", 8, 0)

TOT_AB$Method <- as.factor(TOT_AB$Method)
TOT_AB$abundance <- as.integer(TOT_AB$abundance)
TOT_AB$county <- as.factor(TOT_AB$county)
TOT_AB$month <- as.factor(TOT_AB$month)

# GLM
# Without interaction.   AIC: 262.14
##### Coefficients:
#####                   Estimate  Std.Error z value Pr(>|z|)    
##### (Intercept)         5.2111     0.7910   6.588 4.45e-11 ***
##### sourceIP            1.1623     0.5655   2.055  0.03984 *  
##### countyB            -2.6766     0.8836  -3.029  0.00245 ** 
##### countyC            -3.5107     0.8903  -3.943 8.03e-05 ***
##### countyD            -3.9828     0.8972  -4.439 9.03e-06 ***
##### countyE            -2.5237     0.8829  -2.859  0.00425 ** 
##### as.factor(month)7  -0.3730     0.6942  -0.537  0.59107    
##### as.factor(month)8   0.2737     0.6903   0.396  0.69180 
AB_GLM1 <- glm.nb(formula = abundance ~ Method + county + as.factor(month), data = TOT_AB)
summary(AB_GLM1)

# With interaction.      AIC: 251.54
##### Coefficients:
#####                    Estimate Std. Error z value Pr(>|z|)    
##### (Intercept)       5.672e+00  7.261e-01   7.811 5.69e-15 ***
##### sourceIP         -2.529e+00  9.292e-01  -2.722 0.006488 ** 
##### countyB          -2.892e+00  9.323e-01  -3.102 0.001923 ** 
##### countyC          -4.300e+01  3.875e+07   0.000 0.999999    
##### countyD          -5.426e+00  1.040e+00  -5.220 1.79e-07 ***
##### countyE          -2.951e+00  9.328e-01  -3.163 0.001560 ** 
##### month7           -2.685e-02  5.527e-01  -0.049 0.961254    
##### month8            4.275e-01  5.496e-01   0.778 0.436605    
##### sourceIP:countyB  2.762e+00  1.321e+00   2.090 0.036591 *  
##### sourceIP:countyC  4.323e+01  3.875e+07   0.000 0.999999    
##### sourceIP:countyD  4.954e+00  1.401e+00   3.536 0.000407 ***
##### sourceIP:countyE  3.242e+00  1.320e+00   2.456 0.014055 *  
AB_GLM2 <- glm.nb(formula = abundance ~ Method * county + month, data = TOT_AB)
tidy(AB_GLM2)
summary(AB_GLM2)

# Generate boxplots
### Separated
AB_BOX1 <- ggplot(TOT_AB, aes(x=as.factor(Method), y=as.numeric(abundance), fill=Method)) + 
  geom_boxplot() +
  facet_wrap(~county, nrow = 1, scale="fixed") +
  labs(title="Odonate abundance by method per county", x="Sampling Method", y="Abundance") +
  theme(legend.key.size = unit(2, 'cm')) +
  scale_fill_manual(values=c("#fde725","#2d708e"))
plot(AB_BOX1)
### Combined
AB_BOX2 <- ggplot(TOT_AB, aes(x=county, y=as.numeric(abundance), fill=Method)) + 
  geom_boxplot() +
  scale_fill_manual(values=c("#fde725","#2d708e"))
plot(AB_BOX2)

# T.Test
### t = 1.1867, df = 15.166, p-value = 0.2536
ABUND_TT1 <- t.test(TOT_AB$abundance ~ TOT_AB$Method)
ABUND_TT1




### ---------- SECTION 5. COMPOSITION ---------- ###

# Subset dataset. 
GBIF_COMP <- subset(GBIF_OBS, month == 6 | month == 7 | month == 8, select=c(county, Method, month, species))
GBIF_COMP <- subset(GBIF_OBS, select = c("species", "county", "Method"))
FIELD_COMP <- subset(FIELD_OBS, select = c("species", "county", "Method"))

# Aggregate dataset. 
GBIF_COMP <- count(GBIF_COMP, c("species", "county","Method"))
GBIF_COMP <- reshape(GBIF_COMP, idvar = c("county", "Method"), timevar = "species", direction = "wide")
FIELD_COMP <- count(FIELD_COMP, c("species", "county", "Method"))
FIELD_COMP <- reshape(FIELD_COMP, idvar = c("county", "Method"), timevar = "species", direction = "wide")

# Removing NULL columns.
GBIF_COMP[3] <- NULL

# Combine datasets.
TOTAL_COMP <- rbind.fill(GBIF_COMP,FIELD_COMP)
TOTAL_COMP[is.na(TOTAL_COMP)] = 0     # Replacing NAs with 0s

# Create matrix
ENV.MATRIX <- TOTAL_COMP[c(1:2)]      # Environmental variables
COM.MATRIX <- TOTAL_COMP[c(3:89)]     # Community variables


# NMDS | Stress: 0.03141493 
NMDS <- metaMDS(COM.MATRIX, distance="bray", k=3, autotransform=FALSE, trymax=1000)
stressplot(NMDS)

# NMDS Visualization
plot(NMDS)
summary(NMDS)
# add ellipsoids with ordiellipse
ordiellipse(NMDS, ENV.MATRIX$Method, draw="polygon", col="#fde725",kind="sd", conf=0.95, label=FALSE, show.groups = "Unstructured")
ordiellipse(NMDS, ENV.MATRIX$Method, draw="polygon", col="#2d708e",kind="sd", conf=0.95, label=FALSE, show.groups = "Structured") 
# display sampling Method as shapes. Unstructured = Circle. Structured = Triangle. 
points(NMDS, display="sites", select=which(ENV.MATRIX$Method=="Unstructured"),pch=19, col="#FFC000")
points(NMDS, display="sites", select=which(ENV.MATRIX$Method=="Structured"), pch=17, col="#003976")
# legend
# legend(1.46,1.45, title=NULL, pch=c(19,17,15,25), col=c("#fde725","#2d708e"), cex=.7, legend=c("Citizen Science", "In-Person"))
legend(x = "topright", title=NULL, pch=c(19,17,15,25), col=c("#fde725","#2d708e"), cex=.7, legend=c("Citizen Science", "In-Person"))
?legend

# ANOSIM, p = 0.0679
ano = anosim(COM.MATRIX, ENV.MATRIX$Method, distance = "bray", permutations = 9999)
ano


# ADONIS, p = 0.07
AD_COMP<-adonis2(COM.MATRIX ~ Method, data = ENV.MATRIX, permutations = 999, method="bray")
AD_COMP







### ---------- SECTION 6. BETA DISPERSION ---------- ###

# Convert to presence absence.
TOTAL_COMP_2 <- TOTAL_COMP
TOTAL_COMP_2A <- TOTAL_COMP[2]
TOTAL_COMP_2B <- TOTAL_COMP_2[-1:-2]
TOTAL_COMP_2B[TOTAL_COMP_2B>0]<- 1

source_vec <- TOTAL_COMP_2A$Method 

TOTAL_COMP_2B["Method"] <- source_vec

# Create beta part object
COMP_CORE <- betapart.core(TOTAL_COMP_2B[1:87])

# Returns three dissimilarity matrices containing 
# Pairwise between-site values of each beta-diversity component
COMP_DIST <- beta.pair(COMP_CORE, index.family = "sorensen")
str(COMP_DIST)

# Beta dispersion
#####     Unstructured     Structured 
##### 0.4245 0.3239 
BETA_SOR <- betadisper(COMP_DIST$beta.sor, TOTAL_COMP_2B$Method, type = c("median"))
BETA_SOR

anova(BETA_SOR) # p = 0.4662
plot(BETA_SOR)
boxplot(BETA_SOR, ylab = "Distance to median")

# Tukey multiple comparisons of means
#####             diff        lwr       upr     p adj
##### Structured-Unstructured -0.1006139 -0.4094017 0.2081738 0.4662198
TukeyHSD(BETA_SOR, which = "group", conf.level = 0.95)
