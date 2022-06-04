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
# CM Bullion; 17 May 2022; updated 17 May 2022
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



### ---------- SECTION 1. SET UP ---------- ###

# Import raw data 
### GBIF data contains observations from iNaturalist and Odonata Central. 
### Field observation data contains structured observations at each county site.
setwd("~/Coding/GitHub/odonata-gap")
GBIF_OBS <- read.csv("citsci_obs_2019.csv", header=T)     # GBIF species occurrence.
FIELD_OBS <- read.csv("field_obs_2019.csv", header=T)     # Field observation data.

# Cleaning raw data.
### Subsetting both data sets.
GBIF_OBS <- subset(GBIF_OBS, select = c("species", "County", "month"))
FIELD_OBS <- subset(FIELD_OBS, select = c("county", "month", "scientific_name"))
### Adding source type to both data sets. 'CS' = Citizen science, 'IP' = In-person
GBIF_OBS$source="CS"
FIELD_OBS$source="IP"
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
GBIF_COMP <- subset(GBIF_OBS, month == 6 | month == 7 | month == 8, select=c(county, source, month, species))
GBIF_COMP <- subset(GBIF_OBS, select = c("species", "county", "source"))
FIELD_COMP <- subset(FIELD_OBS, select = c("species", "county", "source"))
### Aggregate dataset. 
GBIF_COMP <- count(GBIF_COMP, c("species", "county","source"))
GBIF_COMP <- reshape(GBIF_COMP, idvar = c("county", "source"), timevar = "species", direction = "wide")
FIELD_COMP <- count(FIELD_COMP, c("species", "county", "source"))
FIELD_COMP <- reshape(FIELD_COMP, idvar = c("county", "source"), timevar = "species", direction = "wide")
GBIF_COMP[3] <- NULL                  # Removing NULL columns.

# Combine datasets.
TOTAL_COMP <- rbind.fill(GBIF_COMP,FIELD_COMP)
TOTAL_COMP[is.na(TOTAL_COMP)] = 0     # Replacing NAs with 0s

# Change County and Source to factors
TOTAL_COMP <- within(TOTAL_COMP, {county <- factor(county)})
TOTAL_COMP <- within(TOTAL_COMP, {source <- factor(source)})

# Check data structure
str(TOTAL_COMP)
head(TOTAL_COMP)
tail(TOTAL_COMP)
dim(TOTAL_COMP)

# Set up color and point vectors for figures.
colvec1 <- c("darkorange1", "darkgoldenrod1", "black", "chartreuse4")
colvec2 <- c("#FDE725FF", "#2D708EFF", "#481567FF", "#73D055FF")
pchvec1 <- c(16, 17, 15, 18)
pchvec2 <- c(19, 15, 4, 9)
ltyvec <- c(1, 2, 3, 4)

# Functions
make.sorted.plot <- function(x){ordered <- sort(x, T)
plot(
  ordered,
  col = colvec1,
  xaxt = "n", pch = 16, cex = 2,
  ylim = c(min(ordered)*0.5, max(ordered)),
  xlim = c(0, length(x)+1),
  ylab = "Diversity measure", xlab = "Sites",
  main = substitute(x))
text(ordered,
     names(ordered),
     srt = -75,
     pos = 4)
}








### ---------- SECTION 2. SHANNON DIVERSITY ---------- ###

# Comparing richness by method.

# Rarefactions.
levels(TOTAL_COMP$source)
RARE_SOURCE_S <- TOTAL_COMP[which(TOTAL_COMP$source == "IP"),]
RARE_SOURCE_U <- TOTAL_COMP[which(TOTAL_COMP$source == "CS"),]

SP_SOURCE_S <- specaccum(RARE_SOURCE_S[3:89], method = "rarefaction", permutations = 100, gamma = "jack2")
SP_SOURCE_U <- specaccum(RARE_SOURCE_U[3:89], method = "rarefaction", permutations = 100, gamma = "jack2")

plot(SP_SOURCE_S, pch = 19, col = "#FDE725FF", xvar = c("individuals"), lty = 1, lwd = 3,
     ylab = "Species Richness", xlab = "Number of Individuals", xlim = c(0, 1000), ylim = c(0, 100))
plot(SP_SOURCE_U, add = TRUE, pch = 15, xvar = c("individuals"), lty = 2, lwd = 3, col = "#2D708EFF")
legend("topleft", legend = c("Structured","Community Sci."), lty = ltyvec, bty = "n", lwd = 3, col = colvec2)

# Species richness for each site
specnumber(TOTAL_COMP[3:89], groups = TOTAL_COMP$county)

# Species richness for each method
specnumber(TOTAL_COMP[3:89], groups = TOTAL_COMP$source)

# Shannon Diversity
### By County
SHANNON_COUNTY <- subset (TOTAL_COMP, select = -c(2))
SHANNON_COUNTY <- rowsum(SHANNON_COUNTY[,c(2:88)],SHANNON_COUNTY$county,na.rm=T)
SHANNON_COUNTY_RES <- diversity(SHANNON_COUNTY)
### By Method
SHANNON_METHOD <- subset (TOTAL_COMP, select = -c(1))
SHANNON_METHOD <- rowsum(SHANNON_METHOD[,c(2:88)],SHANNON_METHOD$source,na.rm=T)
SHANNON_METHOD_RES <- diversity(SHANNON_METHOD)
### Plots
make.sorted.plot(SHANNON_COUNTY_RES)
make.sorted.plot(SHANNON_METHOD_RES)





### ---------- SECTION 3. RICHNESS ---------- ###

# Aggregate each dataset.
GBIF_RICH <- aggregate(species~county+source+month, GBIF_OBS, FUN = function(x) length(unique(x)))
FIELD_RICH <- aggregate(species~county+source+month, FIELD_OBS, FUN = function(x) length(unique(x)))

# Combine GBIF and Field Observation data into a single dataframe. 
TOT_RICH <- rbind.fill(GBIF_RICH,FIELD_RICH)
names(TOT_RICH)[names(TOT_RICH) == "species"] <- "richness"     # Rename species column to richness.

# Subset appropriate months
SUB_RICH <- subset(TOT_RICH, month == 6 | month == 7 | month == 8, select=c(county, source, month, richness))

# Add missing combinations
SUB_RICH[nrow(SUB_RICH) + 1,] = c("A", "IP", 6, 0)
SUB_RICH[nrow(SUB_RICH) + 1,] = c("B", "IP", 6, 0)
SUB_RICH[nrow(SUB_RICH) + 1,] = c("C", "CS", 6, 0)
SUB_RICH[nrow(SUB_RICH) + 1,] = c("C", "CS", 7, 0)
SUB_RICH[nrow(SUB_RICH) + 1,] = c("C", "CS", 8, 0)
SUB_RICH[nrow(SUB_RICH) + 1,] = c("D", "CS", 8, 0)
SUB_RICH[nrow(SUB_RICH) + 1,] = c("E", "IP", 8, 0)

SUB_RICH$source <- as.factor(SUB_RICH$source)
SUB_RICH$richness <- as.integer(SUB_RICH$richness)
SUB_RICH$month <- as.factor(SUB_RICH$month)
SUB_RICH$county <- as.factor(SUB_RICH$county)

# GLMs
### Without interaction.   AIC: 200.75
RICH_GLM1 <- glm.nb(formula = richness ~ source + county + month, data = SUB_RICH)
tidy(RICH_GLM1)
summary(RICH_GLM1)
### With interaction.      AIC: 186.59
RICH_GLM2 <- glm.nb(formula = richness ~ source * county + month, data = SUB_RICH)
tidy(RICH_GLM2)
glance(RICH_GLM2)
summary(RICH_GLM2)

# Generate boxplots
### Separated
RICH_BOX1 <- ggplot(SUB_RICH, aes(x=as.factor(source), y=as.numeric(richness), fill=source)) + 
  geom_boxplot() +
  facet_wrap(~county, scale="free")
plot(RICH_BOX1)
### Combined
RICH_BOX2 <- ggplot(SUB_RICH, aes(x=county, y=as.numeric(richness), fill=source)) + 
  geom_boxplot()
plot(RICH_BOX2)

# T.Test
RICH_TT1 <- t.test(SUB_RICH$richness ~ SUB_RICH$source)
RICH_TT1





### ---------- SECTION 4. ABUNDANCE ---------- ###

# Aggregate each dataset.
GBIF_ABUND <- count(GBIF_OBS, c("county", "month","source"))
FIELD_ABUND <- count(FIELD_OBS, c("county", "month", "source"))

# Combine GBIF and Field Observation data into a single dataframe.
TOT_AB <- rbind.fill(GBIF_ABUND,FIELD_ABUND)
names(TOT_AB)[names(TOT_AB) == "freq"] <- "abundance"

# Subset appropriate months
TOT_AB <- subset(TOT_AB, month == 6 | month == 7 | month == 8, select=c(county, source, month, abundance))

# Add missing combinations
TOT_AB[nrow(TOT_AB) + 1,] = c("A", "IP", 6, 0)
TOT_AB[nrow(TOT_AB) + 1,] = c("B", "IP", 6, 0)
TOT_AB[nrow(TOT_AB) + 1,] = c("C", "CS", 6, 0)
TOT_AB[nrow(TOT_AB) + 1,] = c("C", "CS", 7, 0)
TOT_AB[nrow(TOT_AB) + 1,] = c("C", "CS", 8, 0)
TOT_AB[nrow(TOT_AB) + 1,] = c("D", "CS", 8, 0)
TOT_AB[nrow(TOT_AB) + 1,] = c("E", "IP", 8, 0)

TOT_AB$source <- as.factor(TOT_AB$source)
TOT_AB$abundance <- as.integer(TOT_AB$abundance)
TOT_AB$county <- as.factor(TOT_AB$county)
TOT_AB$month <- as.factor(TOT_AB$month)

# GLM
### Without interaction.   AIC: 262.14
AB_GLM1 <- glm.nb(formula = abundance ~ source + county + as.factor(month), data = TOT_AB)
summary(AB_GLM1)
### With interaction.      AIC: 251.54
AB_GLM2 <- glm.nb(formula = abundance ~ source * county + month, data = TOT_AB)
tidy(AB_GLM2)
summary(AB_GLM2)

# Generate boxplots
### Separated
AB_BOX1 <- ggplot(TOT_AB, aes(x=as.factor(source), y=as.numeric(abundance), fill=source)) + 
  geom_boxplot() +
  facet_wrap(~county, scale="free")
plot(AB_BOX1)
### Combined
AB_BOX2 <- ggplot(TOT_AB, aes(x=county, y=as.numeric(abundance), fill=source)) + 
  geom_boxplot()
plot(AB_BOX2)





### ---------- SECTION 5. COMPOSITION ---------- ###

# Subset dataset. 
GBIF_COMP <- subset(GBIF_OBS, month == 6 | month == 7 | month == 8, select=c(county, source, month, species))
GBIF_COMP <- subset(GBIF_OBS, select = c("species", "county", "source"))
FIELD_COMP <- subset(FIELD_OBS, select = c("species", "county", "source"))

# Aggregate dataset. 
GBIF_COMP <- count(GBIF_COMP, c("species", "county","source"))
GBIF_COMP <- reshape(GBIF_COMP, idvar = c("county", "source"), timevar = "species", direction = "wide")
FIELD_COMP <- count(FIELD_COMP, c("species", "county", "source"))
FIELD_COMP <- reshape(FIELD_COMP, idvar = c("county", "source"), timevar = "species", direction = "wide")

# Removing NULL columns.
GBIF_COMP[3] <- NULL

# Combine datasets.
TOTAL_COMP <- rbind.fill(GBIF_COMP,FIELD_COMP)
TOTAL_COMP[is.na(TOTAL_COMP)] = 0     # Replacing NAs with 0s

# Create matrix
ENV.MATRIX <- TOTAL_COMP[c(1:2)]      # Environmental variables
COM.MATRIX <- TOTAL_COMP[c(3:89)]     # Community variables


# NMDS
NMDS <- metaMDS(COM.MATRIX, distance="bray", k=3, autotransform=FALSE, trymax=1000)
stressplot(NMDS)

# NMDS Visualization
plot(NMDS)
summary(NMDS)
# add ellipsoids with ordiellipse
ordiellipse(NMDS, ENV.MATRIX$source, draw="polygon", col="#FFC000",kind="sd", conf=0.95, label=FALSE, show.groups = "CS")
ordiellipse(NMDS, ENV.MATRIX$source, draw="polygon", col="#003976",kind="sd", conf=0.95, label=FALSE, show.groups = "IP") 
# display sampling source as shapes. CS = Circle. IP = Triangle. 
points(NMDS, display="sites", select=which(ENV.MATRIX$source=="CS"),pch=19, col="#FFC000")
points(NMDS, display="sites", select=which(ENV.MATRIX$source=="IP"), pch=17, col="#003976")
# legend
legend(1.46,1.45, title=NULL, pch=c(19,17,15,25), col=c("#FFC000","#003976"), cex=.7, legend=c("Citizen Science", "In-Person"))


# ANOSIM, p = 0.0679
ano = anosim(COM.MATRIX, ENV.MATRIX$source, distance = "bray", permutations = 9999)
ano


# ADONIS
AD_COMP<-adonis(COM.MATRIX ~ source, data = ENV.MATRIX, permutations = 999, method="bray")
AD_COMP
tidy(AD_COMP)






### ---------- SECTION 6. BETA DISPERSION ---------- ###

# Convert to presence absence.
TOTAL_COMP_2 <- TOTAL_COMP
TOTAL_COMP_2[TOTAL_COMP_2>0]<- 1

# Create beta part object
COMP_CORE <- betapart.core(TOTAL_COMP_2[3:89])

# Returns three dissimilarity matrices containing 
# Pairwise between-site values of each beta-diversity component
COMP_DIST <- beta.pair(COMP_CORE, index.family = "sorensen")
str(COMP_DIST)

# Beta dispersion
BETA_SOR <- betadisper(COMP_DIST$beta.sor, TOTAL_COMP_2$county, type = c("median"))
BETA_SOR
anova(BETA_SOR)
plot(BETA_SOR)
boxplot(BETA_SOR, ylab = "Distance to median")
TukeyHSD(BETA_SOR, which = "group", conf.level = 0.95)

