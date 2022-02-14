# NEEDED PACKAGES.
library(rgbif)
library(plyr)


### ---------- SECTION 1. IMPORTING DATA ---------- ###

# Here, we import source data. 
# GBIF data contains observations from iNaturalist and Odonata Central. 
# Field observation data contains 'expert' in-field observations at each county site.
GBIF_OBS <- read.csv("citsci_obs_2019.csv", header=T)     # GBIF species occurrence.
FIELD_OBS <- read.csv("field_obs_2019.csv", header=T)     # Field observation data.

# Cleaning raw data.
# Subsetting both data sets.
GBIF_OBS <- subset(GBIF_OBS, select = c("species", "County", "month"))
FIELD_OBS <- subset(FIELD_OBS, select = c("county", "month", "scientific_name"))

# Adding source type to both data sets. 'CS' = Citizen science, 'IP' = In-person
GBIF_OBS$source="CS"
FIELD_OBS$source="IP"

# Reorganize columns.
GBIF_OBS <- GBIF_OBS[, c(2, 3, 4, 1)]
FIELD_OBS <- FIELD_OBS[, c(1, 2, 4, 3)]

# Standardize column names. 
names(GBIF_OBS)[names(GBIF_OBS) == "County"] <- "county"
names(FIELD_OBS)[names(FIELD_OBS) == "scientific_name"] <- "species"

# Standardize site names. 
### Citizen Science
GBIF_OBS$county[GBIF_OBS$county == 'Cuyahoga'] <- 'A'
GBIF_OBS$county[GBIF_OBS$county == 'Guernsey'] <- 'B'
GBIF_OBS$county[GBIF_OBS$county == 'Wayne'] <- 'C'
GBIF_OBS$county[GBIF_OBS$county == 'Knott'] <- 'D'
GBIF_OBS$county[GBIF_OBS$county == 'Wise'] <- 'E'
### In-person
FIELD_OBS$county[FIELD_OBS$county == 'Cuyahoga'] <- 'A'
FIELD_OBS$county[FIELD_OBS$county == 'Guernsey'] <- 'B'
FIELD_OBS$county[FIELD_OBS$county == 'Wayne'] <- 'C'
FIELD_OBS$county[FIELD_OBS$county == 'Knott'] <- 'D'
FIELD_OBS$county[FIELD_OBS$county == 'Wise'] <- 'E'





### ---------- SECTION 2. RICHNESS ---------- ###

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

# GLMs
### Without interaction.   AIC: 200.75
RICH_GLM1 <- glm.nb(formula = richness ~ source + county + as.factor(month), data = SUB_RICH)
summary(RICH_GLM1)
### With interaction.      AIC: 186.59
RICH_GLM2 <- glm.nb(formula = richness ~ source * county + as.factor(month), data = SUB_RICH)
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





### ---------- SECTION 3. ABUNDANCE ---------- ###

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

# GLM
### Without interaction.   AIC: 262.14
AB_GLM1 <- glm.nb(formula = abundance ~ source + county + as.factor(month), data = TOT_AB)
summary(AB_GLM1)
### With interaction.      AIC: 251.54
AB_GLM2 <- glm.nb(formula = abundance ~ source * county + as.factor(month), data = TOT_AB)
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





### ---------- SECTION 4. COMPOSITION ---------- ###

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
NMDS <- metaMDS(COM.MATRIX, distance="bray", k=2, autotransform=FALSE, trymax=100)
stressplot(NMDS)

# NMDS Visualization
plot(NMDS)
# add ellipsoids with ordiellipse
ordiellipse(NMDS, ENV.MATRIX$source, draw="polygon", col="#E69F00",kind="sd", conf=0.95, label=FALSE, show.groups = "CS")
ordiellipse(NMDS, ENV.MATRIX$source, draw="polygon", col="#009E73",kind="sd", conf=0.95, label=FALSE, show.groups = "IP") 
# display sampling source as shapes. CS = Circle. IP = Triangle. 
points(NMDS, display="sites", select=which(ENV.MATRIX$source=="CS"),pch=19, col="#E69F00")
points(NMDS, display="sites", select=which(ENV.MATRIX$source=="IP"), pch=17, col="#009E73")
# legend
legend(1.46,1.45, title=NULL, pch=c(19,17,15,25), col=c("#E69F00","#009E73"), cex=.7, legend=c("Citizen Science", "In-Person"))




# ADONIS
AD_COMP<-adonis(COM.MATRIX ~ source, data = ENV.MATRIX, permutations = 999, method="bray")
AD_COMP
