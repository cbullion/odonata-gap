TTILE: Data Gap

ABSTRACT: 
Observations of Odonata (dragonflies and damselflies) across five focal counties for the summer months of 2019. Each observation possesses taxonomic information (genus, scientific name, suborder), spatial information (county and site), and date of the observation. The dataset consists of two files, representing their respective sampling source. 

-----

citsci_obs_2019.csv

Data collected from GBIF (https://doi.org/10.15468/dl.hmnc5b).

gbifID - Numerical ID of GBIF observation
datasetKey - Serial number of GBIF observation
occurrenceID - ID of occurrence 

kingdom - Kingdom-level classification of the insect
phylum - Phylum-level classification of the insect
class - Class-level classification of the insect
order - Order-level classification of the insect
family - Family-level classification of the insect
genus - Genus-level classification of the insect
species - Scientific name of the insect
infraspecificEpithet - Infraspecies information of the insect, if available
taxonrank - Classification level of the observation
scientificName - Full scientific name of the insect, with citation
verbatimScientificName - Full scientific name, without citation
verbatimScientificNameAuthorship - Empty and unused

countryCode - Location code for the country of observation (ex. US)
County - County of observation (ex. Cuyahoga)
stateProvince - State of observation (ex. Ohio)
occurrenceStatus - Presence or absence of insect
individualCount - Optional number of individuals observed
publishingOrgKey - Serial number of the publishing organization
decimalLatitude - Latitude of observation
decimalLongitude - Longitude of observation
coordinateUncertaintyInMeters - Geographical uncertainty (in meters) 
coordinatePrecision - N/A
elevation - N/A
elecationAccuracy - N/A
depth - N/A
depthAccuracy - N/A

eventDate - Full date of observation
day - DD of observation
month - MM of observation
year - YYYY of observation

taxonKey - GBIF ID of taxon
speciesKey - GBIF ID of species

basisOfRecord - Observation source type, human or indirect
institutionCode - Institution source of observation
collectionCode - Collection code
catalogNumber - Catalog number
recordNumber - N/A
identifiedBy - User associated with the observation identification
dateIdentified - Date submitted
license - License type of observation
rightsholder - Rights holder
recordedBy - User associated with the observation record
typeStatus - N/A
establishmentMeans - N/A
lastInterpreted - API processing timestamp
mediaType - Media Type
issue - Data flags

-----

field_obs_2019.csv

Data collected from in-field sampling.

observed_on - Date of observation
county - County of observation
month - Month of observation
year - Year of observation
site - Site label of observation
scientific_name - Scientific name of observation
taxon_genus_name - Genus name of insect
taxon_faily_name - Family name of insect
taxon_suborder_name - Suborder of insect