*************************************************************************************************************************************
ReadMe file accompanying data file for: 
"Logging increases the functional and phylogenetic dispersion of understorey plant communities in tropical lowland rainforest."

Authors: Timm F. Dobert, Bruce L. Webber, John B. Sugau, Katharine J.M. Dickinson & Raphael K. Didham

Published in: Journal of Ecology
DOI: 
*************************************************************************************************************************************
File name: PlotData.csv
Description: This csv file contains the full database of characteristics for each of 180 vegetation plots including logging metrics, environmental variables as well as taxonomic, functional and phylogenetic diversity indices.

Contact authors: Timm Dobert (timm.dobert@gmail.com); Raphael Didham (raphael.didham@csiro.au)

Column headings:
plot.code: Vegetation plot unique identifier
block: One of eight sampling blocks (A-F, LF and OG) across the study landscape 
fragment: One of three fragment designations (1ha, 10ha or 100ha) within a block
north: UTM northing coordinate
east: UTM easting coordinate 
elevation: The elevation of a plot (m) 
forestloss17: The proportion of forest canopy loss at the local scale, i.e. within a 17m radius (%)
forestloss562: The proportion of forest canopy loss at the landscape scale, i.e. within a 562m radius (%) 
roaddensprim: Regional-scale density of continuously-maintained primary logging roads (km.km-2)
roaddenssec: Regional-scale density of occasionally-used secondary logging roads (km.km-2)
roaddistprim: Distance to nearest primary road (m)
roaddistsec: Distance to nearest secondary road (m)
river: Distance to nearest river (> 3m width; m) 
pH: Soil pH
humus: Soil humus depth (cm) 
soilN: Total soil nitrogen content (mg.cm-3)
soilP: Plant available soil phosphorus content (µg.cm-3)
bulkdens: Soil bulk density (g.cm-3)
soilmoist: Soil moisture factor
soilPC1: Principal component axis 1 scores from a PCA of the six soil-biogeochemical variables 
soilPC2: Principal component axis 2 scores from a PCA of the six soil-biogeochemical variables 
ntaxa: Species richness 
FDisLOG: Biomass-weighted functional dispersion (standardized effect sizes)
nriLOG: Biomass-weighted net relatedness index (phylogenetic dispersion based on the branch length between each pair of taxa in a plot)
FRicLOG: Biomass-weighted functional richness (standardized effect sizes)
FDisLOGr: Biomass-weighted functional dispersion of response traits (standardized effect sizes)
FDisLOGe: Biomass-weighted functional dispersion of effect traits (standardized effect sizes)
woody: The proportion of woody species in a plot 
treesapling: The proportion of tree saplings out of all species per plot
nriLOGWoody: Biomass-weighted net relatedness index of woody species only
nriLOGSapling: Biomass-weighted net relatedness index of tree saplings only
nriLOGnonWoody: Biomass-weighted net relatedness index of non-woody species only
************************************************************************************************************************************
File name: SpeciesTraitData.csv 
Description: This csv file contains the complete list of species (691) sampled across 180 vegetation plots (2 x 2m) located at the Stability of Altered Forest Ecosystems (SAFE) project in Sabah, Malaysia, including their allocation to a range of plant functional traits as well as a distinction between native or exotic species origin.

Contact authors: Timm Dobert (timm.dobert@gmail.com); Raphael Didham (raphael.didham@csiro.au)

Column headings:
species.code: Unique code for each plant taxa
family: Plant family name
genus: Plant genus name 
species: Plant species name or unique identifier where species indeterminate
species.name: The genus and species name
tree: Distinction between tree (yes), no tree (no) or indeterminate (na)
woody: Distinction between woody (yes), non-woody (no) or indeterminate (na)
origin: Distinction between native (n) and exotic (e) plant species
pgf: Plant growth form: A = fern, B = graminoid, C = forb, D = herbaceous climber, E = herbaceous shrub, F = tree sapling, G = woody climber, H = woody shrub, na = indeterminate 
height: Maximum plant height (m)
sla: Specific leaf area (m2.kg-1) 
wood.dens: Wood density (g.cm-3)
dispersal: Predominant dispersal mode: J = animal, K = ant, L = ballistic, M = bat, N = bird, O = primate, P = water, Q = wind, na = indeterminate
fruit: The type of fruit: R = achene, S = berry, T = berry-like, U = capsule, V = caryopsis, W = drupe, X = follicle, Y = legume, Z = nut, a = samara, b = schizocarp, na = indeterminate
seed: The number of seeds per fruit: 1 = 1, 2 = <4, 3 = <10, 4 = >10, na = indeterminate
reproduction: The reproduction strategy: s = seed, v = vegetative, sv = seed or vegetative, na = indeterminate 
lifehistory: The lifehistory strategy: a = annual, abp = annual or biennial or perennial, ap = annual or perennial, p = perennial, na = indeterminate
pollination: The pollination syndrome: c = bat, d = bee, e = beetle, f = bird, g = butterfly, h = entomophilous.broad, i = entomophilous.narrow, j = fly, k = moth, l = passive, m = self, n = thrip, o = wasp, p = wind, na = indeterminate 
***********************************************************************************************************************************
File name: PlotSpeciesData.csv
Description: This csv file contains a matrix of biomass values for 691 plant taxa sampled across 180 vegetation plots (2 x 2m) located at the Stability of Altered Forest Ecosystems (SAFE) project in Sabah, Malaysia.

Contact authors: Timm Dobert (timm.dobert@gmail.com); Raphael Didham (raphael.didham@csiro.au)

Column and row headings:
Columns: The unique identifiers of 691 plant taxa
Rows: The 180 vegetation plots across 3 fragments nested within 8 blocks 
Matrix values: The dry weight biomass (g.m-2) of each species in each plot
***********************************************************************************************************************************
