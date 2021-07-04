# Script to follow "Biodiversity Facets tutorial"
# Authors: Matthias Grenié & Marten Winter

# Setup ------------------------------------------------------------------------
# Install packages
install.packages(
  c("ade4",
    "ape",
    "FD",
    "fundiversity",
    "ggplot2",
    "ggspatial",
    "performance",
    "picante",
    "rnaturalearth",
    "sf")
)


# Loading the data -------------------------------------------------------------
# Load data
plot_data         = read.csv("data/doi_10.5061_dryad.f77p7__v1/PlotData.csv",
                             na.strings = c("NA", "na"),
                             stringsAsFactors = TRUE)
plot_species_data = read.csv("data/doi_10.5061_dryad.f77p7__v1/PlotSpeciesData.csv")
species_traits    = read.csv("data/doi_10.5061_dryad.f77p7__v1/SpeciesTraitData.csv",
                             na.strings = c("NA", "na"), stringsAsFactors = TRUE)

# Describe the data
str(plot_data)
summary(plot_data)

str(plot_species_data[, 1:5])
summary(head(plot_species_data)[,1:5])
dim(plot_species_data)

# Transform one column for further analyzed
species_traits$seed = ordered(species_traits$seed)
str(species_traits)
summary(species_traits)

#' **Questions for you:**
#' 
#' * **Q1**: How many plots were sampled?
#' * **Q2**: How many species are there in the dataset?
#' * **Q3**: How many traits are available?
#' * **Q4**: How many of them are continuous? How many of them are discrete?
#' * **Q5**: What is the most numerous family among all observed species?
#' * **Q6**: What is the most numerous genus?
#'

# Description of the environment variables
boxplot(forestloss17 ~ block, data = plot_data,
        xlab = "Block of plot", ylab = "Forest loss (%)",
        main = "Forest loss in funciton of block of data")


# Functional diversity ---------------------------------------------------------
## Data wrangle
# Make site-species data.frame
sp_com           = plot_species_data[, -1]
rownames(sp_com) = plot_species_data$X
sp_com = as.matrix(sp_com)

# Make synthesized trait data.frame
traits = species_traits[, -c(1:5)]
rownames(traits) = species_traits$species.code

## Get CWM
# Get only continuous CWM
quanti_cwm = FD::functcomp(traits[, c("height", "sla", "wood.dens")],
                           sp_com, CWM.type = "dom")
quanti_cwm$plot.code = rownames(quanti_cwm)

## CWM-environment value
# Merge environmental data with CWM
cwm_env = merge(
  quanti_cwm,
  plot_data[, c("plot.code", "block", "forestloss17", "roaddensprim")],
  by = "plot.code"
)

## Plot CWM-environment relationships
par(mfrow = c(2, 2))
plot(cwm_env$forestloss17, cwm_env$height,
     xlab = "Forest loss (%)", ylab = "Biomass-weighted height",
     main = "CWM Height vs. forest loss")
plot(cwm_env$forestloss17, cwm_env$sla,
     xlab = "Forest loss (%)", ylab = "Biomass-weighted SLA",
     main = "CWM SLA vs. forest loss")
plot(cwm_env$forestloss17, cwm_env$wood.dens,
     xlab = "Forest loss (%)", ylab = "Biomass-weighted wood density",
     main = "CWM Wood density vs. forest loss")
plot(cwm_env$roaddensprim, cwm_env$height,
     xlab = "Road density (km.km^-2)", ylab = "Biomass-weighted height",
     main = "CWM Height vs. road density")

#' **Questions for you:**
#' 
#' * **Q7**: How would you describe the relationship between the different CWMs
#'           and forest loss?
#' * **Q8**: Can you test the correlation using the function `cor.test()`
#'           and does it support your previous statements?
#' * **Q9**: How would you describe the understorey vegetation changes
#'           with increasing forest loss?

# Correlation between CWMs and environments
cor.test(cwm_env$forestloss17, cwm_env$height)
cor.test(cwm_env$forestloss17, cwm_env$sla)
cor.test(cwm_env$forestloss17, cwm_env$wood.dens)
cor.test(cwm_env$roaddensprim, cwm_env$height)

# Non-quantitative CWM
non_quanti_cwm = FD::functcomp(traits[, -c(5:7)],
                               sp_com, CWM.type = "all")
non_quanti_cwm$plot.code = rownames(non_quanti_cwm)

non_quanti_cwm = merge(
  non_quanti_cwm,
  plot_data[, c("plot.code", "block", "forestloss17", "roaddensprim")],
  by = "plot.code"
)

# Visualization Categorical CWM-environment relationship
par(mfrow = c(1, 1))
plot(non_quanti_cwm$forestloss17, non_quanti_cwm$woody_no,
     xlab = "Forest loss (%)", ylab = "Sum of biomass of non-woody species",
     main = "Biomass of non-woody species vs. forest loss")

#' **Questions for you**
#' 
#' * **Q10**: How does this observation compare to above description
#'            of the change of understorey vegetation along
#'            the forest loss gradient?

## Building functional space
# Gower's Dissimilarity
gower_dissim = cluster::daisy(traits)

# Perform the PCOA
trait_pcoa = ade4::dudi.pco(ade4::quasieuclid(gower_dissim), nf = 3,
                            scannf = FALSE)
trait_pcoa

# Visualize PCoA
ade4::scatter(trait_pcoa, clab.row = 0)

# Woodiness PCoA
ade4::s.class(trait_pcoa$li[,1:2], fac = traits$pgf)

#' **Questions for you**
#' 
#' * **Q11**: Using the metadata available in the `README.txt` file,
#'            what is the meaning `pgf` column means?
#' * **Q12**: How do you interpret the PCoA results given your answer
#'            to the previous question?

# Compute Functional Richness (FRic)
site_fric = fundiversity::fd_fric(trait_pcoa$li, sp_com, stand = FALSE)

# Compute Functional Evenness (FEve) and Rao's Quadratic Entropy (Q)
site_raoq = fundiversity::fd_raoq(trait_pcoa$li, sp_com)
site_feve = fundiversity::fd_feve(trait_pcoa$li, sp_com)

site_fd = merge(
  merge(site_fric, site_raoq, by = "site"),
  site_feve,
  by = "site"
)
site_fd$plot.code = site_fd$site
site_fd = site_fd[, -1]

# FD-environment relationships
site_env_fd = merge(site_fd,
                    plot_data[, c("plot.code", "forestloss17", "roaddensprim")],
                    by = "plot.code")

par(mfrow = c(2, 2))
plot(site_env_fd$forestloss17, site_env_fd$FRic,
     xlab = "Forest loss (%)", ylab = "Functional Richness (FRic)",
     main = "Functional Richness vs. forest loss")
plot(site_env_fd$forestloss17, site_env_fd$Q,
     xlab = "Forest loss (%)", ylab = "Rao's Quadratic Entropy",
     main = "Q vs. forest loss")
plot(site_env_fd$forestloss17, site_env_fd$FEve,
     xlab = "Forest loss (%)", ylab = "Functional Evenness (FEve)",
     main = "FEve vs. forest loss")
plot(site_env_fd$roaddensprim, site_env_fd$FRic,
     xlab = "Primary Road Density (km.km^-2)", ylab = "Functional Richness (FRic)",
     main = "FRic vs. road density")

#' **Questions for you**
#' 
#' * **Q13**: How would you describe the relationships between
#'            functional diversity and forest loss and road density?
#' * **Q14**: Using the plot generated by the code beneath how could you
#'            describe the relationships between the three different
#'            functional diversity indices we computed?
#' 

# Pairwise relaitonships between functional diversity indices
panel.cor = function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y, use = "complete.obs"))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

pairs(~FEve + Q + FRic, data = site_env_fd, lower.panel = panel.smooth,
      upper.panel = panel.cor, gap = 0, row1attop = FALSE)

# Correlation with taxa richness
site_rich_fd = merge(
  site_fd,
  plot_data[, c("plot.code", "ntaxa")],
  by = "plot.code"
)

pairs(ntaxa ~ FRic + FEve + Q, data = site_rich_fd, upper.panel = panel.cor)

# Relationship between functional diversity and biomass
site_biomass = rowSums(sp_com)
site_biomass = stack(site_biomass)

site_biomass$plot.code = site_biomass$ind
site_biomass$tot_biomass   = site_biomass$values

site_biomass = site_biomass[, c("plot.code", "tot_biomass")]

site_rich_fd = merge(
  site_rich_fd,
  site_biomass,
  by = "plot.code"
)

pairs(tot_biomass ~ FRic + FEve + Q, data = site_rich_fd,
      upper.panel = panel.cor)

#' **Questions for you**
#' 
#' * **Q15**: How does the relationship between indices with species richness
#'            compare with the one observed with total biomass values?
#'            (You can use the function `cor.test()` if you want
#'            to test the association)

### Null modeling

# Compute null traits
# Set random seed so that everybody gets the same null traits
set.seed(20210705)

# Number of null simulations
# CAUTION: increasing this number may increase future computation time by a lot
n_null = 99

# Repeat the operation as many times as set aboev
null_traits = lapply(seq.int(n_null), function(x) {
  null_trait = trait_pcoa$li

  # Shuffle species names
  null_species = sample(rownames(trait_pcoa$li), nrow(trait_pcoa$li))

  # Replace species name in table
  rownames(null_trait) = null_species

  # Do not forget to return the modified table!
  return(null_trait)
})

str(null_traits, max.l = 0)
head(null_traits[[1]])

# Compute null fd indices
# Beware this make take a long time
null_fd = lapply(seq(length(null_traits)), function(y) {

  x = null_traits[[y]]

  null_fric = fundiversity::fd_fric(x, sp_com, stand = FALSE)
  null_raoq = fundiversity::fd_raoq(x, sp_com)
  null_feve = fundiversity::fd_feve(x, sp_com)

  # Combine all null functional diversity values
  null_all = merge(
    merge(null_fric, null_raoq, by = "site"), null_feve, by = "site"
  )

  # Null Index to separate between all null simulations
  null_all$null_id = y

  return(null_all)
})

null_fd_all = do.call(rbind.data.frame, null_fd)
head(null_fd_all)

# Load pre-computed 999 iterations of null FD indices
null_fd_999 = readRDS("data/null_fd_999.Rds")

head(null_fd_999)

## Compare null and observed FD
# The observed value of FRic for the site
subset(site_fd, plot.code == "a100f177r")$FRic

# The null distribution of FRic for the same site
summary(subset(null_fd_999, site == "a100f177r")$FRic)

# Histogram null FRic
par(mfrow = c(1, 1))
# Visualize histogram of null values
hist(subset(null_fd_999, site == "a100f177r")$FRic,
     breaks = 20,
     xlab = "null Functional Richness",
     ylab = "Frequency",
     main = "FRic comparison for site 'a100f177r'")
abline(v = subset(site_fd, plot.code == "a100f177r")$FRic,
       col = "darkred", lwd = 2)

#' 
#' **Question for you**
#' 
#' * **Q16**: How would describe verbally the position of the observed value
#'            of FRic for site "a100f177r" compared to the null distribution?

# Build the Empirical Cumulative Distribution Function
# Build the ECDF
one_null_fric_ecdf = ecdf(subset(null_fd_999, site == "a100f177r")$FRic)

# Then actually use it
obs_fric = subset(site_fd, plot.code == "a100f177r")$FRic

one_null_fric_ecdf(obs_fric)

#' **Question for you**
#' 
#' * **Q17**: What the quantile of the observed FRic value in the end?

# Compute FD Standardized Effect Sizes
# Compute average and standard deviation of null distribution
mean_null_fd = aggregate(
  cbind(mean_FRic = FRic, mean_Q = Q, mean_FEve = FEve) ~ site,
  data = null_fd_999, FUN = mean, na.rm = TRUE
)
sd_null_fd   = aggregate(
  cbind(sd_FRic = FRic, sd_Q = Q, sd_FEve = FEve) ~ site, data = null_fd_999,
  FUN = sd, na.rm = TRUE
)

# Merge null mean & sd with observed values
obs_null_fd = merge(
  site_fd,
  merge(mean_null_fd, sd_null_fd, by = "site"),
  by.x = "plot.code", by.y = "site"
)

# Compute SES
obs_null_fd$ses_FRic = (obs_null_fd$mean_FRic - obs_null_fd$FRic)/obs_null_fd$sd_FRic
obs_null_fd$ses_Q = (obs_null_fd$mean_Q - obs_null_fd$Q)/obs_null_fd$sd_Q
obs_null_fd$ses_FEve = (obs_null_fd$mean_FEve - obs_null_fd$FEve)/obs_null_fd$sd_FEve

# Cleaner table
ses_fd = obs_null_fd[, c("plot.code", "FRic", "Q", "FEve", "ses_FRic", "ses_Q",
                         "ses_FEve")]

#' **Question for you**
#' 
#' * **Q18**: Using the `subset()` function with the greater (or equal) than `>=`
#'            and the lower (or equal) than `<=`, can you determine
#'            how many sites show a significant deviation from the null observation?
#'            (absolute SES >= 2)
#' * **Q19**: Using similar code as used for observed values,
#'            what are the relationships between SES values and forest loss?


## Mapping functional diversity
# Get plot spatial coordinates
head(plot_data[, c(1, 4, 5)])

plot_sf = sf::st_as_sf(
  plot_data[, c(1:7)],
  coords = c("north", "east"),
  crs = sf::st_crs("+proj=utm +zone=50 +datum=WGS84 +units=m +no_defs")
)

# See plots location in the world
library("ggplot2")

ggplot() +
  geom_sf(data = rnaturalearth::ne_countries(returnclass = "sf")) +
  geom_sf(data = plot_sf, aes(color = forestloss17)) +
  scale_color_viridis_c() +
  coord_sf(crs = sf::st_crs("+proj=eck4")) +  # Set projection
  labs(title = "Map of the concerned plots at world scale") +
  theme_bw()

# Focused map on Malaysia
ggplot() +
  geom_sf(data = rnaturalearth::ne_countries(continent = "Asia",
                                             returnclass = "sf")) +
  geom_sf(data = plot_sf, aes(color = forestloss17)) +
  scale_color_viridis_c() +
  coord_sf(crs = sf::st_crs(3376), xlim = c(-1072025.83, 1053446.00),
           ylim = c(85496.43, 767752.41)) +
  labs(title = "Map of plots focused on Malaysia")
ggspatial::annotation_scale() +
  theme_bw()

# Zoomed Map
ggplot() +
  geom_sf(data = rnaturalearth::ne_countries(country = "Malaysia",
                                             returnclass = "sf")) +
  geom_sf(data = plot_sf, aes(color = forestloss17)) +
  scale_color_viridis_c() +
  coord_sf(crs = sf::st_crs(3376), xlim = c(800000, 890000),
           ylim = c(500000, 550000)) +
  labs(title = "Map of plots zoomed-in on Sabah region") +
  ggspatial::annotation_scale() +
  ggspatial::annotation_north_arrow(location = "br") +
  theme_bw()

# Map with context
ggplot() +
  ggspatial::annotation_map_tile(zoomin = -1) +
  geom_sf(data = plot_sf, aes(color = forestloss17)) +
  scale_color_viridis_c() +
  coord_sf(crs = sf::st_crs(3376), xlim = c(800000, 890000),
           ylim = c(500000, 550000)) +
  labs(title = "Map of plots zoomed-in on Sabah region") +
  ggspatial::annotation_scale() +
  ggspatial::annotation_north_arrow(location = "br") +
  theme_bw()

# Better context map
ggplot() +
  ggspatial::annotation_map_tile(zoomin = -1) +
  geom_sf(data = subset(plot_sf, block != "og"),
          aes(color = forestloss17)) +
  scale_color_viridis_c() +
  coord_sf(crs = sf::st_crs(3376), xlim = c(875000, 890000),
           ylim = c(518500, 531000)) +
  labs(title = "Map of all plots but block 'og'") +
  ggspatial::annotation_scale() +
  ggspatial::annotation_north_arrow(location = "br") +
  theme_bw()

# Map of functional diversity
ggplot() +
  geom_sf(
    data = merge(subset(plot_sf, block != "og"), ses_fd, by = "plot.code"),
    aes(color = ses_Q)
  ) +
  scale_color_distiller(type = "div", palette = "RdYlBu",
                        name =  "SES of Rao's Quadratic Entropy") +
  coord_sf(crs = sf::st_crs(3376), xlim = c(875000, 890000),
           ylim = c(518500, 531000)) +
  labs(title = "Map of all plots but block 'og'") +
  ggspatial::annotation_scale() +
  ggspatial::annotation_north_arrow(location = "br") +
  theme_gray()

# Phylogenetic diversity -------------------------------------------------------
## Getting the phylogenetic tree
# Load phylogenetic tree
phylo_tree = ape::read.tree("data/doi_10.5061_dryad.f77p7__v1/phylo_tree.nwk")

phylo_tree
str(phylo_tree)

#' **Questions from you**
#' 
#' * **Q20**: How many taxa are in the phylogenetic tree?
#' * **Q21**: How does this number compare to the number of taxa found in the
#'            dataset?

# Get species in phylogenetic tree
phylo_tree$tip.label

#' * **Q22**: What do you notice with the species names?
#'            Especially compared to the ones available in `species_traits`.

## Get corresponding code from phylogenetic tree
# Create an indexed list of names
phylo_names = species_traits[, c("species.code", "species")]
phylo_names$code_id = seq(nrow(phylo_names))

# Get the first species code based on species epithet
code_id_to_use = aggregate(code_id ~ species, phylo_names,
                           FUN = function(x) head(x, 1))

# Get back the data.frame of species names with the actual species.code
code_species = merge(
  code_id_to_use, phylo_names[, c("code_id", "species.code")], by = "code_id"
)

# Tidying code for edge cases
code_species$species = gsub(" ", "", code_species$species)
code_species$species = paste0(
  tolower(substr(code_species$species, 1, 1)),
  substr(code_species$species, 2, nchar(code_species$species))
)

code_species = code_species[, c("species.code", "species")]

dim(code_species)

# Intersection between phylogenetic tree and species in rest of data
length(intersect(phylo_tree$tip.label, code_species$species))

### Visualizing the phylogenetic tree
# Plot phylogenetic tree
ape::plot.phylo(phylo_tree)

# Better phylogenetic tree
ape::plot.phylo(phylo_tree, type = "fan", show.node.label = TRUE,
                show.tip.label = FALSE, cex = 0.6)


### Computing phylogenetic diversity indices
# Subset phylogenetic matrix
# Initial site-species matrix
head(sp_com[, 1:5])
dim(sp_com)

# Subset of site-species matrix compatible with phylogenetic tree
sub_phylo_com = sp_com[, as.character(code_species$species.code)]
dim(sub_phylo_com)

# Cophenetic distance
# Compute cophenetic distances from the phylogenetic tree
cophen_dist = ape::cophenetic.phylo(phylo_tree)

str(cophen_dist)

# We need to change the names to species codes
corres_codes = data.frame(
  species = rownames(cophen_dist)
)
corres_codes = merge(corres_codes, code_species, by = "species")
rownames(cophen_dist) = corres_codes$species.code
colnames(cophen_dist) = corres_codes$species.code

## Compute Mean Pairwise Distance
# Observed Mean Pairwise Distance
# Unweighted
mpd_val_uw = picante::mpd(sub_phylo_com, cophen_dist, abundance.weighted = FALSE)
# Weighted
mpd_val_w = picante::mpd(sub_phylo_com, cophen_dist, abundance.weighted = TRUE)

# Make a nice data.frame with observed MPD values
obs_mpd = data.frame(
  plot.code = rownames(sub_phylo_com),
  mpd_unweighted = mpd_val_uw,
  mpd_weighted = mpd_val_w
)

# Add forest loss proportion and richness for each site
obs_mpd = merge(obs_mpd, plot_data[, c("plot.code", "forestloss17", "ntaxa")])

#' **Questions for you**
#' 
#' * **Q23**: What are the relationships between weighted and
#'            unweighted version of the MPD?
#' * **Q24**: What are the relationships between MPD and taxa richness?
#'            And with forest loss? Plot these relationships to visualize them
#'            and use the `cor.test()` function to validate your observations.

### Null modeling 
# Compute SES of MPD
# Set random seed for repeatability of analysis
set.seed(20210705)

# Compute null permutation of MPD
ses_mpd = picante::ses.mpd(
  sub_phylo_com, cophen_dist, null.model = "taxa.labels",
  abundance.weighted = TRUE, runs = 99
)
head(ses_mpd)

# Load the pre-computed null MPDs
ses_mpd_999 = readRDS("data/null_mpd_999.Rds")

#' **Questions for you**
#' 
#' * **Q25**: Explain what does the column `mpd.obs.z` means?
#'            How does this compare with the SES values we computed
#'            for functional diversity indices?
#' * **Q26**: How does the standardized value relates with taxa richness?
#' * **Q27**: What are the relationships between MPD values considering
#'            null models and forest loss?
#'            Visualize the relationships with the `plot()` function,
#'            validate your observations with the `cor.test()` function.

## Mapping phylogenetic diversity
# Map of SES MPD
ses_mpd_999$plot.code = rownames(ses_mpd_999)

ggplot() +
  geom_sf(
    data = merge(subset(plot_sf, block != "og"), ses_mpd_999, by = "plot.code"),
    aes(color = mpd.obs.z)
  ) +
  scale_color_distiller(type = "div", palette = "RdYlBu",
                        name =  "SES of MPD") +
  coord_sf(crs = sf::st_crs(3376), xlim = c(875000, 890000),
           ylim = c(518500, 531000)) +
  labs(title = "Map of all plots but block 'og'") +
  ggspatial::annotation_scale() +
  ggspatial::annotation_north_arrow(location = "br") +
  theme_gray()


# Comparing facets -------------------------------------------------------------

# Combine taxonomic, functional, and phylogenetic diversity
all_diversity = merge(
  plot_data[, c("plot.code", "ntaxa")],
  merge(
    ses_fd, ses_mpd_999[, -1], by = "plot.code"
  )
)

# Comparison of observed values
pairs(all_diversity[, c("ntaxa", "FRic", "Q", "FEve", "mpd.obs")],
      upper.panel = panel.cor)

# Comparison of SESs
pairs(all_diversity[, c("ntaxa", "ses_FRic", "ses_Q", "ses_FEve", "mpd.obs.z")],
      upper.panel = panel.cor)

#' **Question for you**
#' 
#' * **Q28**: How are related are observed values of functional diversity and
#'            phylogenetic diversity? What about the SESs?

# Modelling the effect of logging ----------------------------------------------

# Merge diversity and environmental variables
plot_div_env = merge(
  all_diversity,
  plot_data[, c(1, 6:21)],
  by = "plot.code"
)

dim(plot_div_env)
head(plot_div_env)


## Single-predictor models
# Model taxa richness in function of forest loss
mod_taxa_loss = lm(ntaxa ~ forestloss17, data = plot_div_env)

mod_taxa_loss
summary(mod_taxa_loss)

# Plot model taxa-loss
par(mfrow = c(1, 1))
plot(mod_taxa_loss$model$forestloss17, mod_taxa_loss$model$ntaxa,
     xlab = "Forest Loss (%)", ylab = "Taxa Richness")
abline(coef = coef(mod_taxa_loss), col = "darkred", lwd = 1)

#' **Questions for you**
#' 
#' * **Q29**: How would you qualify the effect of forest loss on the taxa richness?
#' * **Q30**: With the same formula build similar models with the other
#'            predictors `roaddensprim` and `roaddistprim`.
#'            How do they compare with forest loss?

# Model taxa richness in function of other disturbances
mod_taxa_dens = lm(ntaxa ~ roaddensprim, data = plot_div_env)
mod_taxa_dist = lm(ntaxa ~ roaddistprim, data = plot_div_env)

# Model of FD and PD in function of forest loss
mod_fd_loss = lm(ses_Q ~ forestloss17, data = plot_div_env)
mod_pd_loss = lm(mpd.obs.z ~ forestloss17, data = plot_div_env)

## Multi-predictors models
# Model taxa richness in function of all predictors
mod_taxa_all = lm(
  ntaxa ~ elevation + forestloss17 + forestloss562 + roaddenssec +
    roaddistprim + soilPC1 + soilPC2,
  data = plot_div_env
)
# Model FD in function of all predictors
mod_fd_all = lm(
  ses_Q ~ elevation + forestloss17 + forestloss562 + roaddenssec +
    roaddistprim + soilPC1 + soilPC2,
  data = plot_div_env
)
# Model PD in function of all predictors
mod_pd_all = lm(
  mpd.obs.z ~ elevation + forestloss17 + forestloss562 + roaddenssec +
    roaddistprim + soilPC1 + soilPC2,
  data = plot_div_env
)

#' **Question for you**
#' 
#' * **Q31**: What can you say about the effect of the disturbances
#'            on the different diversity metrics?
#'            What are the explanatory power of our models?

# Model diagnostics
par(mfrow = c(2, 2))
plot(mod_taxa_all)
# Or even better
performance::check_model(mod_taxa_all)
