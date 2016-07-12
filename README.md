Zika_nmicrobiol_2016
====================

This repository contains code used in the following paper.

Perkins TA, Siraj AS, Ruktanonchai CW, Kraemer MUGK, Tatem AJ. (2016) **Model-based projections of Zika virus infections in childbearing women in the Americas**. *Nature Microbiology* XX:YY-ZZ. doi:[http://www.nature.com/nmicrobiol/](http://www.nature.com/nmicrobiol/)

All code contained within this repository is released under the [CRAPL v0.1 License](http://matt.might.net/articles/crapl/). Because of the large sizes of many files used in the paper, we have included only a subset of data.


====================

### code/models folder

The scripts in this folder produce Figures 2-4 and S1. More generally, they are used to fit relationships between seroprevalence from the 13 sites described in Table S1 and the covariates (i.e., temperature, economic index, *Aedes aegypti* occurrence probability) under both mechanistic and statistical model formulations. The outputs of these models are projected attack rates following the first wave of the epidemic. These scripts were run in the following order and made use of both a personal laptop (Mac OSX) and the University of Notre Dame's Center for Research Computing cluster [http://crc.nd.edu](http://crc.nd.edu).

* `runjob.pbs`
* `script.R`
* `0_relationship_R0_AR.R`
* `1_params_random_draws.R`
* `2_fit_attackrate_random_draws.R`
* `3_fig_infection_distributions.R`
* `4_fit_attackrate_random_draws_stat.R`
* `5_fig_infection_distributions_stat.R`
* `6_fig_attackrate_relationships.R`


### code/maps folder

The scripts in this folder produce the maps used in Figures 1 & S2-S10 and the numbers used in Figures 1-3. These scripts were run in the following order and made use of both a personal laptop (Mac OSX) and the University of Notre Dame's Center for Research Computing cluster [http://crc.nd.edu](http://crc.nd.edu).

* `0_numfunctions.R`
* `1_raster_aligning.R`
* `run2job.pbs` calls `2_fit_attackrates_seroprev.R`
* `run2stat.pbs` calls `2_fit_attackrates_seroprev_stat.R`
* `2_pop_projection.R`
* `3_fill_blank_econs.R`
* `run4job.pbs` calls `4_output_grids_crc.R`
* `run4stat.pbs` calls `4_output_grids_stat_crc.R`
* `run5job.pbs` calls `5_min_max_mean_1st_round.R`
* `5_run_AR_R0.R`
* `run6job.pbs` calls `6_min_max_mean_2nd_round.R`
* `6_fig_maps.R`
* `run7job.pbs` calls `7_min_max_mean_surface.R`
* `run9job.pbs` calls `9_median_1st_round.R`
* `9_median_surface.R`
* `10_country_summary.py`


### data folder

Data included here pertain to the 13 sites listed in Table S1 from which we obtained seroprevalence estimates, 100 replicates of *Aedes aegypti* occurrence probabilities from those sites from [Kraemer et al. (2015)](https://elifesciences.org/content/4/e08347), and a generalized additive model object that describes the relationship between temperature and adult female *Aedes aegypti* mortality from [Brady et al. (2013)](https://parasitesandvectors.biomedcentral.com/articles/10.1186/1756-3305-6-351).


### generated folder

Files included here contain 1,000 replicates of parameterizations of the mechanistic and statistical models, respectively.


### outputs folder

The file included here contains the country-level sums of total infections and infections among childbearing women used in Figures 2 and 3.


### maps folder

This folder contains files that can be used to access raster data shown in the maps in Figures S2-S10. The numbers in each 5x5 km grid cell are only pertinent to that grid cell, and totals across multiple grid cells cannot be interpreted as sums across larger areas. The reason is that, for example, the map of median attack rates contains each grid cell's median across 1,000 replicates rather than a reflection of what might somehow be considered a median spatial layer of attack rates across all grid cells. By contrast, distributions of country- and continent-level totals in Figures 1b, 2, and 3 reflect each such quantity calculated in each of the 1,000 replicates and then examined as a distribution.

