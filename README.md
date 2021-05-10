# AsheEtAl2019
This repository contains the MATLAB code base for select models from [Ashe et al. 2019]:
'Statistical modeling of rates and trends in Holocene relative sea level'
 
 Ashe et al. (2019) is published in Quaternary Science Reviews.

[![DOI](https://zenodo.org/badge/156264585.svg)](https://zenodo.org/badge/latestdoi/156264585)

**Project abstract:**

> Characterizing the spatio-temporal variability of relative sea level (RSL) and estimating local, regional,
and global RSL trends requires statistical analysis of RSL data. Formal statistical treatments, needed to
account for the spatially and temporally sparse distribution of data and for geochronological and elevational
uncertainties, have advanced considerably over the last decade. Time-series models have
adopted more flexible and physically-informed specifications with more rigorous quantification of uncertainties.
Spatio-temporal models have evolved from simple regional averaging to frameworks that
more richly represent the correlation structure of RSL across space and time. More complex statistical
approaches enable rigorous quantification of spatial and temporal variability, the combination of
geographically disparate data, and the separation of the RSL field into various components associated
with different driving processes.We review the range of statistical modeling and analysis choices used in
the literature, reformulating them for ease of comparison in a common hierarchical statistical framework.
The hierarchical framework separates each model into different levels, clearly partitioning measurement
and inferential uncertainty from process variability. Placing models in a hierarchical
framework enables us to highlight both the similarities and differences among modeling and analysis
choices. We illustrate the implications of some modeling and analysis choices currently used in the
literature by comparing the results of their application to common datasets within a hierarchical
framework. In light of the complex patterns of spatial and temporal variability exhibited by RSL, we
recommend non-parametric approaches for modeling temporal and spatio-temporal RSL.

If you have any questions, comments, or feedback on this work or code, please [contact Erica](mailto:ericaashe@gmail.com) 

### Dependencies
All dependencies can be found in MFILES, and all data files needed to run this code are found in IFILES.

## File Descriptions

There are three main files:
1.  runET_GP_CC.m   This code analyzes the continuous core data from either Northern North Carolina ('NNC_CC.csv') or New Jersey (NJ_CC.csv') with the Empirical Temporal Gaussian Process model described in sections 3.2.3, 4.3, and 5.1.
2.  runESTGP_AtlUS.m  This code analyzes proxy data from the Atlantic Coast of the U.S. in the file 'US_Atlantic_Coast_for_ESTGP.csv' with the Empirical Spatio-Temporal Gaussian Process model described in sections 3.3.1, 4.3, and 5.2.
3.  runESTGP_TG.m   This code analyses tide gauges along the Atlantic Coast of the U.S.  These tide guages must be downloaded from https://www.psmsl.org/data/obtaining/complete.php.  In this analysis we used annual RSL averages.  Unzip the rlr_annual file in the IFILES directory to run this code on it.

After running the chosen model, the results can be found within a folder where you are running (or have specified within) the code.


## Authors

### Contributors
* **Erica Ashe, PhD** - *Co-author, Bayesian Statistics* - [GitHub](https://github.com/ericaashe)

### Co-authors
* **Niamh Cahill**
* **Nicole S. Khan**
* **Carling Hay**
* **Andrew Kemp**
* **Simon E. Engelhart**
* **Benjamin P. Horton**
* **Andrew C, Parnell**
* **Robert E. Kopp**

> Copyright (C) 2018 by Erica L. Ashe
> This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
> This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
> A copy of the GNU General Public License comes with this program.  If not, see <http://www.gnu.org/licenses/>.

## Acknowledgments

This work was supported by National Science Foundation grants
OCE-1458904 (ELA, REK, and BPH), OCE-1458903 (SEE), OCE-
1458921 (ACK), OCE-1702587 (ELA and REK), Singapore Ministry of
Education Academic Research Fund Tier 2 MOE218-T2-1-030, the
National Research Foundation of Singapore, and the Singapore
Ministry of Education under the Research Centres of Excellence
initiative (NSK and BPH), Science Foundation Ireland Career
Development Award grant number 17/CDA/4695 (ACP), and is a
contribution to IGCP Project 639, INQUA Project CMP1601P “HOLSEA”
and PALSEA3.
