# AsheEtAl2018
Code for select models from 'Statistical modeling of rates and trends in Holocene relative sea level' in Quaternary Science Reviews.

There are three main files:
1.  runET_GP_CC.m   This code analyzes the continuous core data from either Northern North Carolina ('NNC_CC.csv') or New Jersey (NJ_CC.csv') with the Empirical Temporal Gaussian Process model described in sections 3.2.3, 4.3, and 5.1.
2.  runESTGP_AtlUS.m  This code analyzes proxy data from the Atlantic Coast of the U.S. in the file 'US_Atlantic_Coast_for_ESTGP.csv' with the Empirical Spatio-Temporal Gaussian Process model described in sections 3.3.1, 4.3, and 5.2.
3.  runESTGP_TG.m   This code analyses tide gauges along the Atlantic Coast of the U.S.  These tide guages must be downloaded from https://www.psmsl.org/data/obtaining/complete.php.  In this analysis we used annual RSL averages.  Unzip the rlr_annual file in the IFILES directory to run this code on it.


After running the chosen model, the results can be found within a folder where you are running (or have specified within) the code.

All other dependencies can be found in MFILES, and all data files needed to run this code are found in IFILES.

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
