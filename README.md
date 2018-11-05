# AsheEtAl2018
Code for select models from 'Statistical modeling of rates and trends in Holocene relative sea level' in Quaternary Science Reviews.

There are three main files:
1.  runET_GP_CC.m   This code analyzes the continuous core data from either Northern North Carolina ('NNC_CC.csv') or New Jersey (NJ_CC.csv') with the Empirical Temporal Gaussian Process model described in sections 3.2.3, 4.3, and 5.1.
2.  runESTGP_AtlUS.m  This code analyzes proxy data from the Atlantic Coast of the U.S. in the file 'US_Atlantic_Coast_for_ESTGP.csv' with the Empirical Spatio-Temporal Gaussian Process model described in sectiona 3.3.1, 4.3, and 5.2.
3.  runESTGP_TG.m   This code analyses tide gauges along the Atlantic Coast of the U.S., downloaded from https://www.psmsl.org/data/.

After running the chosen model, the results can be found within a folder where you are running (or have specified within) the code.

All other dependencies can be found in MFILES, and all data files needed to run this code are found in IFILES.
