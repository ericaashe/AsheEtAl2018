%% Empirical Spatio-Temporal Gaussian Process analysis from Ashe et al., 2018

% change the '~/' to the directory where you have downloaded or cloned the main files.
% To run the optimization, uncomment lines 56 - 59 and comment out lines 62 - 63
% To skip the optimization and use the parameters in the published model,
% comment out line 56 - 59 and uncomment line 62 - 63

cd ('~/AsheEtAl2018');
addpath('~/AsheEtAl2018');
addpath('~/AsheEtAl2018/MFILES');
pd=pwd;
CEpd='~/AsheEtAl2018';

%% for example:
%cd ('/Users/ericaashe/AsheEtAl2018');
%addpath('/Users/ericaashe/AsheEtAl2018');
%addpath('/Users/ericaashe/AsheEtAl2018/MFILES');
%pd=pwd;
%CEpd='/Users/ericaashe/AsheEtAl2018';

CEIFILES=[CEpd '/IFILES'];
IFILES=fullfile(CEpd,'IFILES');
date_field='181114';
label1='_EST_GP';
label=label1;
WORKDIR=[date_field label1];
if ~exist(WORKDIR,'dir')
    mkdir(WORKDIR);
end
cd(WORKDIR);
idHolo = 3e4;

defval('firsttime',-12000);

% import proxy data

idHolo = 3e4;
datPX=importdata(fullfile(CEIFILES,'USAtl_ESTGP.csv'));

%%%%%%
% prepare data
%%%%%%

prep_data_ST;

%%%%%%
% define covariance function
%%%%%%

DefCovST;

%%%%%%
% optimize covariance
%%%%%%

%   tic
%   [thetPX,trainsubset,logp] = OptimizeHoloceneCovariance(datasets{1},modelspec(1),[2.4 2.0 3.4 3.0],[],[],0,0,0); 
%   opt_time=toc;
%   thetTGG{1}=thetPX(1:end-1);

 thetPX=[260224.3 59388.15 0.01134 12 5602.32 4247.65 12 500.13 993.57 3 100.01 0];
 thetTGG{1}=thetPX(1:end-1);

%%%%%%
%% do GP regression (condition prior GP with optimized hyperparameters on the data)
%%%%%%

regress_data_ST;

%% set limits for maps similar to published maps

ylims=[ -40e3 1000;      % Full
        -15e3 20e3;     % Non-linear regional
        -15e3 20e3;     % Regional
        -60e3 10e3;     % Global
        -2e3 1.5e3;     % Local
        -40e3 1e3]; 

for iii=1:size(noiseMasks,1)
    if iii==1
        plot_dat=1;
    else
        plot_dat=0;
    end
    makeplot_slrate(datasets{1},f2s{iii},sd2s{iii},V2s{iii},testlocs,lab{iii},2,2,22,meanSL,plot_dat,[],[],[],[],ylims(iii,:));
end

testX=[];
    testreg=testlocs.reg;
    testsites=testlocs.sites;
    testX(:,1:2)=testlocs.X(:,1:2);
    testX(:,3)=1950-testlocs.X(:,3);
    firstyr = [-10000:4000:-2000];
    lastyr = [-6000:4000:2000];
    
%%%%%%
%% Plot Maps
%%%%%%

runPlotMaps;

save;
