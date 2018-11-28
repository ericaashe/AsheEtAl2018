%% Analysis of Tide Gauges along the Atlantic Coast of the United States from Ashe et al., 2018

% If using tide gauges downloaded from https://psmsl.org/data/, place 'rlr_annual' directory with all 
% subdirectories in IFILES directory and uncomment line 11 of load_TGs.m

% To skip the optimization, comment out lines 49 - 51, and uncomment line 53

% change the '~/' to the directory where you have downloaded or cloned the main files:
% cd ('~/AsheEtAl2018/IFILES');
% load('TG_data.mat');
% addpath('~/AsheEtAl2018');
% addpath('~/AsheEtAl2018/MFILES');
% pd=pwd;
% CEpd='~/AsheEtAl2018';

%% example:
cd ('/Users/ericaashe/AsheEtAl2018');
addpath('/Users/ericaashe/AsheEtAl2018');
addpath('/Users/ericaashe/AsheEtAl2018/MFILES');
pd=pwd;
CEpd='/Users/ericaashe/AsheEtAl2018';

date_field = '1801114';
label = 'TG';    

%%%%%%
%% load tide gauge data
%%%%%%

load_TGs;

cd(pd);
WORKDIR=[label date_field];
if ~exist(WORKDIR,'dir')
    mkdir(WORKDIR);
end
cd(WORKDIR);

%%%%%%
%% Define Covariance Functions
%%%%%%

TGDefineCovFuncs;

%%%%%%
%% Optimize Covariance Function
%%%%%%

  tic
  [thetTG,trainsubset,logp] = OptimizeHoloceneCovariance(datasets{1},modelspec(1),[2.4 2.0 3.4 3.0],[],[],0,0,0); 
  opt_time=toc

%thetTG = [185.9522   94.7445   21.8909   55.3918    2.0854    3   0.1553    3.4539    0];

thetTGG{1}=thetTG(1:end-1);

datasets{2}=datasets{1}; testt=[1900:1:2010];

%%%%%%
%% now do a regression
%%%%%%

TGRegression;

%  makeTGplots;
for iii=1:size(noiseMasks,1)
    makeplot_slrate(datasets{1},f2s{iii},sd2s{iii},V2s{iii},testlocs,lab{iii},2,1,4,[],0,[],[],[],[]);
end

makeplot_tg(datasets{1},f2s{1},sd2s{1},V2s{1},testlocs,lab{1},2,1,4,[],0,0,[],[],[],[]);

prepTGrunMaps;
