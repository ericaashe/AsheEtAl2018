%% Empirical Spatio-Temporal Gaussian Process analysis from Ashe et al., 2018

% change the '~/' to the directory where you have downloaded or cloned the
% main files.
% to run the optimization, uncomment line 58 and comment out line 62
% skip the optimization and use the parameters in the published model,
% comment out line 62 and uncomment line 58.

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
date_field='181105';
label1='EST_GP';
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
datPX=importdata(fullfile(CEIFILES,'US_Atlantic_Coast_for_ESTGP.csv'));

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

%  tic
%  [thetPX,trainsubset] = OptimizeHoloceneCovariance(datasets{1},modelspec(1),[2.4 2.0 3.4 3.0],[],[],[],[]);
%  opt_time=toc;
% thetTGG{1}=thetPX(1:end-1);

thetPX=[ 72749 32363 1424 4184 7.13 500.8 2502 0.1136 10.1631 0.1065];

thetTGG{1}=thetPX;

%%%%%%
%% do GP regression
%%%%%%

regress_data_ST;

for iii=1:size(noiseMasks,1)
    makeplot_slrate_ea(datasets{1},f2s{iii},sd2s{iii},V2s{iii},testlocs,lab{iii},2,2,22,[],0,0,[],[],[],[]);
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