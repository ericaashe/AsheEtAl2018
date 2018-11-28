%% Empirical Temporal Gaussian Process Regression of Continuous Core Data from Ashe et al., 2018
% single site for comparison to the spatio-temporal model

cd ('~/AsheEtAl2018');
addpath('~/AsheEtAl2018');
addpath('~/AsheEtAl2018/MFILES');
pd=pwd;
CEpd='~/AsheEtAl2018';

%% example:
% cd ('/Users/ericaashe/AsheEtAl2018');
% addpath('/Users/ericaashe/AsheEtAl2018');
% addpath('/Users/ericaashe/AsheEtAl2018/MFILES');
% pd=pwd;
% CEpd='/Users/ericaashe/AsheEtAl2018';
%%

CEIFILES=[CEpd '/IFILES'];
IFILES=fullfile(pd,'IFILES');
date_field='181114';

label1='_NC_proxy';
label=label1;
WORKDIR=[date_field label1];
if ~exist(WORKDIR,'dir')
    mkdir(WORKDIR);
end
cd(WORKDIR);
idHolo = 3e4;

%%%%%%
% import proxy data, read files
%%%%%%

datPX=importdata(fullfile(IFILES,'NC_proxies.csv'));%

%%%%%%
% prepare data
%%%%%%

prep_data;

%%%%%%
% define covariance function
%%%%%%

DefCovOneSite;

%%%%%%
tic
[thetPX,trainsubset,logp] = OptimizeHoloceneCovariance(datasets{1},modelspec(1),[2.4 2.0 3.4 3.0]);
opt_time=toc

%thetPX= [15623.2 4827.4 19.75 222.52 1.14137 0.03517];

thetTGG={};
thetTGG{1}=thetPX;

%%%%%%
%% do GP regression
%%%%%%

regress_proxy_data;

%%%%%%
% plot results
%%%%%%

plot_ET_proxy_figures;

%%%%%%
save;

