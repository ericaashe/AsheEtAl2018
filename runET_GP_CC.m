%% Empirical Temporal Gaussian Process Regression of Continuous Core Data from Ashe et al., 2018
% for independent runs of the continuous core data
% change label1 and datPX line to switch from New Jersey to Northern North
% Carolina data analysis

cd ('/Users/ericaashe/AsheEtAl2018');
addpath('/Users/ericaashe/AsheEtAl2018');
addpath('/Users/ericaashe/AsheEtAl2018/MFILES');
pd=pwd;
CEpd='/Users/ericaashe/AsheEtAl2018';

CEIFILES=[CEpd '/IFILES'];
IFILES=fullfile(pd,'IFILES');
date_field='181101';
label1='CC_NNC';
%label1='CC_NJ';
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

%datPX=importdata(fullfile(IFILES,'NJ_CC.csv'));%
datPX=importdata(fullfile(IFILES,'NNC_CC.csv'));%

%%%%%%
% prepare data
%%%%%%

prep_data;

%%%%%%
% define covariance function
%%%%%%

DefCovOneSite;

%%%%%%
% optimize covariance
%%%%%%

tic
[thetPX,trainsubset,logp] = OptimizeHoloceneCovariance(datasets{1},modelspec(1),[2.4 2.0 3.4 3.0]); %2.4 2.0 3.4 3.0
opt_time=toc;
thetTGG={};
thetTGG{1}=thetPX;
time_for_opt=toc;

%%%%%%
%% do GP regression
%%%%%%

regress_data;

%%%%%%
% plot results
%%%%%%

plot_figures;

%%%%%%
