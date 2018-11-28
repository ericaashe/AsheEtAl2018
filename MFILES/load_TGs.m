%% load_TGs loads the tide gauges to be used within the model

cd ('IFILES');
load('TG_data.mat');

targcoord= [39.355 -74.4183;    % Altantic City, NJ
34.2267 -77.9533;               % Wilmington, NC
41.36   -72.09];                % New London, CT

%% uncomment the following line of code if you are using tide gauge files from https://psmsl.org/data/
% [TG,TG0,thetL,TGmodellocal] = GPSmoothNearbyTideGauges(targcoord,396,[],2,1,[500 80 78],1.0,'rlr_annual','none',5,[]);

clear datasets;

datasets{1}=TG;
datasets{1}.label='TG';

for ii=1:length(datasets)
    datasets{ii}.Y0=datasets{ii}.Y;
    t1=datasets{ii}.time1; t2=datasets{ii}.time2;
    datasets{ii}.long = mod(datasets{ii}.long,360); sub=find(datasets{ii}.long>180); datasets{ii}.long(sub)=datasets{ii}.long(sub)-360;
    datasets{ii}.meantime=mean([t1 t2],2);
    datasets{ii}.dt = abs(t1-t2)/4;
    datasets{ii}.compactcorr=sparse(datasets{ii}.compactcorr);
    datasets{ii}.GIAproj = 0*datasets{ii}.Y;
    datasets{ii}.Y=datasets{ii}.Y0-datasets{ii}.GIAproj;
end