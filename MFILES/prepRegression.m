%% set up for regression of TGs   
clear Loc cnt Locid oldest youngest testsitedef;
        distfrom = dDist(datasets{2}.sitecoords,datasets{2}.sitecoords);
        for ii=1:length(datasets{2}.sitenames)
            s0=0;
            Loc{ii}=datasets{2}.sitenames{ii};
            sub=strfind(Loc{ii},' ');
            Loc{ii}=Loc{ii}(setdiff(1:length(Loc{ii}),sub));
            sub=strfind(Loc{ii},'/');
            Loc{ii}=Loc{ii}(setdiff(1:length(Loc{ii}),sub));
            sub=find(datasets{2}.datid==datasets{2}.siteid(ii));
            cnt(ii)=length(sub);
            if length(sub)>0
                oldest(ii)=floor(min(datasets{2}.time1(sub)));
                youngest(ii)=ceil(max(datasets{2}.time2(sub)));
            else
                youngest(ii)=2000;
                oldest(ii)=2010;
            end
        end

        clear oldestnear;
        for ii=1:length(datasets{2}.sitenames)
            sub=find(distfrom(ii,:)<2);
            oldestnear(ii) = min(oldest(sub));
        end

[uLoc,ui]=unique(Loc);
trainsub = find(datasets{2}.limiting==0);
        clear testsitedef;
        testsitedef.sites=[];
        testsitedef.names={};
        testsitedef.names2={};
        testsitedef.firstage=[];

        for ii=1:length(Loc)
                si=find(datasets{2}.datid==datasets{2}.siteid(ii)); 
                site_lat=datasets{2}.sitecoords(ii,1);%lat(si(sis));                
                site_long=datasets{2}.sitecoords(ii,2);%(si(sis));
%                 site_lat=mean(datasets{2}.lat(sis));                
%                 site_long=mean(datasets{2}.long(sis));
                                si=si(1);
                testsitedef.sites(end+1,:)=[datasets{2}.datid(si) site_lat site_long];
                    %datasets{2}.lat(si) datasets{2}.long(si)];
                testsitedef.names2={testsitedef.names2{:}, datasets{2}.sitenames{ii}};
                testsitedef.names={testsitedef.names{:}, Loc{ii}};
                testsitedef.firstage = [testsitedef.firstage min(oldest(ii),-17000)];
        end

        trainspecs=1;
        trainlabels={'default'};

        regresssets=[1];
        regressparams=[1];
        clear regresslabels;
        for i=1:length(regresssets)
            regresslabels{i} = [datasets{regresssets(i)}.label '_' trainlabels{regressparams(i)}];
        end
testt=[1900:1:2010];
%for iii=1:length(regresssets)
iii=1;
    ii=regresssets(iii);
    jj=regressparams(iii);

    clear wdataset;

noiseMasks = ones(4,length(thetTGG{trainspecs(jj)}));
noiseMasks(1,[7 9]) = 0;    % full model without white noise
noiseMasks(2,[4 7 9]) = 0;  % low-frequency term only
noiseMasks(3,[1 7 9]) = 0;  % medium-frequency term only
noiseMasks(4,[1 4 7 9]) = 0;  % local offset only

    defaultmask=1;
    wdataset=datasets{2};ii=2;
    labls{iii}=['_' regresslabels{iii}];
    disp(labls{iii});

    trainsub = find(wdataset.limiting==0);
    % find(wdataset.limiting==0 & wdataset.indicator==1);% only index points & normally distr
    wdataset.dY = sqrt(datasets{ii}.dY.^2); % + (thetTGG{jj}(end)*wdataset.compactcorr).^2);
    wdataset.Ycv = datasets{ii}.Ycv;% + diag(thetTGG{jj}(end)*wdataset.compactcorr).^2;
    subtimes=find(testt>=-5000+min(wdataset.time1));
    trainsub = find(datasets{2}.limiting==0);

subsites=[1 5 11 15];
tsitedef.sites=testsitedef.sites(subsites,:);
tsitedef.names=testsitedef.names(:,subsites);
tsitedef.names2=testsitedef.names2(:,subsites);
tsitedef.firstage=testsitedef.firstage(:,subsites);

    collinear =[];f2s={};sd2s={};V2s={};testlocs=[];
for iii=1:size(noiseMasks,1)
    [f2s{iii},sd2s{iii},V2s{iii},testlocs,logp(ii,jj),passderivs,invcv]=RegressHoloceneDataSets(datasets{1},tsitedef,modelspec(trainspecs),thetTGG{1},trainsub,noiseMasks(iii,:),testt,refyear,collinear);
end
lab=[{'Full'} {'Low-Freq'} {'Med-Freq'} {'Local-Offset'}];

 
%     ii=regresssets(iii);
%     jj=regressparams(iii);
%     
%     wdataset=datasets{2};
% 
%     labls{iii}=['_' regresslabels{iii}];
%     disp(labls{iii});
% 
%     trainsub = find((wdataset.limiting==0)); % only index points
%     [fp(subps),sdp(subps),~,testlocp]=RegressHoloceneDataSets(wdataset,testsitedefp,modelspec(trainspecs(jj)),thetTGG{jj},trainsub,noiseMasks(1,:),testtsp,refyear,collinear,passderivs,invcv);
%     
%     clear fp2 sdp2;
