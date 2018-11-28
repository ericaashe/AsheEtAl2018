%% regress_data_ST does a spatio-temporal GP regression (conditions the prior GP on the data)

ii=1;
trainsets=[1];
trainspecs=[1];
trainlabels={'default'};

%% do a regression

        clear Loc cnt Locid oldest youngest testsitedef;
        for ii=1:length(datasets{1}.sitenames)
            s0=0;
            Loc{ii}=datasets{1}.sitenames{ii};
            sub=strfind(Loc{ii},' ');
            Loc{ii}=Loc{ii}(setdiff(1:length(Loc{ii}),sub));
            sub=strfind(Loc{ii},'/');
            Loc{ii}=Loc{ii}(setdiff(1:length(Loc{ii}),sub));
            sub=find(datasets{1}.datid==datasets{1}.siteid(ii));
            cnt(ii)=length(sub);
            if length(sub)>0
                oldest(ii)=floor(min(datasets{1}.time1(sub)));
                youngest(ii)=ceil(max(datasets{1}.time2(sub)));
            else
                youngest(ii)=2000;
                oldest(ii)=2010;
            end
        end

        [uLoc,ui]=unique(Loc);

        clear testsitedef;
        testsitedef.sites=[];
        testsitedef.names={};
        testsitedef.names2={};
        testsitedef.firstage=[];

trainsub = find(datasets{1}.limiting==0);

        for ii=1:length(Loc)
                si=find(datasets{1}.datid==datasets{1}.siteid(ii)); 
                site_lat=datasets{1}.sitecoords(ii,1);%lat(si(sis));                
                site_long=datasets{1}.sitecoords(ii,2);%(si(sis));
                                si=si(1);
                testsitedef.sites(end+1,:)=[datasets{1}.datid(si) site_lat site_long];
                testsitedef.names2={testsitedef.names2{:}, datasets{1}.sitenames{ii}};
                testsitedef.names={testsitedef.names{:}, Loc{ii}};
                testsitedef.firstage = [testsitedef.firstage min(oldest(ii),-17000)];
        end
        
%% add Merritt Island as a prediction site:
    testsitedef.sites(end+1,:)=[1000 28.34 279.33];
    testsitedef.names2={testsitedef.names2{:}, 'Merritt-Island'};
    testsitedef.names={testsitedef.names{:}, 'Merritt-Island'};
    testsitedef.firstage = [testsitedef.firstage -17000];
        
testt=[-14000:50:2000];
if isfield(testsitedef,'GIA')
        if size(testsitedef.GIA,2)>1
            % matrix of projections for each site by times, not just a linear rate
        elseif size(testsitedef.GIA,2)==1
            % find projection at each time
            for ii=1:length(testsitedef.GIA)
                % matrix with a row for each site, column for each time
                testsitedef.GIAproj(ii,:)=-(ii).*(2010-testt);
            end
        end
end
        trainspecs=1;
        trainlabels={'default'};
        regresssets=[1];
        regressparams=[1];
        clear regresslabels;
        for i=1:length(regresssets)
            regresslabels{i} = [datasets{regresssets(i)}.label '_' trainlabels{regressparams(i)}];
        end
iii=1;
    ii=regresssets(iii);
    jj=regressparams(iii);

    clear wdataset;

    noiseMasks = ones(4,length(thetTGG{trainspecs(jj)}));
    noiseMasks(1,[11]) = 0;         % full model without white noise
    noiseMasks(2,[1 8 11]) = 0;     % All Regional only
    noiseMasks(3,[3 5 8 11]) = 0;   % Global only
    noiseMasks(4,[1 3 5 11]) = 0;   % Local

    defaultmask=1;
    wdataset=datasets{1};ii=1;
    labls{iii}=['_' regresslabels{iii}];
    disp(labls{iii});

    trainsub = find(wdataset.limiting==0);
    wdataset.dY = sqrt(datasets{ii}.dY.^2); 
    wdataset.Ycv = datasets{ii}.Ycv;
    subtimes=find(testt>=-5000+min(wdataset.time1));
    trainsub = find(datasets{1}.limiting==0);

    collinear =[];f2s={};sd2s={};V2s={};testlocs=[];
for iii=1:size(noiseMasks,1)
    [f2s{iii},sd2s{iii},V2s{iii},testlocs,logp(ii,jj),passderivs,invcv]=RegressHoloceneDataSets(datasets{1},testsitedef,modelspec(trainspecs),thetTGG{1},trainsub,noiseMasks(iii,:),testt,refyear,collinear);
end
lab=[{'Full'} {'Regional'} {'Global'} {'Local'}];
