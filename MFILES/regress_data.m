%% regress data does a GP regression on the data
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
        
testt=[-500:5:2020];

        trainspecs=1;
        trainlabels={'default'};

        regresssets=[1];
        regressparams=[1];
        clear regresslabels;

for iii=1:length(regresssets)
    ii=regresssets(iii);
    jj=regressparams(iii);

    clear wdataset;

    noiseMasks = ones(1,length(thetTGG{trainspecs(jj)}));

    defaultmask=1;
    wdataset=datasets{1};
     labls{iii}=['_default'];
    disp(labls{iii});

    trainsub = find(wdataset.limiting==0);
    wdataset.dY = sqrt(datasets{ii}.dY.^2); % + (thetTGG{jj}(end)*wdataset.compactcorr).^2);
    wdataset.Ycv = datasets{ii}.Ycv;% + diag(thetTGG{jj}(end)*wdataset.compactcorr).^2;
    subtimes=find(testt>=-5000+min(wdataset.time1));
    trainsub = find(datasets{1}.limiting==0);

 % find average height for each 1000 years -10050 = 12kaBP
    collinear =[];

    [f2s,sd2s,V2s,testlocs,logp(ii,jj),passderivs,invcv]=RegressHoloceneDataSets(datasets{1},testsitedef,modelspec(trainspecs(jj)),thetTGG{jj},trainsub,noiseMasks(1,:),testt,refyear,collinear);

    testX=[];
    testreg=testlocs.reg;
    testsites=testlocs.sites;
    testX(:,1:2)=testlocs.X(:,1:2);
    testX(:,3)=1950-testlocs.X(:,3);
    
% shift to put all plots at zero at 2000
    L = eye(length(testlocs.reg));
    sub_datum=find(testlocs.X(:,3)==2000);
    for i=1:length(testsitedef.names)
        sub=find(testlocs.reg==testsitedef.sites(i));
        subd=intersect(sub,sub_datum);
        L(sub,subd)=L(sub,subd)-1;
    end
    f2_shift=L*f2s;
    shifts=f2s(sub_datum);
    V2_shift=L*V2s*L';
    sd_shift=sqrt(diag(V2_shift));
    
    for jjj=1:length(shifts)
        sub=find(wdataset.datid==wdataset.siteid(jjj));
            wdataset.offset(sub)=shifts(jjj);
    end

 end

%     ii=regresssets(iii);
%     jj=regressparams(iii);
% 
%     noiseMasks = ones(4,length(thetTGG{trainspecs(jj)}));
% 
%     defaultmask=1;
%     
%     wdataset=datasets{1};
% 
%     disp(labls{iii});
% 
%     trainsub = find((wdataset.limiting==0)); % only index points
%     wdataset.dY = sqrt(wdataset.dY.^2 + (thetTGG{jj}(end)*wdataset.compactcorr).^2);
%     wdataset.Ycv = wdataset.Ycv + diag(thetTGG{jj}(end)*wdataset.compactcorr).^2;
%     subtimes=find(testt>=min(union(wdataset.time1,wdataset.time2)));
%     
%     [f2s{ii,jj},sd2s{ii,jj},V2s{ii,jj},testlocs{ii,jj},logp(ii,jj),passderivs,invcv]=RegressHoloceneDataSets(wdataset,testsitedef,modelspec(trainspecs(jj)),thetTGG{jj},trainsub,noiseMasks,testt(subtimes),refyear,[]);

    labl=labls{iii}; disp(labl);

    u=unique(wdataset.datid);
    clear fp sdp;
    clear testsitedefp;
    subps=[];
    for pp=1:length(u)
        subp=find(wdataset.datid==u(pp));
        subq=find(wdataset.siteid==u(pp));
        if length(subp)>0
            testtsp{pp}=wdataset.meantime(subp);
            testsitedefp.sites(pp,:)=[wdataset.siteid(subq) ...
                                wdataset.sitecoords(subq,:)];
            testsitedefp.names(pp)=wdataset.sitenames(subq);
            testsitedefp.names2(pp)=wdataset.sitenames(subq);
            testsitedefp.firstage(pp)=min(wdataset.meantime(subp));
            %            testsitedefp.GISfp(pp)=wdataset.siteGISfp(subq);
            %            testsitedefp.GIA(pp)=wdataset.siteGIA(subq);
        end
        subps=[subps ; subp];
    end

    [fp(subps),sdp(subps),~,testlocp]=RegressHoloceneDataSets(wdataset,testsitedefp,modelspec(trainspecs(jj)),thetTGG{jj},trainsub,noiseMasks(1,:),testtsp,refyear,collinear,passderivs,invcv);
    
    clear fp2 sdp2;
    [fp2(subps),sdp2(subps)]=RegressHoloceneDataSets(wdataset,testsitedefp,modelspec(trainspecs(jj)),thetTGG{jj},trainsub,noiseMasks(1,:),testtsp,refyear,collinear,passderivs,invcv);
    fpoffset = fp-fp2;

%% write predictions to a file
    fid=fopen(['ProxyData' labl '.tsv'],'w');
    fprintf(fid,['Site\tID\ttime1\ttime2\tlimiting\tY\' ...
                 'tdY\tYshifted\tcompactcorr\tistg\tlat\tlong\tmean time\tYposterior\tdYposterior\n']);
    for i=1:size(wdataset.datid,1)
        subq=find(wdataset.siteid==wdataset.datid(i));
        compactcorr=full(wdataset.compactcorr);
        subq=find(wdataset.siteid==wdataset.datid(i));
        fprintf(fid,[wdataset.sitenames{subq} '\t']);
        fprintf(fid,'%d\t',wdataset.datid(i));
        fprintf(fid,'%0.1f\t',wdataset.time1(i));
        fprintf(fid,'%0.1f\t',wdataset.time2(i));
        fprintf(fid,'%d\t',wdataset.limiting(i));
        fprintf(fid,'%0.2f\t',wdataset.Y(i));
        fprintf(fid,'%0.2f\t',wdataset.dY(i));
        fprintf(fid,'%0.2f\t',wdataset.Y(i)-fpoffset(i));       
        fprintf(fid,'%0.2f\t',compactcorr(i));
        fprintf(fid,'%d\t',wdataset.istg(i));
        fprintf(fid,'%0.2f\t',wdataset.lat(i));
        fprintf(fid,'%0.2f\t',wdataset.long(i));
        fprintf(fid,'%0.1f\t',wdataset.meantime(i));
        fprintf(fid,'%0.2f\t',fp(i));
        fprintf(fid,'%0.2f\n',sdp(i));
    end
    fclose(fid)    

xl=[-550 2010];
    wtestlocs=testlocs;
    subsites=[1];
    TG=PX; wdataset.offset=0*Y; Y0s=wdataset.Y0; testreg=testlocs.reg;

testt=[-500:5:2020];
tsitedef=testsitedef;
[f1s,sds,Vs,testlocs,logp,passderivs,invcv]=RegressHoloceneDataSets(datasets{1},tsitedef,modelspec(trainspecs(1)),thetTGG{1},trainsub,noiseMasks(1,:),testt,refyear,collinear);
time_to_run=toc;
