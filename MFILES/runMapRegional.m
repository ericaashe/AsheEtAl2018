
labl = 'Regional';
    ndatapoints =[];
    for j=1:length(wdataset.siteid)
      ndatapoints(j)=length(intersect(find(wdataset.limiting==0),find(wdataset.datid==wdataset.siteid(j))));
    end        
    lim_data_sub = find(ndatapoints>0);

nulldataset=SubsetDataStructure(wdataset,1,1);
nulldataset.meantime=0; nulldataset.dt=0; nulldataset.dY=200e3; nulldataset.limiting=0;

trainsub = find(wdataset.limiting==0);
cc=find(wdataset.istg==2);
tg=find(wdataset.istg==1);
di=find(wdataset.istg==0);
firstyears=[-10000:4000:-2000];
lastyears=[-6000:4000:2000];

for i=1:length(firstyears)
    clf;
    nsitepts{i}=[];
    td = lastyears(i)-firstyears(i);
    addt=(2000-td)/2;
    fy=firstyears(i);
    ly=lastyears(i);
    subt = intersect(find(wdataset.time1<=ly),find(wdataset.time2>=fy));
    
    %find the number of sites with data influencing the model close to this time period
    sub=find(testsites(:,2)<=360);
    wmodelspec = modelspec(1);

    [fslopeFs,sdslopeF,fsF,sdsF,~,~,~,~,passderivs,invcv,fmeanFs,fsdF] = RegressRateField_ea(wdataset,wmodelspec,thetTGG{jj},noiseMasks(3,:),Flat,Flong,firstyears(i),lastyears(i),trainsub,[],[],[]);
    [priorslope,sdpriorslope,priorfsF,priorsdsF] = RegressRateField_ea(nulldataset,wmodelspec,thetTGG{jj},noiseMasks(1,:),Flat,Flong,firstyears(i),lastyears(i),trainsub,[],[],[]);

     fmeanF = (fmeanFs)./1000;
     fmeanSL = 0;
     fslopeF  = fslopeFs;
     fsdF = fsdF./1000;
     sdslope = sdslopeF./1000;
     
    %% Start mapping %%%

    clf;
    ax = worldmap([minlat maxlat],[minlong maxlong]);
    setm(ax,'meridianlabel','off','parallellabel','off','flinewidth',3);
    land = shaperead('landareas', 'UseGeoCoords', true);
    hold on;

    Flat1=min(Flat):max(Flat);
    Flong1=min(Flong):max(Flong);

    [FLONG,FLAT]=meshgrid(Flong,Flat);
    [FLONG1,FLAT1]=meshgrid(Flong1,Flat1);
    
    mapped_Reg = griddata(FLONG(:),FLAT(:),fmeanF(:)+fmeanSL,Flong,Flat(:),'linear'); %+fslopeGSL;
    mapRate = griddata(FLONG(:),FLAT(:),fslopeF(:),Flong,Flat(:),'linear'); %+fslopeGSL;
    sdmapped = griddata(FLONG(:),FLAT(:),fsdF(:),Flong,Flat(:),'linear');
    sdRatemapped = griddata(FLONG(:),FLAT(:),sdslope(:),Flong,Flat(:),'linear');

    u=sdmapped/mean(priorsdsF);
    threshold=quantile(u(:),.85);

    subbad=find(u>threshold);
    mapped(subbad)=NaN;
    mapRate(subbad)=NaN;

    hs1=scatterm(FLAT(:),FLONG(:),250,mapped_Reg(:),'filled','marker','s');
    hold on;
    sublong=find((mod(FLONG1(:),5)==0).*(mod(FLAT1(:),5)==0));
    hold on;    
    box on;
    
    subt = intersect(find(wdataset.time1<=ly),find(wdataset.time2>=fy));
    
    % map the regions where there is data and how much
    for j=1:length(wdataset.siteid)
      nsitepts{i}(j)=length(intersect(intersect(intersect(find(wdataset.limiting==0),find(wdataset.time1<=ly)),find(wdataset.time2>=fy)),find(wdataset.datid==wdataset.siteid(j))));
    end        

    subsite = find(nsitepts{i}>0);
    ndat = nsitepts{i}(subsite)';
    lat_site = wdataset.sitecoords(subsite,1);
    long_site = wdataset.sitecoords(subsite,2);

    uns = unique(ndat);
    
 
    box on;   
    axis tight;
    hcb=colorbar;
    caxis([-4 8]);

    title({['Regiona_Signal' num2str(firstyears(i)) ' to ' num2str(lastyears(i)) ' CE']});
    pdfwrite(['Regional_no_land_' num2str(firstyears(i)) '_to_' num2str(lastyears(i))]);

    scatterm(lat_site,long_site,3*ndat,'MarkerEdgeColor','k','MarkerFaceColor','w');
    pdfwrite(['Regional_no_land_data_' num2str(firstyears(i)) '_to_' num2str(lastyears(i))]);
    box on;

      box on;
      geoshow(ax, land, 'FaceColor', [0.85 0.85 0.85]);
      scatterm(lat_site,long_site,3*ndat,'MarkerEdgeColor','k','MarkerFaceColor','w');

      datehs = (1950-firstyears(i)+1950-lastyears(i))/2;
      datefirst = 1950-firstyears(i);
      datelast = 1950-lastyears(i);

      pdfwrite(['Regional_RSL' labl '_' num2str(firstyears(i)) '_to_' num2str(lastyears(i)) '_no_data']);
end
