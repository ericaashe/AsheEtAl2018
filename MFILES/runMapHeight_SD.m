
labl = 'Full';
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
firstyears=[-10000:4000:-2000];% 1000:200:1800 1900 ]; 
lastyears=[-6000:4000:2000];% 1100:200:2000 2000];

for i=1:length(firstyears)
    clf;
    nsitepts{i}=[];
    td = lastyears(i)-firstyears(i);
    addt=(2000-td)/2;
    fy=firstyears(i);%-addt;
    ly=lastyears(i);%+addt;
    subt = intersect(find(wdataset.time1<=ly),find(wdataset.time2>=fy));
    %find the number of sites with data influencing the model close to this time period
    sub=find(testsites(:,2)<=360);
    wmodelspec = modelspec(1);

    [fslopeFs,sdslopeF,fsF,sdsF,~,~,~,~,passderivs,invcv,fmeanFs,fsdF] = RegressRateField_ea(wdataset,wmodelspec,thetTGG{jj},noiseMasks(1,:),Flat,Flong,firstyears(i),lastyears(i),trainsub,[],[],[]);
    %[fslopeFs,sdslopeF,fsF,sdsF,~,~,~,~,passderivs,invcv,fmeanFs,fsdF] = RegressRateField_noGIA_ea(wdataset,wmodelspec,thetTGG{jj},noiseMasks(1,:),Flat,Flong,firstyears(i),lastyears(i),trainsub,[],[],[]);
    %[priorslope,sdpriorslope,priorfsF,priorsdsF] = RegressRateField_noGIA_ea(nulldataset,wmodelspec,thetTGG{jj},noiseMasks(1,:),Flat,Flong,firstyears(i),lastyears(i),trainsub,[],[],[]);
    [priorslope,sdpriorslope,priorfsF,priorsdsF] = RegressRateField_ea(nulldataset,wmodelspec,thetTGG{jj},noiseMasks(1,:),Flat,Flong,firstyears(i),lastyears(i),trainsub,[],[],[]);

   
     fmeanF = (fmeanFs)./1000;
     fmeanSL = 0;
     fslopeF  = fslopeFs;
     fsdF = fsdF./1000;
     sdslope = sdslopeF./1000;
    % Start mapping %%%

    clf;
    ax = worldmap([minlat maxlat],[minlong maxlong]);
    setm(ax,'meridianlabel','off','parallellabel','off','flinewidth',3);
    land = shaperead('landareas', 'UseGeoCoords', true);
    hold on;

    Flat1=min(Flat):max(Flat);
    Flong1=min(Flong):max(Flong);

    [FLONG,FLAT]=meshgrid(Flong,Flat);
    [FLONG1,FLAT1]=meshgrid(Flong1,Flat1);
    
    %mapped = griddata(FLONG(:),FLAT(:),fmeanF(:)+GIAproj'*.500*((1950-firstyears(i)+(1950-lastyears(i)))),Flong,Flat(:),'linear'); %+fslopeGSL;
    mapped = griddata(FLONG(:),FLAT(:),fmeanF(:)+fmeanSL,Flong,Flat(:),'linear'); %+fslopeGSL;
    mapRate = griddata(FLONG(:),FLAT(:),fslopeF(:),Flong,Flat(:),'linear'); %+fslopeGSL;
    sdmapped = griddata(FLONG(:),FLAT(:),fsdF(:),Flong,Flat(:),'linear');
    sdRatemapped = griddata(FLONG(:),FLAT(:),sdslope(:),Flong,Flat(:),'linear');
    %mapped = griddata(FLONG(:),FLAT(:),fsdF(:),Flong,Flat(:),'linear');

    u=sdmapped/mean(priorsdsF);
    %    threshold=quantile(u(:),.2);
    %    threshold=sqrt(.67);
    threshold=quantile(u(:),.85);

    subbad=find(u>threshold);
    mapped(subbad)=NaN;
    mapRate(subbad)=NaN;

    hs1=scatterm(FLAT(:),FLONG(:),250,mapped(:),'filled','marker','s');
%%% threshold of these points is over the limit
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
if i==1
    caxis([-31 -24]);
elseif i==2
    caxis([-13 -6]);
else
    caxis([-4 -0.5]);
end

title({['SL Height' num2str(firstyears(i)) ' to ' num2str(lastyears(i)) ' CE']});
pdfwrite(['SL_Height_no_land_' num2str(firstyears(i)) '_to_' num2str(lastyears(i))]);

scatterm(lat_site,long_site,3*ndat,'MarkerEdgeColor','k','MarkerFaceColor','w');
pdfwrite(['SL_Height_no_land_data_' num2str(firstyears(i)) '_to_' num2str(lastyears(i))]);
box on;

      box on;
       geoshow(ax, land, 'FaceColor', [0.85 0.85 0.85]);
       scatterm(lat_site,long_site,3*ndat,'MarkerEdgeColor','k','MarkerFaceColor','w');

       datehs = (1950-firstyears(i)+1950-lastyears(i))/2;
       datefirst = 1950-firstyears(i);
       datelast = 1950-lastyears(i);

           %num2str(datefirst) ' to ' num2str(datelast) ' BP']});
       pdfwrite(['Height_Sea_Level' labl '_' num2str(firstyears(i)) '_to_' num2str(lastyears(i)) '_no_data']);

    clf;
    ax = worldmap([minlat maxlat],[minlong maxlong]);
    setm(ax,'meridianlabel','off','parallellabel','off','flinewidth',3);
    land = shaperead('landareas', 'UseGeoCoords', true);
    hold on;

    hs1=scatterm(FLAT(:),FLONG(:),250,mapRate(:),'filled','marker','s');

     hcb=colorbar;
     if i==1
         caxis([5.4 7.4]);
     elseif i==2
        caxis([2 3.8]);
     else
        caxis([0 2]);
     end
     
       title({['Rate of RSL Change (m/ky)' num2str(firstyears(i)) ' to ' num2str(lastyears(i)) ' CE']});
           %num2str(datefirst) ' to ' num2str(datelast) ' BP']});
       pdfwrite(['Rate_no_land_' labl '_' num2str(firstyears(i)) '_to_' num2str(lastyears(i)) '_scaled']);   
     
       geoshow(ax, land, 'FaceColor', [0.85 0.85 0.85]);

       title({['Rate of RSL Change (m/ky)' num2str(firstyears(i)) ' to ' num2str(lastyears(i)) ' CE']});
           %num2str(datefirst) ' to ' num2str(datelast) ' BP']});
       pdfwrite(['Rate_RSL_Change' labl '_' num2str(firstyears(i)) '_to_' num2str(lastyears(i)) '_scaled']);   
    
  %% plot standar deviations
    clf;
    ax = worldmap([minlat maxlat],[minlong maxlong]);
    setm(ax,'meridianlabel','off','parallellabel','off','flinewidth',3);
    land = shaperead('landareas', 'UseGeoCoords', true);
    hold on;

    hs1=scatterm(FLAT(:),FLONG(:),250,sdmapped(:),'filled','marker','s');
    hcb=colorbar;
    caxis([0 4]);
    title({['SD of Sea Level Change ' num2str(firstyears(i)) ' to ' num2str(lastyears(i)) ' BP']});
    pdfwrite(['SL_SD_no_land_' num2str(firstyears(i)) '_to_' num2str(lastyears(i))]);

    hold on;   %     
    geoshow(ax, land, 'FaceColor', [0.85 0.85 0.85]);
    hold on;
    axis tight;
    hcb=colorbar;

    box on;
       box on;

    title({['SD of Sea Level Change ' num2str(firstyears(i)) ' to ' num2str(lastyears(i)) ' BP']});
    pdfwrite(['SL_SD_' num2str(firstyears(i)) '_to_' num2str(lastyears(i))]);

end
