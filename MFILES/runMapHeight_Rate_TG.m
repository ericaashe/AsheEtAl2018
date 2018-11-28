labl = 'AtlCoast';
    ndatapoints =[];
    for j=1:length(wdataset.siteid)
      ndatapoints(j)=length(intersect(find(wdataset.limiting==0),find(wdataset.datid==wdataset.siteid(j))));
    end
    lim_data_sub = find(ndatapoints>0);

nulldataset=SubsetDataStructure(wdataset,1,1);
nulldataset.meantime=0; nulldataset.dt=0; nulldataset.dY=200e3; nulldataset.limiting=0;

trainsub = find(wdataset.limiting==0);
firstyears=1901;
lastyears=2000;
    
GIA=[];
if exist('GIALat')
    GIA.Lat=GIALat;
    GIA.Lon=mod(GIALon,360);
    GIA.GIA=GIA6;
end

for i=1:length(firstyears)
    clf;
    nsitepts{i}=[];
    td = lastyears(i)-firstyears(i);
    addt=(2000-td)/2;
    fy=firstyears(i)-addt;
    ly=lastyears(i)+addt;
    subt = intersect(find(wdataset.time1<=ly),find(wdataset.time2>=fy));
    %find the number of sites with data influencing the model close to this time period
    sub=find(testsitedef.sites(:,2)<=360);
       datehs = (1950-firstyears(i)+1950-lastyears(i))/2;
       datefirst = firstyears(i);
       datelast = lastyears(i);
    wmodelspec = modelspec;

    [fslopeF,sdslopeF,fsF,sdsF,~,~,~,~,passderivs,invcv,fmeanFs,fsdF,fGIAs] = RegressRateField_ea(wdataset,wmodelspec,thetTGG{jj},noiseMasks(1,:),Flat,Flong,firstyears(i),lastyears(i),trainsub,GIA);
    [priorslope,sdpriorslope,priorfsF,priorsdsF] = RegressRateField_ea(nulldataset,wmodelspec,thetTGG{jj},noiseMasks(1,:),mean(Flat),mean(Flong),firstyears(i),lastyears(i));
    [FLONG,FLAT]=meshgrid(Flong,Flat);
       
    fslope = (fslopeF);
    
    % Start mapping %%%
    clf;
    ax = worldmap([minlat maxlat],[minlong maxlong]);
    setm(ax,'meridianlabel','off','parallellabel','off','flinewidth',3);
    hold on;
            
    mapped = griddata(FLONG(:),FLAT(:),fslope(:),Flong,Flat(:),'linear'); 
    
    sdFmapped = griddata(FLONG(:),FLAT(:),fsdF(:),Flong,Flat(:),'linear');
    sdmapped = griddata(FLONG(:),FLAT(:),sdslopeF(:),Flong,Flat(:),'linear');
    u=sdmapped/mean(sdpriorslope);
    threshold=quantile(u(:),.75);
    subbad=find(u>threshold);
    mapped(subbad)=NaN;   
   
    hs1=scatterm(FLAT(:),FLONG(:),250,mapped(:),'filled','marker','s');
    
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(ax, land, 'FaceColor', [0.85 0.85 0.85]);
    colormap(jet);   
    caxis([-5 5]);
    axis tight;
    hcb=colorbar;
    box on;

    title({['Rate of Sea Level Change ' num2str(datefirst) ' to ' num2str(datelast) ' CE']});
    pdfwrite(['SL_Rate_' num2str(datefirst) '_to_' num2str(datelast)]);

    clf;
    ax = worldmap([minlat maxlat],[minlong maxlong]);
    setm(ax,'meridianlabel','off','parallellabel','off','flinewidth',3);
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(ax, land, 'FaceColor', [0.85 0.85 0.85]);
    hold on;

    hs1=scatterm(FLAT(:),FLONG(:),250,sdFmapped(:),'filled','marker','s');

    hold on;   %     
    geoshow(ax, land, 'FaceColor', [0.85 0.85 0.85]);
    colormap(parula);    
    caxis([0.2 1.4]);
    axis tight;
    hcb=colorbar;

    title({['SD of Sea Level Change ' num2str(datefirst) ' to ' num2str(datelast) ' CE']});
    pdfwrite(['SL_SD_' num2str(datefirst) '_to_' num2str(datelast)]);

end
