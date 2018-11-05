%% taking out references to GIA 10-11-17
%     minlat=min(lat)-2;
%     maxlat=max(lat)+2;
%     minlong=min(mod(long,360))-2;
%     maxlong=max(mod(long,360))+2;
% for the Florida model, although we are using more data, we are only
% plotting around Florida
minlat=17.0300;
maxlat=32.7100;
minlong=270.5100;
maxlong=284.8900;
    [minlat maxlat minlong maxlong]
    lat_incr = .2;
    long_incr = .2;
    %     lat_incr = round((maxlat-minlat)/9)/10;
    %     long_incr = round((maxlong-minlong)/7)/10;
        Flat=minlat:lat_incr:maxlat;
        Flong=minlong:long_incr:maxlong;

%labl = 'Full';
labl = 'Local';
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
firstyears=[-4000:1000:1000];% 1000:200:1800 1900 ]; 
lastyears=[-3000:1000:2000];% 1100:200:2000 2000];
% firstyears=[-8050:1000:950]; 
% lastyears=[-7050:1000:1950];
% GIA=[];
% GIA.Lat=GIALat;
% GIA.Lon=mod(GIALon,360);
% GIA.GIA=GIA6;
% %GIA.times=GIAtime;
% GIA.SiteRate=siteGIArate;

for i=1:length(firstyears)
    figure;
    nsitepts{i}=[];
    td = lastyears(i)-firstyears(i);
    addt=(2000-td)/2;
    fy=firstyears(i);%-addt;
    ly=lastyears(i);%+addt;
    subt = intersect(find(wdataset.time1<=ly),find(wdataset.time2>=fy));
    %find the number of sites with data influencing the model close to this time period
    sub=find(testsites(:,2)<=360);
%     ulat=unique(round(testsites(sub,2)));
%     ulong=unique(round(testsites(sub,3)));
%     ulat=unique(bsxfun(@plus,ulat,[-4:2:4]));
%     ulong=unique(mod(bsxfun(@plus,ulong,[-4:2:4]),360));
%     sub=find(abs(ulat)<90); ulat=ulat(sub); ulat=ulat(:)'; ulong=ulong(:)';
%     Flat=union(Flat,ulat);
%     Flong=union(mod(Flong,360),mod(ulong,360));
    wmodelspec = modelspec(1);

    [fslopeFs,sdslopeF,fsF,sdsF,~,~,~,~,passderivs,invcv,fmeanFs,fsdF] = RegressRateField_noGIA_ea(wdataset,wmodelspec,thetTGG{jj},noiseMasks(1,:),Flat,Flong,firstyears(i),lastyears(i),trainsub,[],[],[]);
    [priorslope,sdpriorslope,priorfsF,priorsdsF] = RegressRateField_noGIA_ea(nulldataset,wmodelspec,thetTGG{jj},noiseMasks(1,:),Flat,Flong,firstyears(i),lastyears(i),trainsub,[],[],[]);

    
%(wdataset,wmodelspec,thetTGG{jj},noiseMasks(2,:),Flat,Flong,firstyears(i),lastyears(i),trainsub,[]);
%    [priorslope,sdpriorslope,priorfsF,priorsdsF] = RegressRateField_from_samples(nulldataset,wmodelspec,thetTGG{jj},noiseMasks(1,:),mean(Flat),mean(Flong),firstyears(i),lastyears(i),[],[]);
% N=zeros(0.5*length(fGIAs),length(fGIAs));
% for pp=1:0.5*length(fGIAs)
%     N(pp,2*pp-1:2*pp)=[0.5 0.5];
% end
% meanGIA=N*fGIAs;

     fmeanF = (fmeanFs)./1000;% + meanGIA;
%     fmeanSL = meanSL./1000;
     fmeanSL = 0;
     fslopeF  = fslopeFs./1000;
     fsdF = fsdF./1000;
     sdslope = sdslopeF./1000;
%    fmeanF = (fmeanFs);% + meanGIA;
%    fmeanSL = meanSL;
%    fsdF = fsdF;
    % Start mapping %%%

    clf;
    ax = worldmap([minlat-1 maxlat+1],[minlong-1 maxlong+1]);
    setm(ax,'meridianlabel','off','parallellabel','off','flinewidth',3);
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(ax, land, 'FaceColor', [0.85 0.85 0.85]);
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

  %  hs1=scatterm(FLAT(:),FLONG(:),40,mapped(:),'filled','marker','s');
  %hs1=contourm(FLAT(:),FLONG(:),mapped(:));
 % ,'filled','marker','s');
%  hs1=scatterm(FLAT(:),FLONG(:),60,mapped(:),'filled','marker','s');
  % hs1=scatterm(FLAT(:),FLONG(:),30,fmeanF(:),'filled','marker','s');
    hs1=scatterm(FLAT(:),FLONG(:),250,mapped(:),'filled','marker','s');


%%% threshold of these points is over the limit
%    hs2=plotm(FLAT(subbad),FLONG(subbad),'o','color',[0.99 .99 .99],'markerfacecolor', [.99 .99 .99],'markeredgecolor',[.99 .99 .99],'markersize',12);
%    [c,h] = contourm(Flat,Flong,mapped,[0],'-','Color',[0.99 0.99 0.99]);

            %hs2=plotm(FLAT1(subbad),FLONG1(subbad),'o','color',[1 1 1],'markersize',4,'markerfacecolor','w','markeredgecolor','w')
    hold on;
    sublong=find((mod(FLONG1(:),5)==0).*(mod(FLAT1(:),5)==0));
    %subbad2=intersect(subbad,sublong);
    %hbad=plotm(FLAT1(subbad2),FLONG1(subbad2),'kx');
    %set(hbad,'MarkerSize',3,'MarkerEdgeColor',[.5 .5 .5],'MarkerFaceColor',[.5 .5 .5]);
    %hbad=scatterm(FLAT1(subbad2),FLONG1(subbad2),10,[.5 .5 .5],'x');

    geoshow(ax, land, 'FaceColor', [0.85 0.85 0.85]);
    hold on;
    
    box on;
    
    subt = intersect(find(wdataset.time1<=ly),find(wdataset.time2>=fy));

%     ud=unique(wdataset.datid(find((wdataset.time2>=firstyears(i)).*(wdataset.time1<=lastyears(i)))));
%     sub1=find(ismember(wdataset.siteid,ud));
    %hs1=scatterm(wdataset.sitecoords(sub1,1),wdataset.sitecoords(sub1,2),15,'k','filled','MarkerFaceColor','w','Marker','d','MarkerEdgeColor','k'); hold on;
    % map the regions where there is data and how much
    for j=1:length(wdataset.siteid)
      %is member?
      nsitepts{i}(j)=length(intersect(intersect(intersect(find(wdataset.limiting==0),find(wdataset.time1<=ly)),find(wdataset.time2>=fy)),find(wdataset.datid==wdataset.siteid(j))));
    end        

    subsite = find(nsitepts{i}>0);
    ndat = nsitepts{i}(subsite)';
    lat_site = wdataset.sitecoords(subsite,1);
    long_site = wdataset.sitecoords(subsite,2);

    uns = unique(ndat);
    
%scatterm(wdataset.sitecoords(lim_data_sub,1),wdataset.sitecoords(lim_data_sub,2),2.5*ndatapoints(lim_data_sub),'MarkerEdgeColor','k','MarkerFaceColor',[.7 .7 .7]);
 box on;   
 %% plotting the data
%             ss_1_3=intersect(find(ndat>=1),find(ndat<=4));
%             ss_4_10=intersect(find(ndat>=5),find(ndat<=15));
%             ss_11_20=intersect(find(ndat>=16),find(ndat<=500));
%            ss_21=find(ndat>=21);
% 
%            hgh5=scatterm(lat_site(ss_21),long_site(ss_21),2*ndat(ss_21),'MarkerEdgeColor','k','MarkerFaceColor','w');
%            hgh2=scatterm(lat_site(ss_11_20),long_site(ss_11_20),2*ndat(ss_11_20),'MarkerEdgeColor','k','MarkerFaceColor','w');
%            hgh4=scatterm(lat_site(ss_4_10),long_site(ss_4_10),7*ndat(ss_4_10),'MarkerEdgeColor','k','MarkerFaceColor','w');
%            hgh3=scatterm(lat_site(ss_1_3),long_site(ss_1_3),8*ndat(ss_1_3),'MarkerEdgeColor','k','MarkerFaceColor','w');
%             %hgh3=
%             %hgh2=
%         %     tide_x=find(wdataset.istg==1);
%         %     tgs=unique(wdataset.datid(tide_x));
%         %     tgs=find(testsitedef.sites(:,1)<30000);
%         %     ccs=find(testsitedef.sites(:,1)<30000);
% 
% %         cc=find(wdataset.cont==1);
% %         tg=find(wdataset.istg==1);
% %         di=intersect(find(wdataset.istg==0),find(wdataset.cont==0));
% % %        subt=within this time;
% %     cct=intersect(cc,subt);
% % %    tgt=intersect(tg,subt);
% %     dit=intersect(di,subt);
% %     scatterm(wdataset.lat(dit),wdataset.long(dit),'Marker','x','MarkerEdgeColor','k');%,'MarkerEdgeColor','k','MarkerFaceColor','b');
% %     scatterm(wdataset.lat(cct),wdataset.long(cct),'Marker','d','MarkerEdgeColor','k','MarkerFaceColor','b');
% % %    scatterm(wdataset.lat(tgt),wdataset.long(tgt),'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','w');
% %         %     scatterm(testsitedef.sites(ccs,2),testsitedef.sites(ccs,3),40,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor','b');
% %         %     scatterm(testsitedef.sites(ss_1_3,2),testsitedef.sites(ss_1_3,3),60,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r');
% %         %     scatterm(testsitedef.sites(tgs,2),testsitedef.sites(tgs,3),40,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor','w');
% %          %   hgh4=scatterm(lat_site(ss_1_3),long_site(ss_1_3),3*ndat(ss_1_3),'MarkerEdgeColor','k','MarkerFaceColor','w');
% % 
% % 
%         %  ud=unique(wdataset.datid(find((wdataset.time2>=firstyears(i)).*(wdataset.time1<=lastyears(i)))));
%         % sub1=find(ismember(wdataset.siteid,ud));
%         %
%         %  hs1=scatterm(wdataset.sitecoords(sub1,1),wdataset.sitecoords(sub1,2),20,'k','filled','MarkerFaceColor','w','Marker','d','MarkerEdgeColor','k'); hold on;
% 
% %             hg1=plot(nan,nan,'o','MarkerEdgeColor','k','MarkerFaceColor','w');
% %              hg3=plot(nan,nan,'d','MarkerEdgeColor','k','MarkerFaceColor','b');
% %              hg2=plot(nan,nan,'x','MarkerEdgeColor','k');
% 
% %            legend([hg2,hg3],'Discrete Index Points','Continuous Cores','Location','Southeast');
% %            legend([hg1,hg2,hg3],'Tide Gauges','Discrete Index Points','Continuous Cores','Location','Southeast');
% %             legend([hg1,hg2,hg3],'Tide Gauges','Discrete Index Points','Continuous Cores','Location','Southeast');
%         %     hcbpos=get(hcb,'Position');
%         %     htt = annotation('textbox');
%         %     htt.String=sprintf('GSL: %0.1f \\pm %0.1f mm/yr',[fslopeGSL 2*sdslopeGSL]);
%         %     htt.FontSize=13;
%         %     htt.BackgroundColor=[1 1 1];
%         %     htt.Margin=1;
%         %     htt.VerticalAlignment='bottom';
%         %     pstn=htt.Position;
%         %     pstn(1)=[.30];
%         %     pstn(2)=hcbpos(2)+.005;
%         %     htt.Position=pstn;
% 
% %     s = get(lhgh);
% %     s2=findobj(s.Children,{'type','patch'});
% %    set(s2,'MarkerSize',sqrt(200));
% %    set(s3,'MarkerSize',sqrt(85));
% %   set(s2(3),'MarkerSize',sqrt(38));
% %   set(s2(4),'MarkerSize',sqrt(12));
% 
% %    scatterm(lat_site,long_site,25*ndat,'MarkerEdgeColor','k','MarkerFaceColor','w');
% %colormap(brewermap(11,'RdYlBu'));
% %    colormap(jet);
    axis tight;
    hcb=colorbar;

    box on;
% %    caxis([-32 -12]);
% %   caxis([-25 10]);
% %      if i < 4
% %         caxis([-32 -12]);
% %      elseif i < 7
% %          caxis([-18 -5]);
% %      elseif i < 10
% %          caxis([-8 -2]);
% %      else
% %          caxis([-5 0]);
% %      end
%     if i < 2
%         caxis([-32 -26]);
%      elseif i < 3
%          caxis([-12 -7]);
%      elseif i < 4
%          caxis([-3.2 -1]);
%      else
%          caxis([-5 0]);
%      end
%      %caxis([-2000 100]);
% %     elseif firstyears(i)>1900
% % 
% %     end
% %         caxis([-100 250]);
% %     else
% %         caxis([-50 100]);
% %     end
% %    caxis([-50 250]);
% %   caxis([-25 10]);
% %     if i < 6
% %         caxis([-150 750]);
% %     elseif i < 10
% %         caxis([-150 500]);
% %     elseif i < 13
% %         caxis([-100 250]);
% %     else
% %         caxis([-50 100]);
% %     end
% %    caxis([-50 250]);
      box on;

       datehs = (1950-firstyears(i)+1950-lastyears(i))/2;
       datefirst = 1950-firstyears(i);
       datelast = 1950-lastyears(i);

%        title({['Sea Level Height ' num2str(firstyears(i)) ' to ' num2str(lastyears(i)) ' CE']});
%            %num2str(datefirst) ' to ' num2str(datelast) ' BP']});
%        pdfwrite(['SL_Height' labl '_' num2str(firstyears(i)) '_to_' num2str(lastyears(i))]);

       title({['SL Height' num2str(firstyears(i)) ' to ' num2str(lastyears(i)) ' CE']});
           %num2str(datefirst) ' to ' num2str(datelast) ' BP']});
       pdfwrite(['Height_Sea_Level' labl '_' num2str(firstyears(i)) '_to_' num2str(lastyears(i)) '_no_data']);

%     hsr1=scatterm(FLAT(:),FLONG(:),250,mapRate(:),'filled','marker','s');
%     hsr2=plotm(FLAT(subbad),FLONG(subbad),'o','color',[0.99 .99 .99],'markerfacecolor', [.99 .99 .99],'markeredgecolor',[.99 .99 .99],'markersize',12);
%     [c,h] = contourm(Flat,Flong,mapRate,[0],'-','Color',[0.99 0.99 0.99]);
%     hold on;
%     sublong=find((mod(FLONG1(:),5)==0).*(mod(FLAT1(:),5)==0));
% 
%     geoshow(ax, land, 'FaceColor', [0.85 0.85 0.85]);
%     hold on;
%     box on;
%     subt = intersect(find(wdataset.time1<=ly),find(wdataset.time2>=fy));
% %     ud=unique(wdataset.datid(find((wdataset.time2>=firstyears(i)).*(wdataset.time1<=lastyears(i)))));
% %     sub1=find(ismember(wdataset.siteid,ud));
%     %hs1=scatterm(wdataset.sitecoords(sub1,1),wdataset.sitecoords(sub1,2),15,'k','filled','MarkerFaceColor','w','Marker','d','MarkerEdgeColor','k'); hold on;
%     % map the regions where there is data and how much
%     for j=1:length(wdataset.siteid)
%       %is member?
%       nsitepts{i}(j)=length(intersect(intersect(intersect(find(wdataset.limiting==0),find(wdataset.time1<=ly)),find(wdataset.time2>=fy)),find(wdataset.datid==wdataset.siteid(j))));
%     end        
% 
%     subsite = find(nsitepts{i}>0);
%     ndat = nsitepts{i}(subsite)';
%     lat_site = wdataset.sitecoords(subsite,1);
%     long_site = wdataset.sitecoords(subsite,2);
%     uns = unique(ndat);
%     box on;   
%     axis tight;
% 
%     if i < 2
%         caxis([6.8 7.7 ]);
%     elseif i < 3
%         caxis([2 3.5 ]);
%     elseif i<4
%         caxis([0.5 2 ]);        
%     else
%         caxis([0 15]);
%     end
%         
%     hcb=colorbar;
%     box on;
%        datehs = (1950-firstyears(i)+1950-lastyears(i))/2;
%        datefirst = 1950-firstyears(i);
%        datelast = 1950-lastyears(i);
%        title({['Rate of RSL Change (m/ky)' num2str(firstyears(i)) ' to ' num2str(lastyears(i)) ' CE']});
%            %num2str(datefirst) ' to ' num2str(datelast) ' BP']});
%        pdfwrite(['Rate_RSL_Change' labl '_' num2str(firstyears(i)) '_to_' num2str(lastyears(i)) '_scaled']);   
%        
    figure;
    ax = worldmap([minlat-1 maxlat+1],[minlong-1 maxlong+1]);
    setm(ax,'meridianlabel','off','parallellabel','off','flinewidth',3);
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(ax, land, 'FaceColor', [0.85 0.85 0.85]);
    hold on;

    hs1=scatterm(FLAT(:),FLONG(:),250,sdmapped(:),'filled','marker','s');

    hold on;   %     
    geoshow(ax, land, 'FaceColor', [0.85 0.85 0.85]);
    hold on;
% scatterm(wdataset.sitecoords(lim_data_sub,1),wdataset.sitecoords(lim_data_sub,2),2.5*ndatapoints(lim_data_sub),'MarkerEdgeColor','k','MarkerFaceColor',[.7 .7 .7]);
    %scatterm(lat_site,long_site,25*ndat,'MarkerEdgeColor','k','MarkerFaceColor','w');
%    colormap(jet);
    axis tight;
    hcb=colorbar;

    box on;
      caxis([0 4.5]);
       box on;

    title({['SD of Sea Level Change ' num2str(firstyears(i)) ' to ' num2str(lastyears(i)) ' BP']});
    pdfwrite(['SL_SD_' num2str(firstyears(i)) '_to_' num2str(lastyears(i))]);

end
