function makeplot_tg(dataset,f2s_shifted,sd2s,V2s,testloc,labl,doplots,difftimestep,xlim0,maxdistfrom,meanSL,plot_dat,GIAtime,GIAsite,GIArsl6G,GIArsl7G)
if xlim0==1
    maxdistfrom=1;
    maxerror=700;
    difftimestep=100;
    xl=[0 2010];
    yl=[-4e3 100];
    yl2=[-1.5 4];
    labl ='_Common_Era';
elseif xlim0==2
    maxdistfrom=1;
    maxerror=5000;
    difftimestep=200;
%     xl=[-13000 2010];
%     yl=[-80e3 1000];
%     yl2=[-10 20];
%     difftimestep=100;
      xl=[-10000 2010];
    yl=[-40e3 1000];
    yl2=[-10 20];
%    labl =['_Holocene_dt_' num2str(difftimestep)];
elseif xlim0==4
    maxdistfrom=1;
    maxerror=5000;
    difftimestep=20;
    xl=[1900 2010];
    yl=[-500 200];
    yl2=[-10 20];
%    labl =['_TG_' num2str(difftimestep)];
elseif xlim0==22
    maxdistfrom=1;
    maxerror=5000;
    difftimestep=500;
    xl=[-10000 2010];
    yl=[-1e3 1000];
    yl2=[-10 20];
elseif xlim0==3
    maxdistfrom=1;
    maxerror=5000;
    difftimestep=10;
    xl=[1900 2010];
    yl=[-500 200];
    yl2=[-10 20];
    labl ='_instrumtal_period';
end
minY=round(min(dataset.Y),-4);
maxY=round(max(dataset.Y),-4);
if size(GIArsl6G,1)>0
    GIAtime=1950-GIAtime;
    GIArsl6G=GIArsl6G/1000;
    GIArsl7G=GIArsl7G/1000;
else
    GIAsite=testloc.reg;
    GIArsl6G=testloc.GIAproj/1000;
    GIAtime=1950-testloc.X(:,3);
end

defval('labl','');
defval('doplots',[]);
defval('difftimestep',100);
defval('xlim0',[-10000 2000]);
defval('maxdistfrom',.1);
numrows=1+(difftimestep>0);

angd= @(Lat0,Long0,lat,long) (180/pi)*(atan2(sqrt((cosd(lat).*sind(long-Long0)).^2+(cosd(Lat0).*sind(lat)-sind(Lat0).*cosd(lat).*cosd(long-Long0)).^2),(sind(Lat0).*sind(lat)+cosd(Lat0).*cosd(lat).*cosd(long-Long0))));

dYears=@(years1,years2) abs(bsxfun(@minus,years1',years2));
dDist=@(x1,x2)angd(repmat(x1(:,1),1,size(x2,1)),repmat(x1(:,2),1,size(x2,1)),repmat(x2(:,1)',size(x1,1),1),repmat(x2(:,2)',size(x1,1),1))'+1e6*(bsxfun(@plus,x1(:,1)',x2(:,1))>1000);

testreg=testloc.reg;
testsites=testloc.sites;
testnames=testloc.names;
testnames2=testloc.names2;
testX(:,1:2)=testloc.X(:,1:2);
%testX(:,3)=1950-testloc.X(:,3);
testX(:,3)=testloc.X(:,3);

% change everything back to ages instead of common era dates
istg = dataset.istg;
if isfield(dataset,'GIAproj')
    GIAprojSL=dataset.GIAproj/1000;
end      
    
meantime=dataset.meantime;
%meantime=1950-dataset.meantime;
lat=dataset.lat;
long=dataset.long;
Y=dataset.Y; % original Y from the data before GIA subtracted
%+GIAprojSL)./1000;
%Y=(dataset.Y0+meanSL);
Ycv=dataset.Ycv;
limiting=dataset.limiting;
if max(dataset.datid)>3000
    indicator=dataset.indicator;
    %intercalated=dataset.intercalated;
    compactcorr=dataset.compactcorr;
    %         wtime1=1950-(datHolo.data(sub,6)+datHolo.data(sub,7)) + count/1e5;
    %         wtime2=1950-(datHolo.data(sub,6)-datHolo.data(sub,8)) + count/1e5;
    dY_init=dataset.dY_init;
else
    if plot_dat==1
        indicator=dataset.limiting+1;
    else
        indicator=dataset.limiting;
    end
end
time1=dataset.time1;
%time1 =1950-dataset.time1;
time2=dataset.time2;
%time2 =1950-dataset.time2;
dY=dataset.dY;
datid=dataset.datid;
% shift all sea levels back by the original mean of the data
 f2s=f2s_shifted;

for i=1:size(testsites,1)
	clf; clear hp;
    
%     hp(i)=subplot(numrows,length(js),k+length(js));
%     box on;
% 
%     hold on

	subA = find(testreg == testsites(i,1));
    subB=find(datid==testsites(i,1));    
    subG = find(GIAsite==testsites(i,1)); 
	clf; clear hp;

		offsetA=0;
		box on;
        hold on
        plot(testX(subA,3),f2s(subA)+offsetA,'Color',[.5, .5, .5]);
		plot(testX(subA,3),f2s(subA)+sd2s(subA)+offsetA,'--','Color',[.5, .5, .5]);
		plot(testX(subA,3),f2s(subA)-sd2s(subA)+offsetA,'--','Color',[.5, .5, .5]);
		plot(testX(subA,3),f2s(subA)+2*sd2s(subA)+offsetA,':','Color',[.5, .5, .5]);
		plot(testX(subA,3),f2s(subA)-2*sd2s(subA)+offsetA,':','Color',[.5, .5, .5]);
		hold on
    subB=find(datid==testsites(i,1));    
 plot(meantime(subB),Y(subB),'r');
        xlim(xl);
%        ylim(yl);
%         xlim([xlims(i,1) xlims(i,2)]);
%         ylim([ylim_SLs(i,1) ylim_SLs(i,2)]);
%         set(gca, 'YTick', [ylim_SLs(i,1):ticks_SL(i):ylim_SLs(i,2)]);
%        ylim(ylim_SL);
		ylabel('TG timeseries (mm)','Color','k');
 pdfwrite([testnames{i} '_dat_' labl]); 
end
    
