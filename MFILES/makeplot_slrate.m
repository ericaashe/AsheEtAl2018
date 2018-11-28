function [sl_rate, testnames] = makeplot_slrate(dataset,f2s_shifted,sd2s,V2s,testloc,labl,doplots,difftimestep,xlim0,meanSL,plot_dat,GIAtime,GIAsite,GIArsl6G,GIArsl7G,yl)
%% plotting used for Ashe et al., 2018

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
    xl=[-10000 2010];
    yl=[-40e3 1000];
    yl2=[-10 20];
elseif xlim0==4
    maxdistfrom=1;
    maxerror=5000;
    difftimestep=20;
    xl=[1900 2010];
    yl=[-500 200];
    yl2=[-10 20];
elseif xlim0==22
    maxdistfrom=1;
    maxerror=5000;
    difftimestep=1000;
    xl=[-10000 2000];
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
testX(:,3)=testloc.X(:,3);

istg = dataset.istg;
if isfield(dataset,'GIAproj')
    GIAprojSL=dataset.GIAproj/1000;
end      
    
meantime=dataset.meantime;
lat=dataset.lat;
long=dataset.long;
Y=dataset.Y; 
Ycv=dataset.Ycv;
limiting=dataset.limiting;
if max(dataset.datid)>3000
    indicator=dataset.indicator;
    compactcorr=dataset.compactcorr;
    dY_init=dataset.dY_init;
else
    if plot_dat==1
        indicator=dataset.limiting+1;
    else
        indicator=dataset.limiting;
    end
end
time1=dataset.time1;
time2=dataset.time2;
dY=dataset.dY;
datid=dataset.datid;
f2s=f2s_shifted;

for i=1:size(testsites,1)

	clf; clear hp;
	subA = find(testreg == testsites(i,1));
    subB=find(datid==testsites(i,1));    
    subG = find(GIAsite==testsites(i,1)); 
	js=1;        

	for k=1:length(js)
		j=js(1,k);
		hp(k)=subplot(numrows,length(js),k);
		box on;
		
		hold on
        
        if plot_dat==1
			for nn=1:length(subB)
                if limiting(subB(nn))==0
                    if (indicator(subB(nn))==1)
                        rectangle('Position',[min(time1(subB(nn)),time2(subB(nn))),Y(subB(nn))-2*dY(subB(nn)),dYears(time2(subB(nn)),time1(subB(nn))),4*dY(subB(nn))],'Edgecolor',[.42,.77,1]);  % index points basal, turquoise       
                    elseif indicator(subB(nn))==2 % uniformly distributed data
                        rectangle('Position',[min(time1(subB(nn)),time2(subB(nn))),Y(subB(nn))-2*dY_init(subB(nn)),dYears(time2(subB(nn)),time1(subB(nn))),4*dY_init(subB(nn))],'Edgecolor','k');         
                    end
                elseif limiting(subB(nn))==1  % fresh water limiting points, sea green
                    plot([1 1]*meantime(subB(nn)),Y(subB(nn))+[-2 2]*dY(subB(nn)),'Color',[.05,.7,.5]);
                    plot([time1(subB(nn)) time2(subB(nn))],[1 1]*(Y(subB(nn))-2*dY(subB(nn))),'Color',[.05,.7,.5]);
                elseif limiting(subB(nn))==-1% marine limiting points, navy
                    plot([1 1]*meantime(subB(nn)),Y(subB(nn))+[-2 2]*dY(subB(nn)),'Color',[.1,.35,.7]);
                    plot([time1(subB(nn)) time2(subB(nn))],[1 1]*(Y(subB(nn))+2*dY(subB(nn))),'Color',[.1,.35,.7]);
                end
            end
        end
    end
        hold on
        plot(testX(subA,3),f2s(subA),'Color',[.5, .5, .5]);
		plot(testX(subA,3),f2s(subA)+sd2s(subA),'--','Color',[.5, .5, .5]);
		plot(testX(subA,3),f2s(subA)-sd2s(subA),'--','Color',[.5, .5, .5]);
		plot(testX(subA,3),f2s(subA)+2*sd2s(subA),':','Color',[.5, .5, .5]);
		plot(testX(subA,3),f2s(subA)-2*sd2s(subA),':','Color',[.5, .5, .5]);  
        xlim(xl);
        ylim(yl);
		ylabel('RSL Height (mm)','Color','k');
	
    if difftimestep>0
        Mdiff = bsxfun(@eq,testX(:,3),testX(:,3)')-bsxfun(@eq,testX(:,3),testX(:,3)'+ difftimestep);
        Mdiff = Mdiff .* bsxfun(@eq,testreg,testreg');
        sub=find(sum(Mdiff,2)==0);
        Mdiff=Mdiff(sub,:);
        difftimes=bsxfun(@rdivide,abs(Mdiff)*testX(:,3),sum(abs(Mdiff),2));
        diffreg=bsxfun(@rdivide,abs(Mdiff)*testreg,sum(abs(Mdiff),2));
        Mdiff=bsxfun(@rdivide,Mdiff,Mdiff*testX(:,3));
        subA2 = find(diffreg == testsites(i,1));
        
        clear df2s dV2s dsd2s;
        for n=1:size(f2s,2)
            df2s(:,n)=Mdiff*f2s(:,n);
            dV2s(:,:,n)=Mdiff*V2s(:,:,n)*Mdiff';
            dsd2s(:,n)=sqrt(diag(dV2s(:,:,n)));
        end

        for k=1:length(js)
            hp(k+length(js))=subplot(numrows,length(js),k+length(js));
            box on;

            hold on
            j=js(1,k);
            plot(difftimes(subA2),df2s(subA2),'Color',[.5, .5, .5]);
            plot(difftimes(subA2),df2s(subA2)+dsd2s(subA2),'--','Color',[.5, .5, .5]);
            plot(difftimes(subA2),df2s(subA2)-dsd2s(subA2),'--','Color',[.5, .5, .5]);
            plot(difftimes(subA2),df2s(subA2)+2*dsd2s(subA2),':','Color',[.5, .5, .5]);
            plot(difftimes(subA2),df2s(subA2)-2*dsd2s(subA2),':','Color',[.5, .5, .5]);

            xlim(xl);  
            ylabel('Rate (mm/yr = m/ka)','Color','k');
        end
    end
    title(hp(1),[testnames{i}])

    if plot_dat==1
        pdfwrite([testnames{i} '_with_dat_' labl]); 
    else
        pdfwrite([testnames{i} '_' labl]); 
    end
end

