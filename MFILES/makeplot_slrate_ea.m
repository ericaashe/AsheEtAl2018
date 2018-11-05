function [sl_rate, testnames] = makeplot_slrate_ea(dataset,f2s_shifted,sd2s,V2s,testloc,labl,doplots,difftimestep,xlim0,maxdistfrom,meanSL,plot_dat,GIAtime,GIAsite,GIArsl6G,GIArsl7G)
%% plotting used for the Review paper?
% changed tie-in Barbados early Holocene tie-in to earlier, which changes
% the curves toward the beginning of the plots (3-8-18)

%ylim_SL,ytick_SL,ytick_lab,ylim_Rate,ytick_Rate,ytickRateLab,thetPXholo)
% GIA for each f2s needs to be added back in
% if size(limits,1)>1
%     ylim_SLs=limits(:,1:2);
%     xlims = limits(:,3:4)*1000;
%     ylim_rates = limits(:,5:6);
%     lab_i=limits(:,9);
%     labels{1} = {'20','10','0','-10','-20'};
%     labels{2} = {'15','10','5','0','-5','-10','-15','-20'};
%     labels{3} = {'10','5','0','-5','-10','-15','-20'};
%     labels{4} = {'0','-20','-40','-60','-80'};
%     labels{5} = {'10','0','-10','-20','-30','-40','-50'};
%     ticks = limits(:,7);
%     ticks_SL=limits(:,8);
% else
%     xlims = repmat(xlim0,size(testloc.sites,1),1);
%     defval('ylim_SL',[-15 2]);
% %     defval('ytick_SL',[-60e3:20e3:60e3]);
% %     defval('ytick_lab',[-60:20:60]);
% end
% 

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
%     xl=[-13000 2010];
%     yl=[-80e3 1000];
%     yl2=[-10 20];
%    difftimestep=100;
    xl=[-10000 2010];
    yl=[-1e3 1000];
    yl2=[-10 20];
%    labl =['_Holo_dt_' num2str(difftimestep)];
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
    % defval('ylim_SL',[minY maxY]);
% defval('ytick_SL',[minY:50e3:maxY]);
% defval('ytick_lab',[minY/1000:50:maxY/1000]);
% % % 
% defval('ylim_Rate',[-20 30]);
% defval('ytick_Rate',[-20:10:30]);
% defval('ytickRateLab',{'20','10','0','-10','-20','-30'});
% defval('ylim_Rate',[-20 5]);
% defval('ytick_Rate',[-20:5:5]);
% %defval('ytickRateLab',{'5','0','-5','-10','-15','-20'});%             ylim([-30 40]);
% defval('ytickRateLab',{'20','15','10','5','0','-5'});%             ylim([-30 40]);

%defval('ytickRateLab',{'20','15','10','5','0','-5','-10','-15'});

%              set(gca,'YTick',[-25:5:5]);
%              set(gca,'YDir','reverse');
%              set(gca,'YTickLabel',{'25','20','15','10','5','0','-5'});


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

% testreg=testlocs{1}.reg;
% testsites=testlocs{1}.sites;
% testnames=testlocs{1}.names;
% testnames2=testlocs{1}.names2;
% testX=testlocs{1}.X;

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
 %{1}./1000+GIArsl6G;
% sd2s=sd2s{1}/1000;
% V2s=V2s{1}/1e6;

%figure;

for i=1:size(testsites,1)
    %[36 35 43]
    %
    %
%30:size(testsites,1) % = 46 sites
%i=5;

	clf; clear hp;
    
%     hp(i)=subplot(numrows,length(js),k+length(js));
%     box on;
% 
%     hold on

	subA = find(testreg == testsites(i,1));
    subB=find(datid==testsites(i,1));    
    subG = find(GIAsite==testsites(i,1)); 
	clf; clear hp;
	
js=1;        

	for k=1:length(js)
		offsetA=0;
		j=js(1,k);
		hp(k)=subplot(numrows,length(js),k);
		box on;
		
		hold on
if plot_dat==1
        if k==1
			for nn=1:length(subB)
                if limiting(subB(nn))==0
                    if (indicator(subB(nn))==1)
%                         if (intercalated(subB(nn))==0)
                             rectangle('Position',[min(time1(subB(nn)),time2(subB(nn))),Y(subB(nn))-2*dY(subB(nn)),dYears(time2(subB(nn)),time1(subB(nn))),4*dY(subB(nn))],'Edgecolor',[.42,.77,1]);  % index points basal, turquoise       
                            % plot([1 1]*meantime(subB(nn)),Y(subB(nn))+[-1 1]*dY(subB(nn)),'b');
                            % plot([time1(subB(nn)) time2(subB(nn))],Y(subB(nn))*[1 1],'b');
                            % plot([1 1]*meantime(subB(nn)),Y(subB(nn)),'b.');
%                         elseif (intercalated(subB(nn))==1)
%                             rectangle('Position',[min(time1(subB(nn)),time2(subB(nn))),Y(subB(nn))-2*dY(subB(nn)),dYears(time2(subB(nn)),time1(subB(nn))),4*dY(subB(nn))],'Edgecolor','r');         
%                             % plot([1 1]*meantime(subB(nn)),Y(subB(nn))+[-1 1]*dY(subB(nn)),'g');
%                             % plot([time1(subB(nn)) time2(subB(nn))],Y(subB(nn))*[1 1],'g');
%                             % plot([1 1]*meantime(subB(nn)),Y(subB(nn)),'g.');
%                         end
                    elseif indicator(subB(nn))==2 % uniformly distributed data
                        % rectangle('Position',[x,y,w,h]) draws the rectangle from the point x,y and having a width of w and a height of h
                        % plot([1 1]*meantime(subB(nn)),Y(subB(nn))+[-1 1]*dY(subB(nn)),'b');
                        % plot([time1(subB(nn)) time2(subB(nn))],Y(subB(nn))*[1 1],'b');
                        % plot([1 1]*meantime(subB(nn)),Y(subB(nn)),'b.');
                        rectangle('Position',[min(time1(subB(nn)),time2(subB(nn))),Y(subB(nn))-2*dY_init(subB(nn)),dYears(time2(subB(nn)),time1(subB(nn))),4*dY_init(subB(nn))],'Edgecolor','k');         
                    end
                elseif limiting(subB(nn))==1  % fresh water limiting points, sea green
                    plot([1 1]*meantime(subB(nn)),Y(subB(nn))+[-2 2]*dY(subB(nn)),'Color',[.05,.7,.5]);
                    plot([time1(subB(nn)) time2(subB(nn))],[1 1]*(Y(subB(nn))-2*dY(subB(nn))),'Color',[.05,.7,.5]);
                    % plot([1 1]*meantime(subB(nn)),Y(subB(nn)),'bv');
                elseif limiting(subB(nn))==-1% marine limiting points, navy
                    plot([1 1]*meantime(subB(nn)),Y(subB(nn))+[-2 2]*dY(subB(nn)),'Color',[.1,.35,.7]);
                    plot([time1(subB(nn)) time2(subB(nn))],[1 1]*(Y(subB(nn))+2*dY(subB(nn))),'Color',[.1,.35,.7]);
                    % plot([1 1]*meantime(subB(nn)),Y(subB(nn))+[-1 1]*dY(subB(nn)));
                    % plot([time1(subB(nn)) time2(subB(nn))],Y(subB(nn))*[1 1]);
                    % plot([1 1]*meantime(subB(nn)),Y(subB(nn)),'^');
                end
            end
        end
end
        hold on
        plot(testX(subA,3),f2s(subA)+offsetA,'Color',[.5, .5, .5]);
		plot(testX(subA,3),f2s(subA)+sd2s(subA)+offsetA,'--','Color',[.5, .5, .5]);
		plot(testX(subA,3),f2s(subA)-sd2s(subA)+offsetA,'--','Color',[.5, .5, .5]);
		plot(testX(subA,3),f2s(subA)+2*sd2s(subA)+offsetA,':','Color',[.5, .5, .5]);
		plot(testX(subA,3),f2s(subA)-2*sd2s(subA)+offsetA,':','Color',[.5, .5, .5]);

        %plot(GIAtime(subG),GIArsl6G(subG),'--','Color','r');
        % plot(GIAtime(subG),GIArsl7G(subG),'--','Color','b');

		if j==3
                    %			title('Regional + Local');
		elseif j==4
                    %	title('Regional + Local non-linear');
		elseif j==6
                    %	title('Greenland');
        end
        
        xlim(xl);
%        ylim(yl);
%         xlim([xlims(i,1) xlims(i,2)]);
%         ylim([ylim_SLs(i,1) ylim_SLs(i,2)]);
%         set(gca, 'YTick', [ylim_SLs(i,1):ticks_SL(i):ylim_SLs(i,2)]);
%        ylim(ylim_SL);
		ylabel('RSL Height (mm)','Color','k');
% %        ylim([-40e3 200e3]);
%         ylim(ylim_SL);
%         set(gca,'YTick',ytick_SL);
%         %set(gca,'YTick',[-40000:10000:10000]);
%         set(gca,'YTickLabel',ytick_lab);
%        set(gca,'YTickLabel',{'-40','-20','0','20','40','80','120','160','200',});
        %set(gca,'YTickLabel',{'-35','-30','-25','-20','-15','-10','-5','0','5','10',})
%
    end
	
    if difftimestep>0

        Mdiff = bsxfun(@eq,testX(:,3),testX(:,3)')-bsxfun(@eq,testX(:,3),testX(:,3)'+ difftimestep);
        Mdiff = Mdiff .* bsxfun(@eq,testreg,testreg');
        sub=find(sum(Mdiff,2)==0);
        Mdiff=Mdiff(sub,:);
        difftimes=bsxfun(@rdivide,abs(Mdiff)*testX(:,3),sum(abs(Mdiff),2));
        diffreg=bsxfun(@rdivide,abs(Mdiff)*testreg,sum(abs(Mdiff),2));
        Mdiff=bsxfun(@rdivide,Mdiff,Mdiff*testX(:,3));
    subA2 = find(diffreg == testsites(i,1));
        
%         Mdiff = bsxfun(@eq,testX(:,3),testX(:,3)')-bsxfun(@eq,testX(:,3),testX(:,3)'+ difftimestep);
%         Mdiff = Mdiff .* bsxfun(@eq,testreg,testreg');
%         sub=find(sum(Mdiff,2)==0);
%         Mdiff=Mdiff(sub,:);
%         difftimes=bsxfun(@rdivide,abs(Mdiff)*testX(:,3),sum(abs(Mdiff),2));
%         diffreg=bsxfun(@rdivide,abs(Mdiff)*testreg,sum(abs(Mdiff),2));
%         Mdiff=bsxfun(@rdivide,Mdiff,Mdiff*testX(:,3));

        clear df2s dV2s dsd2s;
        for n=1:size(f2s,2)
%             df2s(:,n)=Mdiff*f2s(:,n);
%             dV2s(:,:,n)=Mdiff*V2s(:,:,n)*Mdiff';
%             dsd2s(:,n)=sqrt(diag(dV2s(:,:,n)));
            df2s(:,n)=Mdiff*f2s(:,n);
            dV2s(:,:,n)=Mdiff*V2s(:,:,n)*Mdiff';
            dsd2s(:,n)=sqrt(diag(dV2s(:,:,n)));
        end


        for k=1:length(js)
            offsetA=0;
            hp(k+length(js))=subplot(numrows,length(js),k+length(js));
            box on;

            hold on
            j=js(1,k);
            plot(difftimes(subA2),df2s(subA2),'Color',[.5, .5, .5]);
            plot(difftimes(subA2),df2s(subA2)+dsd2s(subA2),'--','Color',[.5, .5, .5]);
            plot(difftimes(subA2),df2s(subA2)-dsd2s(subA2),'--','Color',[.5, .5, .5]);
            plot(difftimes(subA2),df2s(subA2)+2*dsd2s(subA2),':','Color',[.5, .5, .5]);
            plot(difftimes(subA2),df2s(subA2)-2*dsd2s(subA2),':','Color',[.5, .5, .5]);

% %              ylim(ylim_Rate);
% %              set(gca,'YTick',ytick_Rate);
xlim(xl);  
%ylim(yl2);

%         xlim([xlims(i,1) xlims(i,2)]);
%         ylim([ylim_rates(i,1) ylim_rates(i,2)]);
%         set(gca, 'YTick', [ylim_rates(i,1):ticks(i):ylim_rates(i,2)]);
%         set(gca, 'YTickLabel', labels{lab_i(i)});
%         set(gca,'YDir','reverse');
%                 set(gca, 'YTickLabel', labels{lab_i(i)});

% %              set(gca,'YTickLabel',ytickRateLab);
%             ylim([-30 40]);
%              set(gca,'YTick',[-20:10:20]);
%              set(gca,'YDir','reverse');
%              set(gca,'YTick',[-25:5:5]);
%              set(gca,'YDir','reverse');
%              set(gca,'YTickLabel',{'25','20','15','10','5','0','-5'});
             ylabel('Rate (mm/yr = m/ka)','Color','k');
    %,'-10','-15','-20'});
            %set(gca,'YTickLabel',{'30','20','10','0','-10','-20','-30','-40'});
%            ylim([-20 5]);
        end
    end

    title(hp(1),[testnames{i}])
%    longticks(hp);
    try; 
    [bh,th]=label(hp,'ul',12,[],0,1,1,1.5,1.5); 
    end

if plot_dat==1
    pdfwrite([testnames{i} '_with_dat_' labl]); 
else
    pdfwrite([testnames{i} '_' labl]); 
end
%    sl_rate{i}=[difftimes(subA2),df2s(subA2),dsd2s(subA2)];
end

