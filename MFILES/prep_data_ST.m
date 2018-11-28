% prep_data_ST prepares spatio-temporal data for analysis

datid=[]; time1=[]; time2=[]; limiting=[]; Y=[]; dY = []; compactcorr = [];
istg = []; lat=[]; long=[]; sitelen=[]; indic = []; limiting = []; intercalated=[];
siteid=[]; sitenames={}; sitecoords=[]; overburden = []; %depthbase = []; thick = []; 
minY = []; dY_init=[]; compactrun=[]; GIA=[]; cont=[];

HoloRegions=unique(datPX.data(:,11));
Holositelat=[]; Holositelong=[]; Holositename={};
for curreg=1:length(HoloRegions)
	sub=find(datPX.data(:,11)==HoloRegions(curreg));
        count=[1:length(sub)]';
        wdatid = ones(length(sub),1)*(idHolo+HoloRegions(curreg)*1e3);
        wtime1=1950-(datPX.data(sub,6)+abs(datPX.data(sub,7))) + count/1e5;
        wtime2=1950-(datPX.data(sub,6)-abs(datPX.data(sub,8))) + count/1e5;
        wlimiting=(datPX.data(sub,9));
        windic=(datPX.data(sub,10));
        wintercalated=windic*0;
        wcompactcorr=wintercalated;
        woverburden=wintercalated;
        wY=datPX.data(sub,3)*1000;
        wdY=(datPX.data(sub,4)+abs(datPX.data(sub,5)))/4*1000;
        wminY = wY - 2*wdY;
        wlat=datPX.data(sub,1);
        wlong=mod(datPX.data(sub,2),360);
        wlat=wlat+1e-8*curreg;
        wcont=datPX.data(sub,12)*2;

                
        datid = [datid ; wdatid ; (idHolo+HoloRegions(curreg)*1e3)];
        time1 = [time1 ; wtime1 ; 2016 ];
        time2 = [time2 ; wtime2 ; 2017 ];
        limiting = [limiting ; wlimiting ; 0];
        intercalated = [intercalated ; wintercalated ; 0];
        indic = [indic ; windic ; 1];
        Y = [Y ; wY ; 0 ];
        dY = [dY ; wdY ; 10];
        minY = [minY ; wminY ; -10];
        dY_init = [dY_init ; wdY ; 10];
        compactcorr = [compactcorr ; wcompactcorr;  0];
%        istg = [istg ; 0 * wY ; 0];
        lat = [lat ; wlat ; mean(wlat) ];
        long = [long ; wlong ; mean(wlong) ];
        overburden=[overburden ; woverburden ; 0];
        istg=[istg ; wcont ; 0];
        sitecoords(curreg,:) = [wlat(1) wlong(1)];
        sitenames{curreg} = ['Holo-' datPX.textdata{sub(1)+1,1}];
        siteid(curreg,1)=wdatid(1);
        sitelen(curreg,1) = length(sub) + 1;
        
%         datid = [datid ; wdatid ; (idHolo+HoloRegions(curreg)*1e3); (idHolo+HoloRegions(curreg)*1e3)];
%         time1 = [time1 ; wtime1 ; 2016 ; 1950-13000-1000];
%         time2 = [time2 ; wtime2 ; 2017 ; 1950-13000+1000];
%         limiting = [limiting ; wlimiting ; 0; 0];
%         intercalated = [intercalated ; wintercalated ; 0; 0];
%         indic = [indic ; windic ; 1; 1];
%         Y = [Y ; wY ; 0 ; -62000];
%         dY = [dY ; wdY ; 10; 5000];
%         minY = [minY ; wminY ; -10; -77000];
%         dY_init = [dY_init ; wdY ; 10; 10000];
%         compactcorr = [compactcorr ; wcompactcorr; 0; 0];
% %        istg = [istg ; 0 * wY ; 0];
%         lat = [lat ; wlat ; mean(wlat); mean(wlat) ];
%         long = [long ; wlong ; mean(wlong) ; mean(wlong) ];
%         overburden=[overburden ; woverburden ; 0; 0];
%         istg=[istg ; wcont ; 0; 0];
%         sitecoords(curreg,:) = [wlat(1) wlong(1)];
%         sitenames{curreg} = ['Holo-' datPX.textdata{sub(1)+1,1}];
%         siteid(curreg,1)=wdatid(1);
%         sitelen(curreg,1) = length(sub) + 2;
end

% convert uniform data to approximately normal
% assuming that dY for the normal data is one sigma, then the new dY for the
% approximation to the uniform:  4 * old-dY = (b-a)
% to convert the (b-a) to normal sigma: sigma^2 = variance = (1/12)*(b-a)^2
% So, the new sigma would be the square root of that, then it needs to be
% multiplied by 2, because we want dY = sigma  = 2*old dY * (1/12)^(1/2) * 2
% = 4/(sqrt(12)) * (old dY) = 2/sqrt(3) * old dY
sub_unif=find(limiting==0 & indic==2);
dY(sub_unif)= (2/(sqrt(3))) * dY_init(sub_unif);

meantime=mean([time1 time2],2);
Ycv = sparse(diag(dY.^2));
dY0=dY;
Ycv0 = Ycv;

sub_interc = find(limiting==0 & intercalated ==1);
Y_init=Y;
if isempty(sub_interc)
else
    Y(sub_interc) =  224.15*overburden(sub_interc) +Y(sub_interc);
    dY(sub_interc) = (Y(sub_interc)-minY(sub_interc))./2;
end

meanSL = 0;
Y = Y - meanSL;

compactcorr=sparse(compactcorr);

PX.datid=round(datid);
PX.time1=time1;
PX.time2=time2;
PX.limiting=limiting;
PX.intercalated=intercalated;
PX.indicator=indic;
PX.indic=indic;
PX.Y0=Y;
PX.Y=Y;
PX.dY = dY;
PX.compactcorr=compactcorr;
PX.istg = istg;
PX.lat=lat;
PX.long=long;
PX.Ycv=sparse(diag(dY.^2));
PX.siteid=round(siteid);
PX.sitenames=sitenames;
PX.meantime=(PX.time1+PX.time2)/2;
PX.sitecoords = sitecoords;
PX.sitelen = sitelen;
PX.overburden=overburden;
PX.dY_init=dY_init;
PX.cont=cont;

clear datasets;

datasets{1}=PX;
datasets{1}.label='PX';
ii=1;
