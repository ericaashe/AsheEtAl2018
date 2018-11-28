function [fslopeavgF,sdslopeavgF,fsF,sdsF,fslopeavgdiffF,sdslopeavgdiffF,diffplusF,difflessF,passderivs,invcv,fmeanF,fsdF,GIAmean] = RegressRateField_ea(PX,modelspec,thetL,noiseMasks,Flat,Flong,firstyears,lastyears,trainsub,GIA0,passderivs,invcv)

defval('firstyears',0);
defval('lastyears',1800);
defval('Flat',35:.5:52);
defval('Flong',232:.5:239);
defval('GIA0',[]);
defval('noiseMasks',ones(length(thetL)));
defval('trainsub',[]);
defval('passderivs',[]);
defval('invcv',[]);

clear testsitedefF;
[FLONG,FLAT]=meshgrid(Flong,Flat);
Fsite=[1:length(FLAT(:))]';
allsites = [Fsite FLAT(:) FLONG(:)];
lpt = length(unique(union(firstyears,lastyears)));
fsF=zeros(lpt*length(Fsite),1);
fsF0=zeros(lpt*length(Fsite),1);
sdsF=zeros(lpt*length(Fsite),1);
VsF=sparse(lpt*length(Fsite),lpt*length(Fsite));
testlocsites=[];
testlocreg=[];
testloct=[];
testlocGIA=[];
counter=1;
for qqq=1:1000:length(Fsite)
    dosub=qqq:min(qqq+999,length(Fsite));
    dosub2=counter:(counter-1+length(dosub)*lpt);
    counter=dosub2(end)+1;
    clear testsitedefF;
    testsitedefF.sites=allsites(dosub,:);
    testsitedef.GIA=[];
       for ii=1:size(testsitedefF.sites,1)
            testsitedefF.names2=[num2str(testsitedefF.sites(ii,1)) '-' num2str(testsitedefF.sites(ii,2)) '-' num2str(testsitedefF.sites(ii,3))];
            testsitedefF.names=testsitedefF.names2;
            if length(GIA0)>0
                 testsitedefF.GIA(ii,:)=interp3(GIA0.Lat,GIA0.Lon,GIA0.times,GIA0.GIA,testsitedefF.sites(ii,2),testsitedefF.sites(ii,3),[firstyears lastyears]);
            end
       end
 
    testtF = union(firstyears,lastyears);
    [fsF0(dosub2),sdsF(dosub2),VsF(dosub2,dosub2),testlocsF,~,passderivs,invcv]=RegressHoloceneDataSets(PX,testsitedefF,modelspec,thetL,trainsub,noiseMasks,testtF,1950,[],passderivs,invcv);   
    testlocsites=[testlocsites ; testlocsF.sites];
    testlocreg=[testlocreg ; testlocsF.reg];
    testloct=[testloct ; testlocsF.X(:,3)];
    testlocGIA=[testlocGIA ; testlocsF.GIAproj];
end

    fsF=fsF0;

 [fslopeavgF,sdslopeavgF,fslopeavgdiffF,sdslopeavgdiffF,diffplusF,difflessF,fmeanF,fsdF,GIAmean]=SLRateCompare_ea(fsF,...
     full(VsF),testlocsites,testlocreg,testloct,testlocGIA,firstyears,lastyears,sdsF);

