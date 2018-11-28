%% plot_figures
% plots RSL curve, curve with data, and residuals plots

% plot RSL curves
clf;

    plotdat.x=testlocs.X(:,3);
    plotdat.y=f1s;
    plotdat.dy=sds*2;
    clf;    
    plot(plotdat.x,plotdat.y,'k','linew',2);
    hold on;
    plot(plotdat.x,plotdat.y-plotdat.dy,'k--','linew',1);
    plot(plotdat.x,plotdat.y+plotdat.dy,'k--','linew',1);
%     xlim(xl);
%     ylim([-4500 100]);
    xlabel('Age (Common Era)');
    ylabel('RSL (mm)')
    pdfwrite([date_field '_' label1 '_RSL_Predictions']);
    
    subD=find(Y);
        for uuu=subD(:)'
        if wdataset.limiting(uuu)==0
            plot([wdataset.time1(uuu) wdataset.time2(uuu)],wdataset.Y0(uuu)-2*wdataset.dY(uuu)*[1 1],'r'); hold on;
            plot([wdataset.time1(uuu) wdataset.time2(uuu)],wdataset.Y0(uuu)+2*wdataset.dY(uuu)*[1 1],'r');
            plot([wdataset.time1(uuu) wdataset.time1(uuu)],wdataset.Y0(uuu)+2*wdataset.dY(uuu)*[-1 1],'r');
            plot([wdataset.time2(uuu) wdataset.time2(uuu)],wdataset.Y0(uuu)+2*wdataset.dY(uuu)*[-1 1],'r');
        elseif wdataset.limiting(uuu)==1
            plot([.5 .5]*(wdataset.time1(uuu)+wdataset.time2(uuu)),wdataset.Y(uuu)+[-2 2]*wdataset.dY(uuu),'Color',[.05,.7,.5]);
            plot([wdataset.time2(uuu) wdataset.time2(uuu)],[1 1]*(wdataset.Y(uuu)-2*wdataset.dY(uuu)),'Color',[.05,.7,.5]);
        elseif wdataset.limiting(uuu)==-1         
            plot([.5 .5]*(wdataset.time1(uuu)+wdataset.time2(uuu)),wdataset.Y(uuu)+[-2 2]*wdataset.dY(uuu),'Color',[.1,.35,.7]);
            plot([wdataset.time2(uuu) wdataset.time2(uuu)],[1 1]*(wdataset.Y(uuu)+2*wdataset.dY(uuu)),'Color',[.1,.35,.7]);
        end
    end
    pdfwrite([date_field '_' label1 '_RSLs_with_data']);
    clf;
    
%% calculate rates averaged over 200 year periods
    difftimestep=1000;        
        Mdiff = bsxfun(@eq,testlocs.X(:,3),testlocs.X(:,3)')-bsxfun(@eq,testlocs.X(:,3),testlocs.X(:,3)'+ difftimestep);
        Mdiff = Mdiff .* bsxfun(@eq,testlocs.reg,testlocs.reg');
        sub=find(sum(Mdiff,2)==0);
        Mdiff=Mdiff(sub,:);
        difftimes=bsxfun(@rdivide,abs(Mdiff)*testlocs.X(:,3),sum(abs(Mdiff),2));
        diffreg=bsxfun(@rdivide,abs(Mdiff)*testlocs.reg,sum(abs(Mdiff),2));
        Mdiff=bsxfun(@rdivide,Mdiff,Mdiff*testlocs.X(:,3));
        
        clear df2s dV2s dsd2s;
        for n=1:size(f1s,2)
            df2s(:,n)=Mdiff*f1s(:,n);
            dV2s(:,:,n)=Mdiff*Vs(:,:,n)*Mdiff';
        end
            dsd2s(:,:)=sqrt(diag(dV2s(:,:,:)));

            k=1;
            offsetA=0;
            box on;

            hold on
            plot(difftimes,df2s,'Color',[.5, .5, .5]);
            plot(difftimes,df2s+dsd2s,'--','Color',[.5, .5, .5]);
            plot(difftimes,df2s-dsd2s,'--','Color',[.5, .5, .5]);
            plot(difftimes,df2s+2*dsd2s,':','Color',[.5, .5, .5]);
            plot(difftimes,df2s-2*dsd2s,':','Color',[.5, .5, .5]);

%             xlim([-550 2010]);
%             yl=[-2 6];
%             ylim(yl);
            ylabel('Rate (mm/yr = m/ka)','Color','k');
        title([wtestlocs.names{1} ' (' num2str(wtestlocs.sites(1,1)) ')']);
        xlim(xl);
        pdfwrite([date_field '_1000' label1 '_Rates']); 


    fid=fopen(['Rates' labl '.tsv'],'w');
    fprintf(fid,['Site\tDate\tRate\t2 SD\n']);
    for i=1:length(difftimes)
        fprintf(fid,[wdataset.sitenames{subq} '\t']);
        fprintf(fid,'%0.3f\t',difftimes(i));
        fprintf(fid,'%0.3f\t',df2s(i));
        fprintf(fid,'%0.3f\n',2*dsd2s(i));
    end
    fclose(fid);  
    
%% calculate residuals
[Xp,Xi]=sort(testlocp.X(:,3));
resids=Y(Xi)-fp(Xi)';

%% plot residuals
figure;
plot(fp(Xi),resids,'*'); 
title('Residuals vs Predicted'); 
xlabel('Predicted RSLs (mm)'); 
ylabel('Residuals (mm)');
pdfwrite([date_field '_' label1 '_Residuals_Plot']);

%% plot autocorrelation plot of residuals
autocorr(resids);
pdfwrite([date_field '_' label1 '_Autocorr_Residuals']);
