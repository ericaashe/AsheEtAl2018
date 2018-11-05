%% runPlotMaps creates the maps in Ashe et al., 2018

    height_min = [-50 -50 -50 -50 -50 -50 -50 -50 -50 -50];
    height_max = [250 250 250 250 250 250 250 250 250 250 ];
    rate_min = [-35 -35 -35 -35 -35 -35 -35 -35 -35 -35 ]; %-5 -5 -5 -5 -5 -5 -5 -5 -5];
    rate_max = [10 10 10 10 10 10 10 10 10 10 ];
    % SD min and max for US Atlantic Coast
    SD_min = [0 0 0 0 0 0 0 0 0 0];
    SD_max = [30 30 30 30 30 30 30 30 30 30 ];
    minlat=23;
    maxlat=47;
    minlong=276;
    maxlong=295;
    [minlat maxlat minlong maxlong]
    lat_incr = .5;
    long_incr = .5;
    %     lat_incr = round((maxlat-minlat)/9)/10;
    %     long_incr = round((maxlong-minlong)/7)/10;
        Flat=minlat:lat_incr:maxlat;
        Flong=minlong:long_incr:maxlong;
runMapHeight_SD;
%runMapHeight_SD_171011_noGIA; 


%runMapHeight_SD_170523;
% %    runMapHeight_SD_170523
% end
