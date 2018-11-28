%prep tide gauge predictions to do maps
    minlat=23;
    maxlat=47;
    minlong=276;
    maxlong=295;
    [minlat maxlat minlong maxlong]
    lat_incr = .5;
    long_incr = .5;
        Flat=minlat:lat_incr:maxlat;
        Flong=minlong:long_incr:maxlong;
runMapHeight_Rate_TG;