%% TGDefineCovFuncs

refyear=2010;
dYears = @(years1,years2) abs(bsxfun(@minus,years1',years2));
dYears0 = @(years1,years2) (bsxfun(@minus,years1',years2));
angd = @(Lat0,Long0,lat,long) (180/pi)*(atan2( sqrt( (cosd(lat) .* sind(long-Long0)).^2 + (cosd(Lat0) .* sind(lat) - sind(Lat0) .* cosd(lat) .* cosd(long-Long0)).^2),(sind(Lat0) .* sind(lat) + cosd(Lat0) .* cosd(lat) .* cosd(long-Long0))));
dDist = @(x1,x2) angd(repmat(x1(:,1),1,size(x2,1)),repmat(x1(:,2),1,size(x2,1)),repmat(x2(:,1)',size(x1,1),1),repmat(x2(:,2)',size(x1,1),1))' + 1e6*(bsxfun(@plus,x1(:,1)',x2(:,1))>1000);

kMat1 = @(dx,thetas) thetas(1).^2 .* (1).*exp(-dx/thetas(2));
kMat3 = @(dx,thetas) thetas(1).^2 .* (1 + sqrt(3)*dx/thetas(2)).*exp(-sqrt(3)*dx/thetas(2));
kMat5 = @(dx,thetas) thetas(1).^2 .* (1 + (sqrt(5)*dx/thetas(2)).*(1 + sqrt(5)*dx/thetas(2)/3)).*exp(-sqrt(5)*dx/thetas(2));
kSE = @(dx,thetas) thetas(1).^2 * exp(-(dx.^2)/(2*thetas(2).^2));
kDELTA = @(dx,thetas) thetas(1).^2 .* (dx==0);
kDP = @(years1,years2,thetas) thetas(1).^2 * bsxfun(@times,(years1-refyear)',(years2-refyear));
kCONST = @(thetas) thetas(1).^2;
kSIN = @(dt,thetas) thetas(1).^2 * exp(-2*sin(pi*dt/(thetas(2))).^2/thetas(3).^2);
kMODSIN = @(dt,thetas) kSIN(dt,thetas(1:3)) .* kSE(dt,[1 thetas(4)*thetas(2)]);
kFREQ = @(dt,thetas) thetas(1).^2 * cos(2*pi*dt/thetas(2));
kRQ = @(dx,thetas) thetas(1).^2 * (1 + dx.^2/(2*thetas(2)*thetas(3))).^-thetas(3);
kMatG = @(dx,thetas) thetas(1).^2 .* 2.^(1-thetas(3))./gamma(thetas(3)) .* (sqrt(2*thetas(3))*(dx+eps)/thetas(2)).^thetas(3) .* besselk(thetas(3),sqrt(2*thetas(3))*(dx+eps)/thetas(2));
kDELTAG = @(ad,thetas)thetas(1).^2.*(abs(ad)<1e-4).*(ad<360);

kMat3d = @(years1,years2,dx,thetas)thetas(1).^2.*(-3/(thetas(2).^2)).*dx.*exp(-sqrt(3)*dx/thetas(2)).*(-1+2*bsxfun(@ge,years1',years2));
kDPd = @(years1,years2,thetas) thetas(1).^2 * repmat((years1-refyear)',length(years2),1);
kMat3dd = @(dx,thetas)thetas(1).^2.*(3/(thetas(2).^2)).*(1-sqrt(3)*dx/thetas(2)).*exp(-sqrt(3)*dx/thetas(2));
kDPdd = @(years1,years2,thetas) thetas(1).^2 * ones(length(years2),length(years1));
kMat5d = @(years1,years2,dx,thetas)thetas(1).^2.*(-5*dx.*exp(-sqrt(5)*dx/thetas(2))).*(thetas(2)+sqrt(5)*dx)/(3*thetas(2).^3);
kMat5dd = @(dx,thetas)thetas(1).^2.*-(5*exp(-(sqrt(5)*dx)/thetas(2)).*(thetas(2).^2+sqrt(5)*thetas(2).*dx-5*dx.^2))/(3*thetas(2).^4);

%%%

clear modelspec;

modelspec(1).cvfunc = @(t1,t2,dt1t2,thetas,ad,fp1fp2)  (kMat3(dt1t2,thetas(1:2))).*(kMat1(ad,[1 thetas(3)])).*(ad<360) +...
 kMat3(dt1t2,thetas(4:5)).*(kMat1(ad,[1 thetas(6)])).*(ad<360) + kDELTAG(ad,1).*(kDELTA(dt1t2,thetas(7))+thetas(8).^2);
modelspec(1).dcvfunc = @(t1,t2,dt1t2,thetas,ad,fp1fp2)(kMat3d(t1,t2,dt1t2,thetas(1:2))).*(kMat1(ad,[1,thetas(3)])).*(ad<360)...
    +kMat3d(t1,t2,dt1t2,thetas(4:5)).*(kMat1(ad,[1,thetas(6)])).*(ad<360);
modelspec(1).ddcvfunc = @(t1,t2,dt1t2,thetas,ad,fp1fp2)  (kMat3dd(dt1t2,thetas(1:2))).*(kMat1(ad,[1 thetas(3)])).*(ad<360)...
    + kMat3dd(dt1t2,thetas(4:5)).*(kMat1(ad,[1 thetas(6)])).*(ad<360);
modelspec(1).traincv = @(t1,t2,dt1t2,thetas,errcv,ad,fp1fp2) modelspec(1).cvfunc(t1,t2,dt1t2,thetas,ad,fp1fp2) + errcv;


tluTGG = [

    100 10 1e3 % Matern LF amplitude
    300 30 3e2 % LF temporal parameters
    50 3 200
 
    20 10 1e3 % Matern MF amplitude
    15 1 30 % temporal parameters
    3 3 3 % geographic length scale
    
%     50 1 1e4 % Matern HF regional amplitude
%     10 0.1 15  % temporal parameters
    
    100 1e-2 1e4 % white noise
    100 1e-2 1e4 % local offset   
         ];

modelspec(1).thet0=tluTGG(:,1)';
modelspec(1).lb = tluTGG(:,2)';
modelspec(1).ub = tluTGG(:,3)';
modelspec(1).subfixed=[6 ];

modelspec(1).sublength=[3 6];
modelspec(1).subamp = [1 4 7 8 ];
modelspec(1).subamplinear = [];
modelspec(1).subampnonlinear = [];
modelspec(1).subampoffset = [8];
modelspec(1).subampnoise = [7];
modelspec(1).subampHF = [];
modelspec(1).label='default';

