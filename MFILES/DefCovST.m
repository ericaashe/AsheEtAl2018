% DefCovST defines the covariance functions and bounds on the
% hyperparameters for the analysis

refyear=2000;

%define covariance functions
kMat1 = @(dx,thetas) thetas(1).^2 .* (1).*exp(-dx/thetas(2));
kMat3 = @(dx,thetas) thetas(1).^2 .* (1 + sqrt(3)*dx/thetas(2)).*exp(-sqrt(3)*dx/thetas(2));
kMat5 = @(dx,thetas) thetas(1).^2 .* (1 + (sqrt(5)*dx/thetas(2)).*(1 + sqrt(5)*dx/thetas(2)/3)).*exp(-sqrt(5)*dx/thetas(2));
kDELTA = @(dx,thetas) thetas(1).^2 .* (dx==0);
kDP = @(years1,years2,thetas) thetas(1).^2 * bsxfun(@times,(years1-refyear)',(years2-refyear));
kDELTAG = @(ad,thetas)thetas(1).^2.*(abs(ad)<1e-4).*(ad<360);

kMat3d = @(years1,years2,dx,thetas)thetas(1).^2.*(-3/(thetas(2).^2)).*dx.*exp(-sqrt(3)*dx/thetas(2)).*(-1+2*bsxfun(@ge,years1',years2));
kDPd = @(years1,years2,thetas) thetas(1).^2 * repmat((years1-refyear)',length(years2),1);
kMat3dd = @(dx,thetas)thetas(1).^2.*(3/(thetas(2).^2)).*(1-sqrt(3)*dx/thetas(2)).*exp(-sqrt(3)*dx/thetas(2));
kDPdd = @(years1,years2,thetas) thetas(1).^2 * ones(length(years2),length(years1));
kMat5d = @(years1,years2,dx,thetas)thetas(1).^2.*(-5*dx.*exp(-sqrt(5)*dx/thetas(2))).*(thetas(2)+sqrt(5)*dx)/(3*thetas(2).^3);
kMat5dd = @(dx,thetas)thetas(1).^2.*-(5*exp(-(sqrt(5)*dx)/thetas(2)).*(thetas(2).^2+sqrt(5)*thetas(2).*dx-5*dx.^2))/(3*thetas(2).^4);

clear modelspec;

% % defines global, high-frequency and white noise covariance functions
 cvfunc.G = @(dt1t2,thetas) kMat3(dt1t2,thetas(1:2));
 cvfunc.R = @(t1,t2,dt1t2,ad,thetas) (kMat3(dt1t2,thetas(1:2))).*(kMat1(ad,[1 thetas(3)])).*(ad<360);
 cvfunc.M = @(dt1t2,ad,thetas) (kMat3(dt1t2,thetas(1:2))).*(kMat1(ad,[1 thetas(3)])).*(ad<360);
 cvfunc.W = @(dt1t2,ad,thetas) kDELTAG(ad,1).*kDELTA(dt1t2,thetas(1));

% % defines first derivatives
 dcvfunc.G = @(t1,t2,dt1t2,thetas) kMat3d(t1,t2,dt1t2,thetas(1:2));
 dcvfunc.R = @(t1,t2,dt1t2,ad,thetas) (kMat3d(t1,t2,dt1t2,thetas(1:2))).*(kMat1(ad,[1 thetas(3)])).*(ad<360);
 dcvfunc.M = @(t1,t2,dt1t2,ad,thetas) (kMat3d(t1,t2,dt1t2,thetas(1:2))).*(kMat1(ad,[1 thetas(3)])).*(ad<360);
 dcvfunc.W = 0;

% % defines sedcond derivatives
 ddcvfunc.G = @(dt1t2,thetas) kMat3dd(dt1t2,thetas(1:2));
 ddcvfunc.R = @(t1,t2,dt1t2,ad,thetas) (kMat3dd(dt1t2,thetas(1:2))).*(kMat1(ad,[1 thetas(3)])).*(ad<360);
 ddcvfunc.M = @(dt1t2,ad,thetas) (kMat3dd(dt1t2,thetas(1:2))).*(kMat1(ad,[1 thetas(3)])).*(ad<360);
 ddcvfunc.W = 0;
 
modelspec(1).cvfunc = @(t1,t2,dt1t2,thetas,ad,fp1fp2)  cvfunc.G(dt1t2,thetas(1:2))+...
    cvfunc.R(t1,t2,dt1t2,ad,thetas(3:5)) +...
    cvfunc.M(dt1t2,ad,thetas(6:8)) +...
    cvfunc.W(dt1t2,ad,thetas(9));

modelspec(1).dcvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) dcvfunc.G(t1,t2,dt1t2,thetas(1:2)) +...
    dcvfunc.R(t1,t2,dt1t2,ad,thetas(3:5)) + dcvfunc.M(t1,t2,dt1t2,ad,thetas(6:8));

modelspec(1).ddcvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) ddcvfunc.G(dt1t2,thetas(1:2))...
    +ddcvfunc.R(t1,t2,dt1t2,ad,thetas(3:5)) + ddcvfunc.M(dt1t2,ad,thetas(6:8));

% modelspec(1).cvfunc = @(t1,t2,dt1t2,thetas,ad,fp1fp2)  kMat3(dt1t2,thetas(1:2))+...
%      (kMat3(dt1t2,thetas(3:4))).*(kMat1(ad,[1 thetas(5)])).*(ad<360) +...
%      (kMat3(dt1t2,thetas(6:7))).*(kMat1(ad,[1 thetas(8)])).*(ad<360) + ...
%      kDELTAG(ad,1).*(kDELTA(dt1t2,thetas(9))+thetas(10).^2);

% modelspec(1).cvfunc = @(t1,t2,dt1t2,thetas,ad,fp1fp2) kMat3(dt1t2,thetas(1:2))+...
%         (kMat3(dt1t2,thetas(3:4))).*(kMat3(ad,[1 thetas(5)])).*(ad<360) +...
%         (kMat3(dt1t2,thetas(6:7))).*(kMat3(ad,[1 thetas(8)])).*(ad<360) + ...
%         (kDELTAG(ad,1)).*(kDELTA(dt1t2,thetas(9))+...
%         thetas(10).^2);
% modelspec(1).cvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) cvfunc.G(dt1t2,thetas(1:2)) + cvfunc.H(dt1t2,thetas(3:4)) + cvfunc.W(dt1t2,ad,thetas(5));
% 

modelspec(1).traincv = @(t1,t2,dt1t2,thetas,errcv,ad,fp1fp2) modelspec(1).cvfunc(t1,t2,dt1t2,thetas,ad,fp1fp2) + errcv;

tluTGG = [
 
    2e4 2000 10e6 % Matern global amplitude
    7e3 10e3 1e5 % Temporal parameter
 
    1000 1000 2e4   % Matern Regional amplitude
    1500 1000 5e3   % Temporal parameters
    10 5 20         % Geographic length scale
    
    %0 0 0
    1000 500 1e4 % Matern HF regional amplitude
    500 300 10000 % temporal parameters
    2 .1 4.5 % geographic
 
    150 1e-2 1e4 % white noise
    0 1e-2 1e4 % local offset   
         ];

modelspec(1).thet0=tluTGG(:,1)';
modelspec(1).lb = tluTGG(:,2)';
modelspec(1).ub = tluTGG(:,3)';
modelspec(1).subfixed=[];

modelspec(1).sublength=[5 8]; % these are ones that will be tuned based only on tide gauges if use Optimize level 1
modelspec(1).subamp = [1 3 6 9 10];
modelspec(1).subamplinear = [];
modelspec(1).subampnonlinear = [];
modelspec(1).subampoffset = [];
modelspec(1).subampnoise = [9];
modelspec(1).subampHF = [];
modelspec(1).label='default';

