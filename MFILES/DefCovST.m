% DefCovST defines the covariance functions and bounds on the
% hyperparameters for the analysis of the spatio-temporal dataset in Ashe et al., 2018

refyear=1950;

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

% % defines global, regional-linear, regional non-linear, local, and white noise covariance functions
 cvfunc.G = @(dt1t2,thetas) kMat3(dt1t2,thetas(1:2));
 cvfunc.RL = @(t1,t2,dt1t2,ad,thetas) (kDP(t1,t2,thetas(1))).*(kMat5(ad,[1 thetas(2)])).*(ad<360); %(kMat5(dt1t2,thetas(1:2))).*(kMat5(ad,[1 thetas(3)])).*(ad<360);
 cvfunc.RN = @(t1,t2,dt1t2,ad,thetas) (kMat5(dt1t2,thetas(1:2))).*(kMat3(ad,[1 thetas(3)])).*(ad<360); %(kMat5(dt1t2,thetas(1:2))).*(kMat5(ad,[1 thetas(3)])).*(ad<360);
 cvfunc.L = @(dt1t2,ad,thetas) (kMat3(dt1t2,thetas(1:2))).*(kMat3(ad,[1 thetas(3)])).*(ad<360);
 cvfunc.W = @(dt1t2,ad,thetas) kDELTAG(ad,1).*kDELTA(dt1t2,thetas(1));

% % defines first derivatives
 dcvfunc.G = @(t1,t2,dt1t2,thetas) kMat3d(t1,t2,dt1t2,thetas(1:2));
 dcvfunc.RL = @(t1,t2,dt1t2,ad,thetas) 0;
 dcvfunc.RN = @(t1,t2,dt1t2,ad,thetas)(kMat5d(t1,t2,dt1t2,thetas(1:2))).*(kMat3(ad,[1 thetas(3)])).*(ad<360);
 dcvfunc.L = @(t1,t2,dt1t2,ad,thetas) (kMat3d(t1,t2,dt1t2,thetas(1:2))).*(kMat3(ad,[1 thetas(3)])).*(ad<360);
 dcvfunc.W = 0;

% % defines second derivatives
 ddcvfunc.G = @(dt1t2,thetas) kMat3dd(dt1t2,thetas(1:2));
 ddcvfunc.RL = @(t1,t2,dt1t2,ad,thetas) 0;
 ddcvfunc.RN = @(dt1t2,ad,thetas) (kMat5dd(dt1t2,thetas(1:2))).*(kMat3(ad,[1 thetas(3)])).*(ad<360);
 ddcvfunc.L = @(dt1t2,ad,thetas) (kMat3dd(dt1t2,thetas(1:2))).*(kMat3(ad,[1 thetas(3)])).*(ad<360);
 ddcvfunc.W = 0;
 
modelspec(1).cvfunc = @(t1,t2,dt1t2,thetas,ad,fp1fp2)  cvfunc.G(dt1t2,thetas(1:2))+...
    cvfunc.RL(t1,t2,dt1t2,ad,thetas(3:4)) +...
    cvfunc.RN(t1,t2,dt1t2,ad,thetas(5:7)) +...
    cvfunc.L(dt1t2,ad,thetas(8:10)) +...
    cvfunc.W(dt1t2,ad,thetas(11));

modelspec(1).dcvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) dcvfunc.G(t1,t2,dt1t2,thetas(1:2)) +...
    dcvfunc.RL(t1,t2,dt1t2,ad,thetas(3:4)) + dcvfunc.RN(t1,t2,dt1t2,ad,thetas(5:7)) +...
    dcvfunc.L(t1,t2,dt1t2,ad,thetas(8:10));

modelspec(1).ddcvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) ddcvfunc.G(dt1t2,thetas(1:2))...
    + ddcvfunc.RL(t1,t2,dt1t2,ad,thetas(3:4)) + ddcvfunc.RN(dt1t2,ad,thetas(5:7))...
    + ddcvfunc.L(dt1t2,ad,thetas(8:10));

modelspec(1).traincv = @(t1,t2,dt1t2,thetas,errcv,ad,fp1fp2) modelspec(1).cvfunc(t1,t2,dt1t2,thetas,ad,fp1fp2) + errcv;

tluTGG = [
 
    2e5 4e4 10e6    % Matern global amplitude
    70e3 30e3 10e5  % Temporal parameter
 
    1 1e-3 10       % Linear amplitude
    12 12 12        % Geographic length scale

    2500 100 4e4    % Matern regional amplitude
    4500 1000 5e3   % Temporal parameter
    12 12 12        % Geographic length scale
    
    2500 500 4e4    % Matern local amplitude
    250 100 1e3     % Temporal parameters
    3 3 3           % Geographic length scale
 
    1500 100 1e4    % white noise
         ];

modelspec(1).thet0=tluTGG(:,1)';
modelspec(1).lb = tluTGG(:,2)';
modelspec(1).ub = tluTGG(:,3)';
modelspec(1).subfixed=[4 7 10];

modelspec(1).sublength=[5 8]; 
modelspec(1).subamp = [1 3 5 8 11];
modelspec(1).subamplinear = [];
modelspec(1).subampnonlinear = [];
modelspec(1).subampoffset = [];
modelspec(1).subampnoise = [11];
modelspec(1).subampHF = [];
modelspec(1).label='default';

