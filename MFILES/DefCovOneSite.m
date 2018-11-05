refyear=2010;

%define covariance functions
kMat3 = @(dx,thetas) thetas(1).^2 .* (1 + sqrt(3)*dx/thetas(2)).*exp(-sqrt(3)*dx/thetas(2));
kDELTA = @(dx,thetas) thetas(1).^2 .* (dx==0);

kMat3d = @(years1,years2,dx,thetas)thetas(1).^2.*(-3/(thetas(2).^2)).*dx.*exp(-sqrt(3)*dx/thetas(2)).*(-1+2*bsxfun(@ge,years1',years2));
kMat3dd = @(dx,thetas)thetas(1).^2.*(3/(thetas(2).^2)).*(1-sqrt(3)*dx/thetas(2)).*exp(-sqrt(3)*dx/thetas(2));

% defines global, high-frequency and white noise covariance functions
cvfunc.G = @(dt1t2,thetas) kMat3(dt1t2,thetas(1:2));
cvfunc.H = @(dt1t2,thetas) kMat3(dt1t2,thetas(1:2));
cvfunc.W = @(dt1t2,ad,thetas) kDELTA(dt1t2,thetas(1));

% defines temporal derivative of covariance function
dcvfunc.G = @(t1,t2,dt1t2,thetas) kMat3d(t1,t2,dt1t2,thetas(1:2));
dcvfunc.H = @(t1,t2,dt1t2,thetas) kMat3d(t1,t2,dt1t2,thetas(1:2));
dcvfunc.W = @(dt1t2,ad,thetas) 0;

% second derivative
ddcvfunc.G = @(dt1t2,thetas) kMat3dd(dt1t2,thetas(1:2));
ddcvfunc.H = @(dt1t2,thetas) kMat3dd(dt1t2,thetas(1:2));
ddcvfunc.W = @(dt1t2,ad,thetas) 0;

clear modelspec;

%define model specifications
modelspec(1).label = 'Full';

%model is sum of functions we defined above
modelspec(1).cvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) cvfunc.G(dt1t2,thetas(1:2)) + cvfunc.H(dt1t2,thetas(3:4)) + cvfunc.W(dt1t2,ad,thetas(5));
modelspec(1).dcvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) dcvfunc.G(t1,t2,dt1t2,thetas(1:2)) +dcvfunc.H(t1,t2,dt1t2,thetas(3:4));% + dcvfunc.R(t1,t2,dt1t2,ad,thetas(3:5));
modelspec(1).ddcvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) ddcvfunc.G(dt1t2,thetas(1:2))+ddcvfunc.H(dt1t2,thetas(3:4));% + ddcvfunc.R(dt1t2,ad,thetas(3:5));

%adds defined error covariance to the covariance function
modelspec(1).traincv = @(t1,t2,dt1t2,thetas,errcv,ad,fp1fp2) modelspec(1).cvfunc(t1,t2,dt1t2,thetas,ad,fp1fp2) + errcv;

%%%%%
%define starting values and lower and upper bounds of the hyperparameters
%%%%%
tluTGG = [
60e3 1e3 10e6      % global amplitude (starting, lower, and uppper); SL in mm; how much values can go up or down over a particular time period  
12e3 2e3 2e6    % temporal scale in years, how much it can very over a particular time scale
1e3 10 20e3      % high frequency amplitude (starting, lower, and uppper); SL in mm; how much values can go up or down over a particular time period  
1e3 100 2e3    % HF temporal scale in years, how much it can very over a particular time scale
1e3 1 2e3     % white noise    
];

modelspec(1).thet0=tluTGG(:,1)';
modelspec(1).lb = tluTGG(:,2)';
modelspec(1).ub = tluTGG(:,3)';
modelspec(1).subfixed=[ ]; 
modelspec(1).sublength=[]; % these are ones that will be tuned based only on tide gauges if use Optimize level 1
modelspec(1).subamp = [1 3 5];
modelspec(1).subampnoise = [3 5];
modelspec(1).subampoffset = [];

