%badpaca
%Rough data smoothing function using CGM and IRLS

N = 128; h = 1/N; t = zeros(N+1,1); b = zeros(N+1,1); %some parameters
beta = 10^-3; noise = 0.1; %some more parameters
for ii=1:N+1
    t(ii) = (ii-1)*h;
end
for jj=1:N+1
    if t(jj)<0.25, b(jj) = 1; end
    if t(jj)>=0.25 && t(jj)<0.5, b(jj) = 2; end
    if t(jj)>=0.5 && t(jj)<0.7, b(jj) = 2-100*(t(jj)-0.5)*(0.7-t(jj)); end
    if t(jj)>=0.7 && t(jj)<=1, b(jj) = 4; end
end
noisev = randn(size(b))*mean(abs(b))*noise;
data = b + noisev;

D = sqrt(h)/sqrt(2).*eye(N+1); endv = vertcat(zeros(N-1,1),1); %D mtx
W = 1/sqrt(h).*(horzcat(-eye(N),endv)+horzcat(zeros(N,1),eye(N))); %W mtx
lside = transpose(D)*D+beta*transpose(W)*W; rside = transpose(D)*D*data;

%Conjugate Gradient Method
u=ones(N+1,1); tol = 10^-5; %size, tol
r=rside-lside*u; delt=transpose(r)*r; 
bdelt=transpose(rside)*rside; p=r;
while delt > tol^2*bdelt
    s = lside*p;
    alph = delt/(transpose(p)*s);
    u = u+alph*p;
    r = r-alph*s; delt0 = delt;
    delt = transpose(r)*r;
    p = r+delt/delt0*p;
end
ufinal = u;

%IRLS
u2=u; tol2 = 10^-2; eps = 10^-6; gamma = 10^-2; %params
rside2 = transpose(D)*D*data; bdelt2=transpose(rside2)*rside2; 
%init stuff ==========
coeffvect = zeros(N,1);
for kk=1:N
    if kk==1, coeffvect(kk) = h/sqrt(eps*h^2); end
    if kk~=1
        coeffvect(kk) = h/sqrt((u2(kk)-u2(kk-1))^2+eps*h^2);
    end
end
Dhat = diag(coeffvect,0); 
lside2 = transpose(D)*D+gamma*transpose(W)*Dhat*W; 
ufinal2 = lside2\rside2;

%plotting stuff below =========================
plot(t,data,'r','LineWidth',2); hold on
plot(t,ufinal2,'b','LineWidth',2); title('Gamma=10e-2, Noise=0.1')
legend('Raw Data','LSQ Solved Data','Location','NorthWest');