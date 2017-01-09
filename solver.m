%badpaca
%Solving a minimization problem with CG, SD, and HLSD

format compact
%% Universal Constants
luni = 2; nuni = (2^(luni)-1)^2; toluni = 10^-5;
disp(strcat('l=',num2str(luni),'_____','N=',num2str(sqrt(nuni)),...
    '_____','n=',num2str(nuni),'_____','tol=',num2str(toluni)));
%% Conjugate Gradient Method
n=nuni; N=sqrt(n); tol=toluni; x=ones(n,1); %size, tol & init
%INPUTS
h=1/(N+1); A=delsq(numgrid('S',N+2)); b=h^2.*ones(n,1); ct=1;
%INITIALIZATION
r=b-A*x; delt=transpose(r)*r; bdelt=transpose(b)*b; k=0; p=r; maxa=0; 
%CG
while delt > tol^2*bdelt
    s = A*p;
    alph = delt/(transpose(p)*s);
    x = x+alph*p;
    r = r-alph*s; delt0 = delt;
    delt = transpose(r)*r;
    p = r+delt/delt0*p;
    if alph > maxa
        maxa = alph;
    end
    ct = ct+1;
end
disp('Solution (x*)'); disp(transpose(x)); %the soln is identical 
                                            %no matter the method
cond1 = condest(A);
disp(strcat('Largest Step (CG)=',num2str(maxa),'_____',...
    '# Steps (CG)=',num2str(ct),'_____','Cond # (CG)=',num2str(cond1))); 

%% Steepest Descent Method
n2=nuni; N2=sqrt(n2); tol2=toluni; x2=ones(n2,1); %size, tol & init
%INPUTS
h2=1/(N2+1); A2=delsq(numgrid('S',N2+2)); b2=h2^2.*ones(n2,1); ct2=1;
%INITIALIZATION
r2=b2-A2*x2; alph2=transpose(r2)*r2/(transpose(r2)*A2*r2);
delt2=transpose(r2)*r2; bdelt2=transpose(b2)*b2; maxa2=0;
%SD
while delt2 > tol2^2*bdelt2
    x2 = x2+alph2*r2;
    r2 = b2-A2*x2;
    alph2 = transpose(r2)*r2/(transpose(r2)*A2*r2);
    delt2 = transpose(r2)*r2;
    if alph2 > maxa2
        maxa2 = alph2;
    end
    ct2 = ct2+1;
end
cond2 = condest(A2);
disp(strcat('Largest Step (SD)=',num2str(maxa2),'_____',...
    '# Steps (SD)=',num2str(ct2),'_____','Cond # (SD)=',num2str(cond2))); 

%% HLSD Method
n3=nuni; N3=sqrt(n3); tol3=toluni; x3=ones(n3,1); %size, tol & init
%INPUTS
h3=1/(N3+1); A3=delsq(numgrid('S',N3+2)); b3=h3^2.*ones(n3,1); ct3=1;
%INITIALIZATION
r3=b3-A3*x3; alph3=transpose(r3)*r3/(transpose(r3)*A3*r3); iter=1;
delt3=transpose(r3)*r3; bdelt3=transpose(b3)*b3; maxa3=0;
%HLSD
while delt3 > tol3^2*bdelt3
    x3 = x3+alph3*r3;
    r3 = b3-A3*x3;
    if mod(iter,2)==0 %only change alpha every other step
        alph3 = transpose(r3)*r3/(transpose(r3)*A3*r3);
    end
    delt3 = transpose(r3)*r3;
    iter = iter + 1;
    if alph3 > maxa3
        maxa3 = alph3;
    end
    ct3 = ct3+1;
end
cond3 = condest(A3);
disp(strcat('Largest Step (HLSD)=',num2str(maxa3),'_____',...
    '# Steps (HLSD)=',num2str(ct3),'_____',...
    'Cond # (HLSD)=',num2str(cond3))); 