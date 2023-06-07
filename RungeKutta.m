function SM=RungeKutta(XB,XW,DelX)

dim=size(XB,2);
C=randi([1 2])*(1-rand);
r1=rand(1,dim);
r2=rand(1,dim);

K1 = 0.5*(rand*XW-C.*XB);
K2 = 0.5*(rand*(XW+r2.*K1.*DelX/2)-(C*XB+r1.*K1.*DelX/2));
K3 = 0.5*(rand*(XW+r2.*K2.*DelX/2)-(C*XB+r1.*K2.*DelX/2));
K4 = 0.5*(rand*(XW+r2.*K3.*DelX)-(C*XB+r1.*K3.*DelX));

XRK = (K1+2.*K2+2.*K3+K4);
SM=1/6*XRK;
end

