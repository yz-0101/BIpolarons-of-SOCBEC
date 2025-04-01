

function [omega,zk]=Dispersion(kr,kx,a,x00)


global k0  w0 g11 g22 g12  

k0=1;
w0=2;
g11=0.5;
g22=g11;
g12=g11;

% a=-0.1; %interaction strength



e0=0.001*1i;

% e1=0.0001*1i;

E1=@(kr,kx) PSOCBEC(kr,kx);          
uk=@(kr,kx) PSOCBECu(kr,kx);
% vk=@(kr,kx) PSOCBECv(kr,kx);
% 
% E11=@(kx,ky,kz) PSOCBECxyz(kx,ky,kz);
% uk1=@(kx,ky,kz) PSOCBECuxyz(kx,ky,kz);
% vk1=@(kx,ky,kz) PSOCBECvxyz(kx,ky,kz);



% GG11=@(z,kr,kx) uk(kr,kx)^2./(z - E1(kr,kx))-vk(kr,kx)^2./(z + E1(kr,-kx));
% GG12=@(z,kr,kx) uk(kr,kx)*vk(kr,kx).*(1./(z - E1(kr,kx))-1./(z +E1(kr,-kx)));

% G11=@(z,kx,ky,kz) uk1(kx,ky,kz).^2./(z - E11(kx,ky,kz))-vk1(kx,ky,kz).^2./(z + E11(-kx,-ky,-kz));
% G12=@(z,kx,ky,kz) uk1(kx,ky,kz).*vk1(kx,ky,kz)./(z - E11(kx,ky,kz))-uk1(kx,ky,kz).*vk1(kx,ky,kz)./(z +E11(-kx,-ky,-kz));

c=@(z,qx,qr,kx,kr) z-E1(qr,qx)-(kx-qx).^2-kr.^2-qr.^2;
b=@(kr,qr) 2*kr.*qr;
% f=@(z,qx,qr,kx,kr) 2*1i*(-1).^floor((-2*angle(c(z,qx,qr,kx,kr)-b(kr,qr))+angle(b(kr,qr).^2-c(z,qx,qr,kx,kr).^2))/2/pi)...
%     *pi./sqrt(b(kr,qr).^2-c(z,qx,qr,kx,kr).^2);
f=@(z,qx,qr,kx,kr) 2*pi./(sqrt(c(z,qx,qr,kx,kr)+b(kr,qr)).*sqrt(c(z,qx,qr,kx,kr)-b(kr,qr)));

F=@(z,qx,qr,kx,kr) (uk(qr,qx).^2.*f(z,qx,qr,kx,kr)+pi./((qx.^2+qr.^2))).*qr;



P=@(z,kr,kx) -1/(8*pi^3)*integral2(@(qr,qx)  F(z,qx,qr,kx,kr) ,0,Inf,-Inf,Inf,'AbsTol',1e-8,'RelTol',1e-3);


Ga=@(z,kr,kx,a) 1./((a./(8*pi))+P(z+e0,kr,kx));

Si=@(z,kr,kx,a) real(Ga(z,kr,kx,a));


ff=@(x,kr,kx,a) x-(kr.^2+kx.^2)-Si(x,kr,kx,a);

fun=@(x) ff(x,kr,kx,a) ;

% x00=-30;

omega=fsolve(fun, x00);
% kr*2+(kx+0.2)^2/2+x00

zk=1/(1-(Si(omega+0.05,kr,kx,a)-Si(omega-0.05,kr,kx,a))/0.1);



end
