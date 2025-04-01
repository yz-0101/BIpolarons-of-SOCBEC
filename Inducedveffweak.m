
function Veff=Inducedveffweak(kr1,kx1,kr2,kx2,a)

global k0 w0 g11 g22 g12  

k0=1;
w0=2;
g11=0.5;
g22=g11;
g12=g11;



% e0=0.02*1i;
% e1=0.0001*1i;

% % E1=@(kr,kx) PSOCBEC(kr,kx);          
% % uk=@(kr,kx) PSOCBECu(kr,kx);
% vk=@(kr,kx) PSOCBECv(kr,kx);

E11=@(kx,ky,kz) PSOCBECxyz(kx,ky,kz);
uk1=@(kx,ky,kz) PSOCBECuxyz(kx,ky,kz);
vk1=@(kx,ky,kz) PSOCBECvxyz(kx,ky,kz);

% GG11=@(z,kx,kr) uk(kr,kx).^2./(z - E1(kr,kx))-vk(kr,kx).^2./(z + E1(kr,-kx));
% GG12=@(z,kx,kr) uk(kr,kx).*vk(kr,-kx)./(z-E1(kr,kx))-uk(kr,-kx).*vk(kr,kx)./(z+E1(kr,-kx));

G11=@(z,kx,ky,kz) uk1(kx,ky,kz).^2./(z - E11(kx,ky,kz))-vk1(kx,ky,kz).^2./(z + E11(-kx,-ky,-kz));
G12=@(z,kx,ky,kz) uk1(kx,ky,kz).*vk1(-kx,-ky,-kz)./(z - E11(kx,ky,kz))-uk1(-kx,-ky,-kz).*vk1(kx,ky,kz)./(z +E11(-kx,-ky,-kz));

% c=@(z,qx,qr,kx,kr) z-E1(qr,qx)-(kx-qx).^2-kr.^2-qr.^2;
% b=@(kr,qr) 2*kr.*qr;
% f=@(z,qx,qr,kx,kr) 2*1i*(-1).^floor((-2*angle(c(z,qx,qr,kx,kr)-b(kr,qr))+angle(b(kr,qr).^2-c(z,qx,qr,kx,kr).^2))/2/pi)...
%     *pi./sqrt(b(kr,qr).^2-c(z,qx,qr,kx,kr).^2);
% F=@(z,qx,qr,kx,kr) (uk(qr,qx).^2.*f(z,qx,qr,kx,kr)+pi./((qx.^2+qr.^2))).*qr;

% f1=@(z,qx,qr,kx,kr) 2*pi./(sqrt(c(z,qx,qr,kx,kr)+b(kr,qr)).*sqrt(c(z,qx,qr,kx,kr)-b(kr,qr)));
% 
% F=@(z,qx,qr,kx,kr) (uk(qr,qx).^2.*f1(z,qx,qr,kx,kr)+pi./((qx.^2+qr.^2))).*qr;

% 
% P=@(z,kr,kx) -1/(8*pi^3)*integral2(@(qr,qx)  F(z,qx,qr,kx,kr) ,0,Inf,-Inf,Inf);

% a=-0.2;



T=8*pi/a;

% veff11=@(kr,kx) T^2*(GG11(0,kx,kr)+GG11(0,-kx,kr)+GG12(0,kx,kr)+GG12(0,-kx,kr));


 veff=@(kr1,kx1,kr2,kx2,theta) T^2*(G11(0,kx1-kx2,kr1-kr2.*cos(theta),-kr2.*sin(theta))+G11(0,kx2-kx1,kr2.*cos(theta)-kr1,kr2.*sin(theta))...
    +G12(0,kx2-kx1,kr2.*cos(theta)-kr1,kr2.*sin(theta))+G12(0,kx1-kx2,kr1-kr2.*cos(theta),-kr2.*sin(theta)));

Veff= integral(@(theta) veff(kr1,kx1,kr2,kx2,theta),0,2*pi); 

end