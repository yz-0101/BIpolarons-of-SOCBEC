
function Veff=Inducedveff(kr1,kx1,kr2,kx2,w00,a)


global k0  w0 g11 g22 g12  

k0=1;
w0=2;
g11=0.5;
g22=g11;
g12=g11;

% a=-1.5; %interaction strength


% e0=0.001*1i;


E1=@(kr,kx) PSOCBEC(kr,kx);          
uk=@(kr,kx) PSOCBECu(kr,kx);
% vk=@(kr,kx) PSOCBECv(kr,kx);

E11=@(kx,ky,kz) PSOCBECxyz(kx,ky,kz);
uk1=@(kx,ky,kz) PSOCBECuxyz(kx,ky,kz);
vk1=@(kx,ky,kz) PSOCBECvxyz(kx,ky,kz);



% GG11=@(z,kr,kx) uk(kr,kx)^2./(z - E1(kr,kx))-vk(kr,kx)^2./(z + E1(kr,-kx));
% GG12=@(z,kr,kx) uk(kr,kx)*vk(kr,kx).*(1./(z - E1(kr,kx))-1./(z +E1(kr,-kx)));

G11=@(z,kx,ky,kz) uk1(kx,ky,kz).^2./(z - E11(kx,ky,kz))-vk1(kx,ky,kz).^2./(z + E11(-kx,-ky,-kz));
G12=@(z,kx,ky,kz) uk1(kx,ky,kz).*vk1(kx,ky,kz)./(z - E11(kx,ky,kz))-uk1(kx,ky,kz).*vk1(kx,ky,kz)./(z +E11(-kx,-ky,-kz));

c=@(z,qx,qr,kx,kr) z-E1(qr,qx)-(kx-qx).^2-kr.^2-qr.^2;
b=@(kr,qr) 2*kr.*qr;

% f=@(z,qx,qr,kx,kr) 2*1i*(-1).^floor((-2*angle(c(z,qx,qr,kx,kr)-b(kr,qr))+angle(b(kr,qr).^2-c(z,qx,qr,kx,kr).^2))/2/pi)...
%     *pi./sqrt(b(kr,qr).^2-c(z,qx,qr,kx,kr).^2);

f=@(z,qx,qr,kx,kr) 2*pi./(sqrt(c(z,qx,qr,kx,kr)+b(kr,qr)).*sqrt(c(z,qx,qr,kx,kr)-b(kr,qr)));

F=@(z,qx,qr,kx,kr) (uk(qr,qx).^2.*f(z,qx,qr,kx,kr)+pi./((qx.^2+qr.^2))).*qr;



P=@(z,kr,kx) -1/(8*pi^3)*integral2(@(qr,qx)  F(z,qx,qr,kx,kr) ,0,Inf,-Inf,Inf,'AbsTol',1e-9,'RelTol',1e-3);


Ga=@(z,kr,kx,a) 1./((a./(8*pi))+P(z,kr,kx));



veff=@(kr1,kx1,kr2,kx2,theta) Ga(w00,kr1,kx1,a)*Ga(w00,kr2,-kx2,a)*G11(0,kx1-kx2,kr1-kr2*cos(theta),-kr2*sin(theta))...
     +Ga(w00,kr1,-kx1,a)*Ga(w00,kr2,kx2,a)*G11(0,kx2-kx1,kr2*cos(theta)-kr1,kr2*sin(theta))...
    +Ga(w00,kr2,kx2,a)*Ga(w00,kr2,-kx2,a)*G12(0,kx2-kx1,kr2*cos(theta)-kr1,kr2*sin(theta))+...
    Ga(w00,kr1,kx1,a)*Ga(w00,kr1,-kx1,a)*G12(0,kx1-kx2,kr1-kr2*cos(theta),-kr2*sin(theta));


 
% veff=@(kr1,kx1,kr2,kx2,theta) Z0^2*(Ga(0,kr1,kx1,a).*Ga(0,kr2,kx2,a).*G11(0,kx1-kx2,kr1-kr2.*cos(theta),-kr2.*sin(theta))...
%      +Ga(0,kr1,kx1,a).*Ga(0,kr2,kx2,a).*G11(0,kx2-kx1,kr2.*cos(theta)-kr1,kr2.*sin(theta))...
%     +Ga(0,kr1,kx1,a).^2.*G12(0,kx1-kx2,kr1-kr2.*cos(theta),-kr2.*sin(theta))+...
%     Ga(0,kr2,kx2,a).^2.*G12(0,kx1-kx2,kr1-kr2.*cos(theta),-kr2.*sin(theta)));
% 
Veff= integral(@(theta) veff(kr1,kx1,kr2,kx2,theta),0,2*pi,'AbsTol',1e-3,'RelTol',1e-3);
 


end