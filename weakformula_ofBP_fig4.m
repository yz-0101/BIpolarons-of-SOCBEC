


global k0  w0 g11 g22 g12  

k0=1;
w0=2;
g11=0.5;
g22=g11;
g12=g11;

e0=0.02*1i;

E1=@(kr,kx) PSOCBEC(kr,kx);          
uk=@(kr,kx) PSOCBECu(kr,kx);
vk=@(kr,kx) PSOCBECv(kr,kx);

E11=@(kx,ky,kz) PSOCBECxyz(kx,ky,kz);
uk1=@(kx,ky,kz) PSOCBECuxyz(kx,ky,kz);
vk1=@(kx,ky,kz) PSOCBECvxyz(kx,ky,kz);

GG11=@(z,kx,kr) uk(kr,kx).^2./(z - E1(kr,kx))-vk(kr,kx).^2./(z + E1(kr,-kx));
GG12=@(z,kx,kr) uk(kr,kx).*vk(kr,-kx)./(z-E1(kr,kx))-uk(kr,-kx).*vk(kr,kx)./(z+E1(kr,-kx));

G11=@(z,kx,ky,kz) uk1(kx,ky,kz).^2./(z - E11(kx,ky,kz))-vk1(kx,ky,kz).^2./(z + E11(-kx,-ky,-kz));
G12=@(z,kx,ky,kz) uk1(kx,ky,kz).*vk1(-kx,-ky,-kz)./(z - E11(kx,ky,kz))-uk1(-kx,-ky,-kz).*vk1(kx,ky,kz)./(z +E11(-kx,-ky,-kz));

a=-0.1;

T=8*pi*a;

veff=@(kr,kx) T^2*(GG11(0,kx,kr)+GG11(0,-kx,kr)+2*GG12(0,kx,kr));

% veff1=@(kr,kx) -T^2*(uk(kr,kx)+vk(kr,-kx)).^2./E1(kr,kx);



 
fs = 10;            % Sampling frequency  in kx                  
t = 1/fs;             % Sampling period       
L0= 16;             % Length of signal
L=L0*fs;             % total points
kx= (-L/2:L/2)*t;  % Coordinate vector

fs1 = 10;            % Sampling frequency   in kr                 
t1 = 1/fs1;             % Sampling period       
L01=16;             % Length of signal
L1=L01*fs1;             % total points
ky=(-L1/2:L1/2)*t1; 

fs2= 10;            % Sampling frequency   in kr                 
t2= 1/fs2;             % Sampling period       
L02=16;             % Length of signal
L2=L02*fs2;             % total points
kz= (-L2/2:L2/2)*t2; 

[KX,KY,KZ]=meshgrid(kx,ky,kz);

Veff=zeros(L1+1,L+1,L2+1);

for i=1:L2+1
    for j=1:L1+1
        for k=1:L+1
            Veff(i,j,k)=veff(sqrt(KY(i,j,k)^2+KZ(i,j,k)^2),KX(i,j,k));
        end
    end
end
% figure
% mesh(KX(:,:,L2/2+1),KY(:,:,L2/2+1),Veff(:,:,L2/2+1));
%  xlabel('$kx/k_{r}$','interpreter','latex','FontSize',20);
% 
%  ylabel('$ky/k_{r}$','interpreter','latex','FontSize',20);
% colorbar

Veff11=ifftshift(ifftshift(ifftshift(Veff,1),2),3);

VR=fftn(Veff11)/fs/fs1/fs2;

x=2*pi/L0*(-L/2:L/2);
y=2*pi/L01*(-L1/2:L1/2);
z=2*pi/L02*(-L2/2:L2/2);
[X,Y,Z]=meshgrid(x,y,z);

Vr=fftshift(fftshift(fftshift(VR,1),2),3);

figure        % Veff  in r space
s=pcolor(X(:,:,L2/2+1),Y(:,:,L2/2+1),real(Vr(:,:,L2/2+1))/(8*pi^3));
s.EdgeColor='none';
 xlabel('$x$','interpreter','latex');
 ylabel('$y$','interpreter','latex');
   tit='$V_{\rm{eff}}(x,y,z=0)$';
title(tit,'interpreter','latex')
 set(gca,'FontSize',30,'xlim',[-5,5],'xtick',[-4,-2,0,2,4],'ylim',[-5,5],'ytick',[-4,-2,0,2,4])
   colorbar('FontSize',20)
hold on
plot(kx,0.2+zeros(1,length(kx)),'--','linewidth',1,'color','r')
%    
% figure         % Veff  in r space
% mesh(X(:,:,L2/2+1),Y(:,:,L2/2+1),real(Vr(:,:,L2/2+1))/(8*pi^3));
%  colorbar
%  xlabel('$xk_{r}$','interpreter','latex','FontSize',20);
%  ylabel('$yk_{r}$','interpreter','latex','FontSize',20);
%   tit='$V_{\rm{eff}}(\bf{R})/\it{E_{r}}$';
% title(tit,'interpreter','latex')
%  set(gca,'xlim',[-5,5],'FontSize',24,'ylim',[-5,5],'zlim',[-7.5,2])
  
 figure          % Veff  in r space  y=0 line 
 plot(X(L/2+1,:,L2/2+1),real(Vr(L/2+1,:,L2/2+1))/(8*pi^3),'linewidth',2)
 xlabel('$x$','interpreter','latex');
 set(gca,'xlim',[-5,5],'FontSize',30,'xtick',[-4,-2,0,2,4],'ylim',[-8,1],'ytick',[-8,-4,0])
  ylabel('$V_{\rm{eff}}(x,y=z=0)$','interpreter','latex','FontSize',28);
 
 figure          % Veff  in r space  x=0 line 
plot(Y(:,L1/2+1,L2/2+1),real(Vr(:,L1/2+1,L2/2+1))/(8*pi^3),'linewidth',2);
 xlabel('$yk_{r}$','interpreter','latex');
 ylabel('$V_{\rm{eff}}(y,x=z=0)/\it{E_{r}}$','interpreter','latex');
 set(gca,'xlim',[-5,5],'FontSize',25,'xtick',[-4,-2,0,2,4],'ylim',[-7.5,1])

kr=0:0.1:3;
kx=-5:0.1:5;
[KX,KR]=meshgrid(kx,kr);
VEFF=zeros(length(kr),length(kx));
 

for i=1:length(kr)
    for j=1:length(kx)
        VEFF(i,j)=veff(KR(i,j),KX(i,j));
    end
end

 
figure                 % Veff  in k space 
s=pcolor(KX,KR,VEFF);
s.EdgeColor='none';
 xlabel('$k_{x}$','interpreter','latex');
 ylabel('$k_{r}$','interpreter','latex');
 tit='$V_{\rm{eff}}(\bf{\tilde{k}})$';
title(tit,'interpreter','latex')
 set(gca,'xlim',[-5,5],'FontSize',30,'xtick',[-4,-2,0,2,4])
 colorbar('Ticks',[-50,-30,-10],'FontSize',20)
hold on
plot(kx,zeros(1,length(kx)),'--','linewidth',1,'color','r')
 
 figure
 plot(kx,VEFF(1,:),'-','linewidth',2)
  xlabel('$k_{x}$','interpreter','latex');
 set(gca,'xlim',[-5,5],'FontSize',30,'xtick',[-4,-2,0,2,4],'ytick',[-60,-30,0])
  ylabel('$V_{\rm{eff}}(\tilde{k_{x}},\tilde{k_{r}}=0)$','interpreter','latex','FontSize',28);

 
 

% plot(kx,zeros(1,length(kx)),'linewidth',1,'color','r')
% 
%  set(gca,'xlim',[-5,5],'FontSize',24,'xtick',[-4,-2,0,2,4])




