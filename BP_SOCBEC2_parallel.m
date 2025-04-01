

global k0  w0 g11 g22 g12  

k0=1;
w0=2;
g11=0.5;
g22=g11;
g12=g11;

dkr=0.25;
dkx=0.25;
kxmax=5.5;
krmax=4;

kr=0:dkr:krmax;
kx=-kxmax:dkx:kxmax;
[KX,KR]=meshgrid(kx,kr);
L=length(kx)*length(kr);
% 
% 
% % A=[-5,-4.5,-4,-3.5,-3,-2.5,-2,-1.5,-1,-0.5,-0.3,-0.2,-0.15,-0.1];
% 
% % OME1=zeros(length(kr),length(kx),length(A));
% % ZZk1=zeros(length(kr),length(kx),length(A));
% % HH=zeros(L,L,length(A));
% % DD=zeros(L,length(A));
% % PSI=zeros(L,L,length(A));
% 
% for k=1:1
% 
% 
% ome1=zeros(length(kr),length(kx));
% Zk1=zeros(length(kr),length(kx));
% 
% H=zeros(L,L);
% 
% parfor i=1:L
% 
%     [a1,a2]=Dispersion(KR(i),KX(i),-0.75,aa(i)-1);
%     ome1(i)=a1;
%     Zk1(i)=a2;
% 
% end
% 
% 
% OME1(:,:,k)=ome1;
% ZZk1(:,:,k)=Zk1;


 w00=min(ome1(1,:));



parfor i=2:L
    v=zeros(1,L);
    for j=1:i-1

        v(j)=Inducedveff(KR(i),KX(i),KR(j),KX(j),w00,-0.75);

    end
    H(i,:)=v;

end
H1=H;
% HH(:,:,k)=H;



ome2=fliplr(ome1);
figure
s=pcolor(KX,KR,ome1);
s.EdgeColor='none';
 colorbar

Zk2=fliplr(Zk1);
figure
s=pcolor(KX,KR,Zk1);
s.EdgeColor='none';
 colorbar

 ZK=Zk1.*Zk2;

 
 
% end


% H=1/(8*pi^3)*dkr*dkx*(H+H.');
% 
% 
% for i=1:L
%     H(:,i)=H(:,i)*KR(i);
% end
% 
% for i=1:L
% 
%     H(i,:)=ZK(i)*H(i,:);
% 
% end
% 
% for i=1:L
%     H(i,i)=ome1(i)+ome2(i);
% end
% 
% 
% 
% [V,D]=eig(H);
% [d,k1]=sort(diag(D));
% V=V(:,k1);
% 
% 
% DD(:,k)=d;
% PSI(:,:,k)=V;
% 


 psi=zeros(length(kr),length(kx));
 for i=1:L
 
    psi(i)=V(i,1);
 end

 figure
s=pcolor(KX,KR,psi);
s.EdgeColor='none';
 xlabel('$kx/k_{r}$','interpreter','latex')
 ylabel('$ky/k_{r}$','interpreter','latex');
colorbar

 figure
mesh(KX,KR,psi);
 colorbar


% 
% k=[288,492];
%  L1=length(k);
% 
% for i=1:L1
%     v=zeros(1,L);
%     parfor j=1:k(i)-1
% 
%         v(j)=Inducedveff(KR(k(i)),KX(k(i)),KR(j),KX(j),w00,-4.5);
% 
%     end
%     H(k(i),:)=v;
% end
% 
% for i=1:L1
%     v=zeros(L,1);
%     parfor j=k(i)+1:L
% 
%         v(j)=Inducedveff(KR(j),KX(j),KR(k(i)),KX(k(i)),w00,-4.5);
% 
%     end
%     H(:,k(i))=v;
% end

% 
% ome1=OME1(:,:,10);
% Zk1=ZZk1(:,:,10);
% H=HH(:,:,10);


















