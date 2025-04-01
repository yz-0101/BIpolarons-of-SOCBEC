
function u=PSOCBECu(kr,kx)
%kr=10*(1/a0) ,  1/a0 is the unit in 13PRA
% hbar^2/(m*a^2) is the energy unit  Er is 100/2=50

global k0  w0 g11 g22 g12 


p0=k0*cos(asin(w0/4));
B0=sqrt((4*p0)^2+w0^2)/2;
nz0=2*p0/B0;
uk0=sqrt((1+nz0)/2);
vk0=sqrt((1-nz0)/2);
e0=p0^2+1-B0;


w=sqrt((4*(p0+kx)).^2+w0^2)./2;
nzk=2*(p0+kx)./w;
uk=sqrt((1+nzk)/2);
vk=sqrt((1-nzk)/2);
e=(p0+kx).^2+kr.^2+1-w;
mu=g11*vk0^4+g22*uk0^4+2*g12*uk0^2*vk0^2+e0;
SIN=2*(g11*vk0^2*vk.^2+g22*uk0^2*uk.^2)+2*g12*(uk0*vk0*uk.*vk)+g12*(uk0^2*vk.^2+vk0^2*uk.^2);

w1=sqrt((4*(p0-kx)).^2+w0^2)./2;
nzk1=2*(p0-kx)./w1;
uk1=sqrt((1+nzk1)/2);
vk1=sqrt((1-nzk1)/2);
SIN1=2*(g11*vk0^2*vk1.^2+g22*uk0^2*uk1.^2)+2*g12*(uk0*vk0*uk1.*vk1)+g12*(uk0^2*vk1.^2+vk0^2*uk1.^2);

SIA=(g11*vk0^2*vk.*vk1+g22*uk0^2*uk.*uk1)+g12*(vk0*uk0*(uk.*vk1+uk1.*vk));
e1=(p0-kx).^2+kr.^2+1-w1;


% L=length(e-mu+SIN);
% u=zeros(1,L);
% h11=e-mu+SIN;
% h12=SIA;
% h21=-SIA;
% h22=-e1+mu-SIN1;
% 
% for i=1:L
%     h=zeros(2);
%     h(1,1)=h11(i);
%     h(1,2)=h12(i);
%     h(2,1)=h21(i);
%     h(2,2)=h22(i);
% [V,D]=eig(h);
% [~,k1]=sort(diag(D),'descend');
% V=V(:,k1);
% u(i)=V(1,1);
% end

h11=e-mu+SIN;
h12=SIA;
h22=-e1+mu-SIN1;
E=(h11+h22+sqrt((h11-h22).^2-4*h12.^2))/2;

v=[(h22-E)./h12;ones(1,length(h22-E))];
norm=sqrt(((h22-E)./h12).^2-1);
v1=v./norm;
u=v1(1,:);

% v=[h22-E;h12];
% vv=normalize(v,"norm");
% u=vv(1,:);

% v=[h12;E-h11];
% vv=normalize(v,"norm");
% u=vv(1,:);

end

