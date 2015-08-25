%%Clear workspace
clear
clc
close all
%%Define the anon functions

%Gaussian burst
u=@(a,taui,tauo,sig,t) a*exp(-((t-taui-tauo).^2)./(2*sig.^2));

%Function
bates=@(t,t1,t2) (exp(-t./t1)-exp(-t./t2)).*(t>=0);
N=@(t,t0,sig) (t>=t0).*exp(-(t-t0).^2/2*(sig.^2))./sqrt(2*pi*sig);
G=@(t,t0,h,l) (((l.^h).*((t-t0).^(h-1)).*exp(-l.*(t-t0)))/gamma(h)).*(t>=t0);

%%Parameters: [stim amp taui tauo sig]

b=[0 0.4 0 1 1.2; ...
    4 0.78 0 0 0.3; ...
    15 0.33 2.7 0 1.5; ...
    19 0.62 0 0 0.3];
sfa=[.34 .77 .23];
scl=.05;
h=0.9445;
l=0.1037;
t0=2.9971;

t=0:.1:30;
scrni=zeros(size(t));
sf=scrni;
SCR=sf;
SF=sf;

%Create RF
rf=G(t,t0,h,l);
sfcheck=diff(b(:,1));
k=1;

for n=1:size(b,1)
    scrni=scrni+u(b(n,2),b(n,3),b(n,1)+b(n,4),b(n,5),t);
    if n<size(b,1) & sfcheck(n)>=5
        time=b(n,1)+5;
        sf=sf+u(sfa(k),time,0,0.3,t);
        k=k+1;
        while time+2<b(n+1,1)
            time=time+2;
            sf=sf+u(sfa(k),time,0,0.3,t);
            k=k+1;
        end
    end
end

SCR=conv(scrni,rf);
SF=conv(sf,rf);

figure(1)
subplot(3,2,1)
plot(t,scrni,'r','LineWidth',2)
title('SCR Neural Input')
ylim([0 1])
subplot(3,2,2)
plot(t,SCR(1:length(t)),'g','LineWidth',2)
title('SCR')
ylim([0 2])
subplot(3,2,3)
plot(t,sf,'r','LineWidth',2)
ylabel('Conductance (\muS)','FontSize',14)
title('SF Neural Input')
ylim([0 1])
subplot(3,2,4)
plot(t,SF(1:length(t)),'g','LineWidth',2)
title('SF')
ylim([0 2])
subplot(3,2,5)
plot(t,u(scl,mean(b(:,1)),0,0.3,t),'r','LineWidth',2)
title('SCL Neural Input')
xlabel('Time (s)')
ylim([0 1])
subplot(3,2,6)
SCL=conv(u(scl,mean(b(:,1)),0,0.3,t),ones(size(t)));
plot(t,SCL(1:length(t)),'g','LineWidth',2)
title('SCL')
xlabel('Time (s)')
ylim([0 2])

figure
SCfinal=SCR+SCL+SF;
plot(t,SCfinal(1:length(t)),'k','LineWidth',2)
xlabel('Time (s)')
ylabel('Conductance (\muS)')
title('SC')

figure
RF=conv(rf,N(t,0,0.5645))
plot(t,RF(1:length(t)),'b','LineWidth',2)
title('Response Function')
xlabel('Time (s)')

figure
plot(t,bates(t,6,3),'LineWidth',2)
title('Bateman Function')
xlabel('Time (s)')