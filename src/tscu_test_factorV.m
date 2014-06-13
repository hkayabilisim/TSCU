clear all
close all
clc


addpath('lib/export_fig/');
hm = scfread('HMD01M0101.scf');
wt = scfread('WTD01M0101.scf');

n = 100;

hmpoint = 759;
wtpoint = 814;

hm = [hm.A(hmpoint-n:hmpoint+n) hm.C(hmpoint-n:hmpoint+n) hm.G(hmpoint-n:hmpoint+n) hm.T(hmpoint-n:hmpoint+n)];
wt = [wt.A(wtpoint-n:wtpoint+n) wt.C(wtpoint-n:wtpoint+n) wt.G(wtpoint-n:wtpoint+n) wt.T(wtpoint-n:wtpoint+n)];
wt=wt/2000;
hm=hm/2000;
df=wt-hm;
t = linspace(0,8,size(wt,1))';

p1baseA = [t wt(:,1)];
save('tscu_test_factorV_p1baseA.txt','p1baseA','-ascii');

figure
plot(t,wt);
figure
plot(t,hm);
figure
plot(t,df);

%%


all = [t wt hm df];
col = {'t','p1A','p1C','p1G','p1T','p2A','p2C','p2G','p2T','dfA','dfC','dfG','dfT'};
for i=2:length(col)
    fd=fopen(sprintf('factorV%s.dat',col{i}),'w');
    for j=1:size(all,1)
        fprintf(fd,'%-8.6f %-8.6f\n',all(j,1),all(j,i));
    end
    fclose(fd);
end

