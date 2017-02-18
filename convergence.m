trials=1;
StartSize=1;
Step=4;
utrack=zeros(trials,1);
htrack=zeros(trials,1);
exact=-31.4*20^3/(3*1E6/12);
for qwert=1:trials
   Lx=20; Ly=1; 
   Dy=floor((StartSize+Step*qwert)/2)*2;
   Dx=StartSize+qwert*Step;
   htrack(qwert)=max(Lx/Dx);
   Pilot
   utrack(qwert)=abs(exact-uMin);
end

%% Plot
R=[0.5*min(htrack),1.5*max(htrack),0.5*min(utrack),1.5*max(utrack)];
loglog(htrack,utrack,'-s')
set(gca, 'xdir','reverse')
axis(R)
xlabel('meshsize')
ylabel('Error')