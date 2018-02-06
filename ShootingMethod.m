%Shooting Method
%Problem 2 (35pts)
function myAns = ShootingMethod(At,Bt,t0,A0,tfinal,stepSize,G1,G2,DR)

FirstRun = RKsystemof2Plot(At, Bt, t0, A0, G1, tfinal, stepSize);
R1 = FirstRun(end,2);

SecondRun = RKsystemof2Plot(At, Bt, t0, A0, G2, tfinal, stepSize);
R2 = SecondRun(end,2);

Answer = G2 + (DR - R2)*(G1 - G2)/(R1 - R2);

ArrayPlot = RKsystemof2Plot(At, Bt, t0, A0, Answer, tfinal, stepSize);

plot(ArrayPlot(:,1),ArrayPlot(:,2),'k')
title('Temperature Profile')
xlabel('Distance, r (m)')
ylabel('Temperature (C)')
myAns = Answer;
