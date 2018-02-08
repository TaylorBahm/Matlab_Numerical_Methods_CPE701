%Taylor Bahm Fall 2017
%This function solves the PDE for a cylindrical 2 dimensional time dependent heat transfer problem
%that models the baking of a pie where convection occurs on all sides

function myAns = ModelingABakingPie(TInitial,Tair,tN,dt,radius,dr,height,dz,k,rho,Cp,Upan,Upie)
numTimes = round(tN/dt); %number of times steps
numRs = round(2*radius/dr); %number of radius steps
numZs = round(height/dz); %number of height steps
T = ones(numTimes,numZs,numRs)*TInitial; %Initalizes array with all TInital
alpha = k/(rho*Cp);
RadiusValues = linspace(-radius,radius,numRs); %Creating an array for r data
%Triple nested for loop to go through every r, z, and t
for t = 1:numTimes - 1
    for z = 1:numZs
        for r = 1:numRs
            if (z==1 && r==1) %flux at top left corner
                TZM1 = (-2*dz*Upie/k)*(T(t,z,r)-Tair) + T(t,z+1,r); %value of T(t,0,r)
                TRM1 = (-2*dr*Upan/k)*(T(t,z,r)-Tair) + T(t,z,r+1); %value of T(t,z,0)
                part1 = alpha*dt*((T(t,z,r+1)-2*T(t,z,r)+TRM1)/(dr^2));
                part2 = alpha*dt*(1/RadiusValues(r))*((T(t,z,r+1)-TRM1)/(2*dr));
                part3 = alpha*dt*((T(t,z+1,r)-2*T(t,z,r)+TZM1)/(dz^2));
                T(t+1,z,r) = part1 + part2 + part3 + T(t,z,r);
            elseif (z==1 && r==numRs) %flux at top right corner
                TZM1 = (-2*dz*Upie/k)*(T(t,z,r)-Tair) + T(t,z+1,r); %value of T(t,0,r)
                TRP1 = (-2*dr*Upan/k)*(T(t,z,r)-Tair) + T(t,z,r-1); %value of T(t,z,numRs+1)
                part1 = alpha*dt*((TRP1-2*T(t,z,r)+T(t,z,r-1))/(dr^2));
                part2 = alpha*dt*(1/RadiusValues(r))*((TRP1-T(t,z,r-1))/(2*dr));
                part3 = alpha*dt*((T(t,z+1,r)-2*T(t,z,r)+TZM1)/(dz^2));
                T(t+1,z,r) = part1 + part2 + part3 + T(t,z,r);
            elseif (z==numZs && r==1) %flux at bottom left corner
                TZP1 = (-2*dz*Upan/k)*(T(t,z,r)-Tair) + T(t,z-1,r); %value of T(t,numZs+1,r)
                TRM1 = (-2*dr*Upan/k)*(T(t,z,r)-Tair) + T(t,z,r+1); %value of T(t,z,0)
                part1 = alpha*dt*((T(t,z,r+1)-2*T(t,z,r)+TRM1)/(dr^2));
                part2 = alpha*dt*(1/RadiusValues(r))*((T(t,z,r+1)-TRM1)/(2*dr));
                part3 = alpha*dt*((TZP1-2*T(t,z,r)+T(t,z-1,r))/(dz^2));
                T(t+1,z,r) = part1 + part2 + part3 + T(t,z,r);
            elseif (z==numZs && r==numRs) %flux at bottom right corner
                TZP1 = (-2*dz*Upan/k)*(T(t,z,r)-Tair) + T(t,z-1,r); %value of T(t,numZs+1,r)
                TRP1 = (-2*dr*Upan/k)*(T(t,z,r)-Tair) + T(t,z,r-1); %value of T(t,z,numRs+1)
                part1 = alpha*dt*((TRP1-2*T(t,z,r)+T(t,z,r-1))/(dr^2));
                part2 = alpha*dt*(1/RadiusValues(r))*((TRP1-T(t,z,r-1))/(2*dr));
                part3 = alpha*dt*((TZP1-2*T(t,z,r)+T(t,z-1,r))/(dz^2));
                T(t+1,z,r) = part1 + part2 + part3 + T(t,z,r);
            elseif z == 1 %flux at top edge
                TZM1 = (-2*dz*Upie/k)*(T(t,z,r)-Tair) + T(t,z+1,r); %value of T(t,0,r)
                part1 = alpha*dt*((T(t,z,r+1)-2*T(t,z,r)+T(t,z,r-1))/(dr^2));
                part2 = alpha*dt*(1/RadiusValues(r))*((T(t,z,r+1)-T(t,z,r-1))/(2*dr));
                part3 = alpha*dt*((T(t,z+1,r)-2*T(t,z,r)+TZM1)/(dz^2));
                T(t+1,z,r) = part1 + part2 + part3 + T(t,z,r);
            elseif z == numZs %flux at bottom edge
                TZP1 = (-2*dz*Upan/k)*(T(t,z,r)-Tair) + T(t,z-1,r); %value of T(t,numZs+1,r)
                part1 = alpha*dt*((T(t,z,r+1)-2*T(t,z,r)+T(t,z,r-1))/(dr^2));
                part2 = alpha*dt*(1/RadiusValues(r))*((T(t,z,r+1)-T(t,z,r-1))/(2*dr));
                part3 = alpha*dt*((TZP1-2*T(t,z,r)+T(t,z-1,r))/(dz^2));
                T(t+1,z,r) = part1 + part2 + part3 + T(t,z,r);
            elseif r == 1 %flux at left edge
                TRM1 = (-2*dr*Upan/k)*(T(t,z,r)-Tair) + T(t,z,r+1); %value of T(t,z,0)
                part1 = alpha*dt*((T(t,z,r+1)-2*T(t,z,r)+TRM1)/(dr^2));
                part2 = alpha*dt*(1/RadiusValues(r))*((T(t,z,r+1)-TRM1)/(2*dr));
                part3 = alpha*dt*((T(t,z+1,r)-2*T(t,z,r)+T(t,z-1,r))/(dz^2));
                T(t+1,z,r) = part1 + part2 + part3 + T(t,z,r);
            elseif r == numRs %flux at right edge
                TRP1 = (-2*dr*Upan/k)*(T(t,z,r)-Tair) + T(t,z,r-1); %value of T(t,z,numRs+1)
                part1 = alpha*dt*((TRP1-2*T(t,z,r)+T(t,z,r-1))/(dr^2));
                part2 = alpha*dt*(1/RadiusValues(r))*((TRP1-T(t,z,r-1))/(2*dr));
                part3 = alpha*dt*((T(t,z+1,r)-2*T(t,z,r)+T(t,z-1,r))/(dz^2));
                T(t+1,z,r) = part1 + part2 + part3 + T(t,z,r);
            else %continuity equation for all center points of the pie
                part1 = alpha*dt*((T(t,z,r+1)-2*T(t,z,r)+T(t,z,r-1))/(dr^2));
                part2 = alpha*dt*(1/RadiusValues(r))*((T(t,z,r+1)-T(t,z,r-1))/(2*dr));
                part3 = alpha*dt*((T(t,z+1,r)-2*T(t,z,r)+T(t,z-1,r))/(dz^2));
                T(t+1,z,r) = part1 + part2 + part3 + T(t,z,r);
            end
        end
    end
end
figure(1)
HeightValues = linspace(0,height,numZs); %Creating an array for z data
numPlots = 20;
TPlot = squeeze(T(1,:,:)); %converts from a 3D arrary to a 2D array
imagesc(RadiusValues,HeightValues,TPlot)
%pbaspect([3 1 1]) %use this if you don't want a square pie
colorbar
colormap hot
drawnow %shows the plot now, doesn't wait until code is over
for k = 1:numPlots
    pause(.01) %to get 10 frames/s
    TPlot = squeeze(T(numTimes/numPlots*k,:,:));
    imagesc(RadiusValues,HeightValues,TPlot);
    %pbaspect([3 1 1]) %use this if you don't want a square pie
    colorbar
    colormap hot
    drawnow
    title('Time Dependent 2D Temp Profile - Pie')
    xlabel('Radius (m)')
    ylabel('Height (m)')
end
%gives temp of very center and corner of the pie
timedata = linspace(0,tN,numTimes);
for i = 1:numTimes
    TempInCenterPie(i) = T(i,round(numZs/2),round(numRs/2)); %temp in center of pie at time i - coldest part
    TempAtCornerofPie(i) = T(i,numZs,numRs/2); %temp at the corner of a pie - hottest part
end
figure(2)
plot(timedata,TempInCenterPie,[0,tN],[347,347],timedata,TempAtCornerofPie,[0,tN],[430 430])
title('Temp in of Pie')
xlabel('Time')
ylabel('Temp')
legend('Center of Pie','Required Cook Temp','Bottom edge of Pie','Max Pie Temp')
%To ensure a convergent solution the alpha value must be kept small so
%slight changes in dr will cause the solution to diverge

%for alumnium pan
%ModelingABakingPie(293.15,450,6000,.25,0.1016,0.003,0.0381,0.0004,0.5,1140,3500,10,5)

%for glass pan
%ModelingABakingPie(293.15,450,3000,.25,0.1016,0.003,0.0381,0.0004,0.5,1140,3500,100,5)

          