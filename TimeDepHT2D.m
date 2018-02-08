%Taylor Bahm
%This function solves a time dependent 2D heat transfer problem in a
%rectangular geometry with convection on one edge

function myAns = TimeDepHT2D(TCold,THot,tN,dt,length,h,K,rho,Cp) %Only 1 length because this is a square

numTimes = tN/dt; %total number of times
numPos = length/h; %total number of positions

TArray = ones(numTimes,numPos,numPos)*TCold; %Creating initial array of 25C because everything starts at 25C

%Boundary Condition
TArray(:,:,numPos) = THot; %Puts hot temperature at every time, every x, and always at the last y

X = K/rho/Cp;
Tau = X*dt/h^2;

%Triple nested for loop to go through every x, every y, every time
for t = 1:numTimes-1
    for x = 2:numPos-1 %Only want to do middle, not boundaries
        for y = 2:numPos-1 %only want to do middle, not boundaries
            TArray(t+1,x,y) = Tau*(TArray(t,x+1,y) +  TArray(t,x-1,y) + TArray(t,x,y+1) + TArray(t,x,y-1)) + (1 - 4*Tau)*TArray(t,x,y);
        end
    end
end

xData = linspace(0,length,numPos); %Creating array for x data
yData = xData; %Only because it is a square
numPlots = 50;
TPlot = squeeze(TArray(1,:,:)); %Gives 2D array that can be plotted
imagesc(xData,yData,TPlot);
colorbar
colormap hot
drawnow %show the plot now, dont wait until code is over
for k = 1:numPlots
    pause(0.1) %to get 10 frames/s
    TPlot = squeeze(TArray(numTimes/numPlots*k,:,:));
    imagesc(xData,yData,TPlot);
    colorbar
    colormap hot
    drawnow
end

%TimeDepHT2D(25,100,500,0.1,1,0.01,238,2700,900)

