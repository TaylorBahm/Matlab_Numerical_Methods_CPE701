%Taylor Bahm Fall 2017
%This function solves the system of PDEs for time dependent heat transfer
%and flow in a cylindrical pipe

function myAns = TempProfileInPipe(tN,dt,radius,dr,dP,rho,mu,Length,dz,Tcold,Thot,k,Cp)
%note that dP is the pressure drop not the pressure gradient
gradP = -dP/Length; %b/c dP should be -1000Pa 
numRows = round(tN/dt); %using the round function so that get an interger value
numCols = round(2*radius/dr);
v = zeros(numRows,numCols); %array of velocities
%need to put in boundary conditions for edges of pipe wall
%redundant b/c have no slip and alreay initialized w/ 0 but need if wall was moving
v(:,1) = 0; %BC for top of pipe - first column
v(:,numCols) = 0;  %BC for bottom of pip - last column
radiusvalues = linspace(-radius,radius,numCols); %array of radius values 

for t = 1:numRows - 1 %for loop for time steps
    for r = 2:numCols - 1
        v(t+1,r) = (mu*dt/rho)*((v(t,r+1)-2*v(t,r)+v(t,r-1))/dr^2+(1/radiusvalues(r))*(v(t,r+1)-v(t,r-1))/(2*dr))+dt*(-gradP)/rho + v(t,r);
    end
end
%start of first figure
figure(1)
plot(v(1,:),radiusvalues);
hold on
numPlots = 20;
for j = 2:numPlots
    plot(v(numRows/numPlots*j,:),radiusvalues);
    title('Time Dependent Velocity Profile - Horizontal Pipe')
    xlabel('Velocity (m/s)')
    ylabel('Distance (m)')
end
numTimes = tN/dt; %number of times steps
numRs = round(2*radius/dr); %number of radius steps
numZs = round(Length/dz); %number of height steps
T = ones(numTimes,numZs,numRs)*Tcold; %Initalizes array with all TInital
T(:,1,:) = Thot; %boundary test
alpha = k/(rho*Cp);
%Triple nested for loop to go through every r, z, and t
for t = 1:numTimes-1
    for z = 2:numZs-1
        for r = 2:numRs-1
                part1 = alpha*dt*(1/radiusvalues(r))*((T(t,z,r+1)-T(t,z,r-1))/(2*dr));
                part2 = alpha*dt*((T(t,z,r+1)-2*T(t,z,r)+T(t,z,r-1))/(dr^2));
                part3 = alpha*dt*((T(t,z+1,r)-2*T(t,z,r)+T(t,z-1,r))/(dz^2));
                part4 = -dt*v(t,r)*((T(t,z+1,r)-T(t,z-1,r))/(2*dz)); %added velocity term
                T(t+1,z,r) = part1 + part2 + part3 + part4 + T(t,z,r);
        end
    end
end
%start of second figure
figure(2)
zData = linspace(0,Length,numZs); %Creating an array for z data
numPlots = 50;
TPlot = squeeze(T(1,:,:)); %converts from a 3D arrary to a 2D array
imagesc(radiusvalues,zData,TPlot);
pbaspect([0.1 1 10])
colorbar;
colormap hot
drawnow; %shows the plot now, doesn't wait until code is over
for k = 1:numPlots
    pause(0.1); %to get 10 frames/s
    TPlot = squeeze(T(numTimes/numPlots*k,:,:));
    imagesc(radiusvalues,zData,TPlot);
    pbaspect([0.1 1 10])
    colorbar;
    colormap hot
    drawnow;
    xlabel('Radius (m)')
    ylabel('Length (m)')
    title('Time Dependent Heat Transfer With Velocity')
end

%TempProfileInPipe(10,0.0005,0.02,0.0005,10,1000,0.01,1,0.01,298,368,1,4000)




