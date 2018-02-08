%Taylor Bahm Fall 2017
%This function uses the 4th order Runge Kutta Method for system of 4 ODEs 
function myAns = RKsystemof4Plot(F1, F2, F3, F4, t0, A0, B0, C0, D0, tfinal, stepSize)
numSteps = (tfinal - t0)/stepSize + 1; %need to include one extra element in order to include both boundaries
myArray = zeros(numSteps,5); %five column array of zeros that has numStep rows
A_n = A0; %initlaize variables
B_n = B0;
C_n = C0;
D_n = D0;
t_n = t0;

for i = 1:numSteps
    myArray(i,1) = t_n; %saves t_n value to the ith row and 1st column
    myArray(i,2) = A_n; %saves A_n value to the ith row and 2nd column
    myArray(i,3) = B_n; %saves B_n value to the ith row and 3nd column
    myArray(i,4) = C_n; %saves C_n value to the ith row and 4th column
    myArray(i,5) = D_n; %saves D_n value to the ith row and 5th column
    
    A = A_n; %assigns A_n value to A so can evaluate K values and not lose the value of A_n
    B = B_n; %^same but for B_n
    C = C_n; %^same but for C_n
    D = D_n; %^same but for D_n
    t = t_n; %^same but for t_n
    K1A = stepSize * eval(F1); %evaluates the value of K1A
    K1B = stepSize * eval(F2);
    K1C = stepSize * eval(F3);
    K1D = stepSize * eval(F4);
   
    A = A_n + K1A/2;  %new value of A to find K2A
    B = B_n + K1B/2;  
    C = C_n + K1C/2;
    D = D_n + K1D/2;
    t = t_n + stepSize/2;  
    K2A = stepSize * eval(F1); %evaluates the value of K2A
    K2B = stepSize * eval(F2);
    K2C = stepSize * eval(F3);
    K2D = stepSize * eval(F4);
    
    
    A = A_n + K2A/2;  %new value of A to find K3A
    B = B_n + K2B/2;  
    C = C_n + K2C/2;
    D = D_n + K2D/2; %t stays the same so don't need to change
    K3A = stepSize * eval(F1); %evaluates the value of K3A
    K3B = stepSize * eval(F2);
    K3C = stepSize * eval(F3);
    K3D = stepSize * eval(F4);
    
    A = A_n + K3A;  %new value of A to find K4A
    B = B_n + K3B;  
    C = C_n + K3C;
    D = D_n + K3D; 
    t = t_n + stepSize;
    K4A = stepSize * eval(F1); %evaluates the value of K4A
    K4B = stepSize * eval(F2);
    K4C = stepSize * eval(F3);
    K4D = stepSize * eval(F4);
    
    A_n = A_n + (1/6)*(K1A+2*K2A+2*K3A+K4A); %calculates the next A_n value using Runge-Kutta equation
    B_n = B_n + (1/6)*(K1B+2*K2B+2*K3B+K4B);
    C_n = C_n + (1/6)*(K1C+2*K2C+2*K3C+K4C);
    D_n = D_n + (1/6)*(K1D+2*K2D+2*K3D+K4D);
    t_n = t_n + stepSize; %calculates the next t_n value based on the step size
end
plot(myArray(:,1),myArray(:,2),'k',myArray(:,1),myArray(:,3),'k--',myArray(:,1),myArray(:,4),'k:',myArray(:,1),myArray(:,5),'k-.'); %the colon tells matlab to use every value in the array
title('Species in a PFR')
xlabel('Reactor Volume (m^3)')
ylabel('Concentration (mol/m^3)')
legend('A','B','C','D')
myArray