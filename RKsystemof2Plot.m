%4th order Runge Kutta Method for a system of 2 ODEs
%Problem 2 (35pts)
function myAns = RKsystemof2Plot(F1, F2, t0, A0, B0, tfinal, stepSize)
numSteps = (tfinal - t0)/stepSize + 1; %need to include one extra in order to include both boundaries
myArray = zeros(numSteps,3); %3 column array of zeros that is numSteps rows
A_n = A0; %initialize values
B_n = B0;
t_n = t0;

for i = 1:numSteps
    myArray(i,1) = t_n; %saves t_n to the ith row and 1st column
    myArray(i,2) = A_n;
    myArray(i,3) = B_n;
    
    t = t_n; %assings t_n to t so can evaluate K values and not lose the value of t_n
    A = A_n;
    B = B_n;
    K1A = stepSize * eval(F1); %evaluates the value of K1A
    K1B = stepSize * eval(F2);

    t = t_n + stepSize/2; 
    A = A_n + K1A/2;
    B = B_n + K1B/2;
    K2A = stepSize * eval(F1); 
    K2B = stepSize * eval(F2);
    
    A = A_n + K2A/2;
    B = B_n + K2B/2;
    K3A = stepSize * eval(F1); 
    K3B = stepSize * eval(F2);
    
    t = t_n + stepSize;
    A = A_n + K3A;
    B = B_n + K3B;
    K4A = stepSize * eval(F1); 
    K4B = stepSize * eval(F2);
    
    A_n = A_n + (1/6)*(K1A + 2*K2A + 2*K3A + K4A);
    B_n = B_n + (1/6)*(K1B + 2*K2B + 2*K3B + K4B);
    t_n = t_n + stepSize;
        
end
myAns = myArray;
plot(myArray(:,1),myArray(:,2),'k',myArray(:,1),myArray(:,3),'k--')

