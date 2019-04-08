%% Assignmnet 4 Part 1: Andrew Paul 100996250
% The first section of this assignment invloved circuit analysis of
% different nodes and components to define a G and C matrix for modelling
% the voltage and gain parameters at different regions in the circuit
%% G matrix
% 
% <<Gmatrix1.PNG>>
% 
%% C matrix
% 
% <<Cmatrix1.PNG>>
% 

close all;
clear;

%Initialize variables and matricies
G = zeros(8,8); 
C = zeros(8,8); 
F = zeros(8,1); 
R1 =1;
R2 = 2;
R3 = 10;
R4 = 0.1;
R0 = 1000;
cap = 0.25;
L = 0.2;
alpha = 100;

G1 = 1/R1;
G2 = 1/R2;
G3 = 1/R3;
G4 = 1/R4;
G0 = 1/R0;

% The rows of the G matrix are classified as [V1 Iin V2 V3 V4 V5 IL I4]

G(1, 1) = -G1;
G(1, 2) =  G1;
G(2, 1) =  G1;
G(1, 3) =  G1;
G(2, 3) = -G1-G2;
G(3, 4) = -G3;
G(2, 7) = -1;
G(3, 7) = 1;
G(4, 3) = 1;
G(4, 4) = -1;
G(5, 6) = G4;
G(5, 7) = -alpha*G4;
G(5, 8) = 1;
G(6, 6) = -G4-G0;
G(6, 7) = alpha*G4;
G(7, 1) = 1;
G(8, 5) = 1; 
G(8, 7) = -alpha;

% Create C matrix
C(1,1)= -cap;
C(2,1)= cap;
C(1,3)= -cap;
C(2,3)= -cap;
C(4,7)= -L;

% Solve F matrix
Vin = linspace(-10,10,100);
V3 = zeros(length(Vin),1);
V0 = zeros(length(Vin),1);
for i = 1:length(Vin)
    F(7,1) = Vin(i);
    V = G\F;
    V3(i) = V(4);
    V0(i) = V(6);
end

figure(1)
plot(Vin,V3);
xlabel('Vin')
ylabel('V3')
title('V3 vs Vin sweep')

figure(2)
plot(Vin,V0);
xlabel('Vin')
ylabel('V0')
title('V0 vs Vin sweep')

w = 2*pi*linspace(0,80,100);
V0 = zeros(length(w),1);
gain = zeros(length(w),1);

for i = 1:length(w)
    s = 1i*w(i);
    M = inv((G +((s).*C)))*F; 
    V0(i) = abs(M(5));
    gain(i) = 20*log10(abs(V0(i))/abs(M(1)));
end

figure(3)
plot(w,V0);
xlabel('w (rads/sec)')
ylabel('V0')
title('AC plot for V0')

figure(4)
semilogx(w,gain);
xlabel('w (rads/sec)')
ylabel('V0/V1 (dB)')
title('Gain');

V0 = zeros(length(w),1);
gain = zeros(length(w),1);

for i = 1:length(gain)
    pert = randn()*0.05;
    C(1,1) = cap*pert;
    C(1,3) = -cap*pert;
    C(2,1) = -cap*pert;
    C(2,3) = cap*pert;

    s = 1i*2*pi*pi;
    M = inv((G +((s).*C)))*F; 
    V0(i) = abs(M(5));
    gain(i) = 20*log10((V0(i))/abs(M(1)));
end

figure(5);
hist(gain,100);
xlabel('Gain')
ylabel('Counts')
title('Hist C')
