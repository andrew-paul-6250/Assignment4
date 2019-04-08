%% Assignment 4 Part 2 - Andrew Paul 100996250
% The second section of this assingment involved doing a transient analysis
% on the given circuit.
% 
% After inspecting the circuit it was found that this is a low pass filter
% and that is a linear cirucit.
%
% The frequency response would be that low frequency signals would be
% passed (lower than the corner frequency) and higher frequency singals
% would be attenuated.

close all;
clear;

%Initialize variables and matricies

G = zeros(8,8); 
C = zeros(8,8); 

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
C(1,3)= cap;
C(2,3)= -cap;
C(4,7)= -L;

time_step = 0.001;

Vout = zeros(1000,1);
Vin = zeros(1000,1);
Vsolp = zeros(8,1);

A = C/(time_step)+G;
time = zeros(1000,1);

i = 1;

F = zeros(1,8);

for t=0:time_step:1
    if(t<0.03)
        F(1,7) = 0;
    else
        F(1,7) = 1;
    end
    
    time(i) = t;
    Vsol = inv(A)*(C*Vsolp/time_step + F');
    Vout(i) = Vsol(6);
    Vin(i) = Vsol(1);
    Vsolp = Vsol;
    i = i+1; 
end

figure(1)
subplot(2,1,2)
plot(time,Vout)
title('Step-Vout vs Time')
grid on

subplot(2,1,1)
plot(time,Vin)
title('Step-Vin vs Time')
grid on

freq = 1000;
freqOut = fft(Vout);
x = length(Vout);
y = fftshift(freqOut);
freqShift = (-x/2:x/2-1)*(freq/x); 
shift = abs(y).^2/x;

figure(2)
semilogy(freqShift,shift)
title('Step frequency spectrum output')
grid on

%Vin

freqIn = fft(Vin);
x = length(Vin);
y = fftshift(freqIn);
freqShift = (-x/2:x/2-1)*(freq/x); 
shift = abs(y).^2/x;

figure(3)
semilogy(freqShift,shift)
title('Step frequency spectrum input')
grid on

%Sine wave input

Vt = @(t) sin(2*pi*(1/0.03)*t);
Vout = zeros(1000,1);
Vin = zeros(1000,1);
Vsolp = zeros(8,1);

A = C/(time_step)+G;
time = zeros(1000,1);

i = 1;

F = zeros(1,8);

for t=0:time_step:1
    
    F(1,7) = Vt(t);
    
    time(i) = t;
    Vsol = inv(A)*(C*Vsolp/time_step + F');
    Vout(i) = Vsol(6);
    Vin(i) = Vsol(1);
    Vsolp = Vsol;
    i = i+1; 
end

figure(4)
subplot(2,1,2)
plot(time,Vout)
title('Sine-Vout vs Time')
grid on

subplot(2,1,1)
plot(time,Vin)
title('Sine-Vin vs Time')
grid on

freq = 1000;
freqOut = fft(Vout);
x = length(Vout);
y = fftshift(freqOut);
freqShift = (-x/2:x/2-1)*(freq/x); 
shift = abs(y).^2/x;    

figure(5)
semilogy(freqShift,shift)
title('Sine frequency spectrum output')
grid on

freqIn = fft(Vin);
x = length(Vin);
y = fftshift(freqIn);
freqShift = (-x/2:x/2-1)*(freq/x); 
shift = abs(y).^2/x;  

figure(6)
semilogy(freqShift,shift)
title('Sine frequency spectrum input')
grid on

% Gaussian input

Vt = @(t) exp(-(1/2)*((t-0.06)/(0.03))^2);
Vout = zeros(1000,1);
Vin = zeros(1000,1);
Vsolp = zeros(8,1);

A = C/(time_step)+G;
time = zeros(1000,1);

i = 1;

F = zeros(1,8);

for t=0:time_step:1
    
    F(1,7) = Vt(t);
    
    time(i) = t;
    Vsol = inv(A)*(C*Vsolp/time_step + F');
    Vout(i) = Vsol(6);
    Vin(i) = Vsol(1);
    Vsolp = Vsol;
    i = i+1; 
end

figure(7)
subplot(2,1,2)
plot(time,Vout)
title('Gaussian-Vout vs Time')
grid on

subplot(2,1,1)
plot(time,Vin)
title('Gaussian-Vin vs Time')
grid on

freq = 1000;
freqOut = fft(Vout);
x = length(Vout);
y = fftshift(freqOut);
freqShift = (-x/2:x/2-1)*(freq/x); 
shift = abs(y).^2/x;     

figure(8)
semilogy(freqShift,shift)
title('Gaussian frequency spectrum output')
grid on

freqIn = fft(Vin);
x = length(Vin);
y = fftshift(freqIn);
freqShift = (-x/2:x/2-1)*(freq/x);
shift = abs(y).^2/x;

figure(9)
semilogy(freqShift,shift)
title('Gaussian frequency spectrum input')
grid on

%Larger time step for Sine wave input

time_step = 0.1;

Vt = @(t) sin(2*pi*(1/0.03)*t);
Vout = zeros(1000,1);
Vin = zeros(1000,1);
Vsolp = zeros(8,1);

A = C/(time_step)+G;
time = zeros(1000,1);

i = 1;

F = zeros(1,8);

for t=0:time_step:1
    
    F(1,7) = Vt(t);
    
    time(i) = t;
    Vsol = inv(A)*(C*Vsolp/time_step + F');
    Vout(i) = Vsol(6);
    Vin(i) = Vsol(1);
    Vsolp = Vsol;
    i = i+1; 
end

figure(10)
subplot(2,1,2)
plot(time,Vout)
title('Sine-Vout vs Time: Larger Time Step')
grid on

subplot(2,1,1)
plot(time,Vin)
title('Sine-Vin vs Time: Larger Time Step')
grid on

freq = 1000;
freqOut = fft(Vout);
x = length(Vout);
y = fftshift(freqOut);
freqShift = (-x/2:x/2-1)*(freq/x); 
shift = abs(y).^2/x;    

figure(11)
semilogy(freqShift,shift)
title('Sine frequency spectrum output: Larger Time Step')
grid on

freqIn = fft(Vin);
x = length(Vin);
y = fftshift(freqIn);
freqShift = (-x/2:x/2-1)*(freq/x); 
shift = abs(y).^2/x;    

figure(12)
semilogy(freqShift,shift)
title('Sine frequency spectrum input: Larger Time Step')
grid on

%%
% From the time domain analyiss of the circuit it was found that as the time step was increased 
% the output waveform would appear to have more distorsion, thus a smaller time
% step would give a cleaner signal and should be used for analysis.

