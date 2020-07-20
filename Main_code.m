
clc;
clear all;
close all;

index=170407;  %input the index number

%Take A,B and C values
A = mod(floor(index/100),10); 
B = mod(floor(index/10),10);
C = mod(index,10);

%Bandstop Filter specifications
Ap = 0.03+(0.01*A);   
Aa = 45+B;
Op1 = (C*100)+400;
Op2 = (C*100)+950;
Oa1 = (C*100)+500;
Oa2 = (C*100)+800;
Os = 2*((C*100)+1300);
T = (2*pi)/Os;

%Kaiser window specifications
Bt = min((Oa1-Op1), (Op2 - Oa2));   
Oc1 = Op1 + (Bt/2);
Oc2 = Op2 - (Bt/2);
dp = (10^(0.05*Ap) - 1)/(10^(0.05*Ap) + 1);
da = 10^(-0.05*Ap);
delta = min(dp, da);

%Stopband attenuation
A_a = (-20)*log10(delta);

if A_a<=21                                             
    alpha = 0;
elseif A_a>21 && Aa<= 50
    alpha = 0.5842*(A_a-21)^0.4 + 0.07886*(Aa-21);
else
    alpha = 0.1102*(A_a-8.7);
end

 
%A parameter D is chosen inorder to obtain N
if A_a<=21
    D = 0.9222;
else
    D = (A_a-7.95)/14.36;
end

%N is chosen such that it is the smallest odd integer value satisfying the
%inequality N ? ?sD Bt +1

N = ceil((Os*D)/Bt + 1);

if mod(N, 2) == 0
    N = N + 1;
end

N

% Length of the filter
n = -(N-1)/2:1:(N-1)/2;   

% Calculating beta
beta = alpha*sqrt(1-(2*n/(N-1)).^2);

% generating Io function
Io_beta = 1;
for k = 1:100
    Io_beta = Io_beta + ((1/factorial(k))*(beta/2).^k).^2;
end

Io_alpha = 1;
for k = 1:100
    Io_alpha = Io_alpha + ((1/factorial(k))*(alpha/2).^k).^2;
end

%Generating kaiser window transfer function
WnT=Io_beta/Io_alpha;

figure
stem(n, WnT);

title('Kaiser Window (Time Domain)');
ylabel('Amplitude');
xlabel('n');

%calculating hnT
n1 = -(N-1)/2:1:-1;
n2 = 1:1:(N-1)/2;
h1nT = 1./(n1*pi).*(sin(Oc1*n1*T)-sin(Oc2*n1*T));
h2nT = 1./(n2*pi).*(sin(Oc1*n2*T)-sin(Oc2*n2*T));
h0nT = 1+ 2/Os*(Oc1-Oc2);
hnT = [h1nT,h0nT,h2nT];
n_causal = 0:1:N-1;
HnT = hnT.*WnT;

figure
stem(n_causal, HnT);

title('Impulse Response - Time Domain');
ylabel('Amplitude');
xlabel('n');

%plotting the magnitude responce of the filter
[H,f] = freqz(HnT);
Hn = 20*log10(abs(H));
w = f*(Os/(2*pi));
figure
plot(w,Hn);

title('Magnitude Response of the filter');
xlabel('Frequency (rad/s)');
ylabel('Amplitude (dB)');

%plotting magnitude response of lower passband
ml = round(length(w)*(2*Oc1/Os));
wl = w(1:ml);
Hl = Hn(1:ml);
figure
plot(wl,Hl);

title('Magnitude Response of the lower passband')
xlabel('Frequency (rad/s)');
ylabel('Amplitude (dB)');


%plotting magnitude response of upper passband
mu = round(length(w)*(2*Oc2/Os));
wu = w(mu:length(w));
Hu = Hn(mu:length(w));
figure
plot(wu,Hu);

title('Magnitude Response of the upper passband')
xlabel('Frequency (rad/s)');
ylabel('Amplitude (dB)');

%generate the input signal
nx = 0:1:300;
wx1 = Oc1/2; %the middle frequency of the lower passband
wx2 = (Oc1+Oc2)/2; %the middle frequency of the stopband
wx3 = (Oc2+Os/2)/2; % the middle frequency of the upper passband

xn = cos(wx1*T*nx)+cos(wx2*T*nx)+cos(wx3*T*nx);

%plotting the time domain representation of xn and tome domain 
%representation offilterd signal
s = length(fft(xn)) + length(fft(HnT)) - 1;
Yn = fft(xn,s).*fft(HnT,s);

yn = ifft(Yn,s);
yn = yn(1:300);
figure
subplot(2,1,1)
stem(nx,xn);
axis tight;
title('Input signal / Time domain representation'); 
xlabel('n');
ylabel('Amplitude');

subplot(2,1,2)
stem(yn);
axis tight;
title('Filtered Output signal/Time domain representation'); 
xlabel('n');
ylabel('Amplitude');

%plotting the frequency domain representation of Xn and 
%frquency representation of corresponding output signal.
Xn = abs(fft(xn));
wx = linspace(-pi,pi, numel(Xn));
figure

subplot(2,1,1)
plot(wx,Xn);

title('Frequency domain representation of Xn');
xlabel('Frequency (rad/s)');
ylabel('Amplitude (dB)');

subplot(2,1,2)

Yn = abs(fft(yn));
wy = linspace(-pi,pi, numel(Yn));

plot(wy,Yn);

title('Frequency domain representation of filtered Output');
xlabel('Frequency (rad/s)');
ylabel('Amplitude (dB)');

%comparing expected output signal from bandpass filter theoritically (only 
%wx1 and wx3 are in the passband wx2 is in the stopband) and actual output

y_n = cos(wx1*T*nx) + cos(wx3*T*nx);
figure
subplot(2,1,1)
stem(y_n);
axis tight;
title('Expected output signal'); 
xlabel('n');
ylabel('Amplitude');

subplot(2,1,2)
stem(yn)
axis tight;
title('Actual output signal'); 
xlabel('n');
ylabel('Amplitude');













