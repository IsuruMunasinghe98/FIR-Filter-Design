clc;
close all;

A=4; 
B=0; 
C=9;

Ap_ =0.03+(0.01*A);  % Maximum passband ripple (Ap ) = 0.07 dB 
Aa_ = 45+B;		    % Minimum stopband attenuation (Aa_) = 45 dB 
wp1 =(C*100)+300;	% Lower passband edge (wp1)	= 1200 rad/s 
wp2 =(C*100)+700;	% Upper passband edge (wp2) = 1600 rad/s 
wa1 =(C*100)+150;	% Lower stopband edge (wa1) = 1050 rad/s 
wa2 =(C*100)+800;	% Upper stopband edge (wa2) = 1700 rad/s 
ws  =2*((C*100)+1200);%Sampling frequency    (ws) = 4200 rad/s

disp("Ap_ :" + num2str(Ap_));
disp("Aa_ :" + num2str(Aa_));
disp("wp1 :" + num2str(wp1));
disp("wp2 :" + num2str(wp2));
disp("wa1 :" + num2str(wa1));
disp("wa2 :" + num2str(wa2));
disp("ws  :" + num2str(ws));

%finding the two cut-off frequencies of the filter using transition bands 
tr1=wp1-wa1;
tr2=wa2-wp2;
Bt = min(tr1,tr2);
wc1 = wp1-(Bt/2);% lower cut-off frequency 
wc2 = wp2+(Bt/2);% upper cut-off frequency 
disp("tr1 :" + num2str(tr1));
disp("tr2 :" + num2str(tr2));

del_p=(10^(0.05*Ap_)-1)/(10^(0.05*Ap_)+1);
del_a=10^(-0.05*Aa_); 
del = min(del_p,del_a);

Aa = -20*log10(del); % actual minimum stopband attenuation
Ap = 20*log10((1+del)/(1-del)); % actual maximum passband ripple

fprintf('\n');
disp("Lower cut-off frequency :" + num2str(wc1)); 
disp("Upper cut-off frequency :" + num2str(wc2)); 
fprintf('\n');
disp("Delta p :" + num2str(del_p)); 
disp("Delta a :" + num2str(del_a)); 
disp("Delta :" + num2str(del));
disp("Actual minimum stopband attenuation(Aa) :" + num2str(Aa)); 
disp("Actual maximum passband ripple(Ap) :" + num2str(Ap));

alpha = get_Alpha(Aa); 
D=get_D(Aa); 
N=get_N(D,Bt,ws);

fprintf('\n'); 
disp("Alpha                  :"+ num2str(alpha));
disp("D                      :"+ num2str(D));
disp("Order of the Filter(N) :" + num2str(N));
fprintf('\n');
 
I0_alpha = get_Bessel(alpha);
disp("I0_alpha :" + num2str(I0_alpha));

wk_nT = zeros(1, ((N-1)/2)+1); 
for n = 0 : ((N-1)/2)
    beta = alpha*sqrt(1-((2*n)/(N-1))^2); 
    I0_beta = get_Bessel(beta); 
    wk_nT(n+1) = I0_beta/I0_alpha;
end

h_nT = zeros(1,((N-1)/2)+1); 
T = (2*pi)/ws;
for n = 0 : ((N-1)/2) 
    if n == 0
        h_nT(n+1) = 2*(wc2-wc1)/ws;
    else
        h_nT(n+1) =(sin(wc2*n*T)-sin(wc1*n*T))/(n*pi);
    end 
end

 

%Concatenate two arrys and get function values for negative n 
wk_nT = cat(2, fliplr(wk_nT), wk_nT(2:((N-1)/2)+1));
h_nT = cat(2, fliplr(h_nT), h_nT(2:((N-1)/2)+1));

hwnT = zeros(1, (2*(N-1)/2)+1); 
for i = 1:(2*(N-1)/2)+1
    hwnT(i) = h_nT(i)*wk_nT(i);
end
%----------------------------------------------------------%

%Plot Kaiser window
figure, stem((-((N-1)/2):(N-1)/2), wk_nT); 
grid on;
title('Impulse response of Kaiser window'); 
ylabel('Amplitude');
xlabel('n')

%Plot causal impulse response 
figure, stem(0:N-1,hwnT); 
grid on;
title('Causal Impulse response — Kaiser Window'); 
ylabel('Amplitude');
xlabel('n')

[HwnT, f] = freqz(hwnT); 
w = f*ws/(2*pi);


%Plot magnitude response in 0:ws/2 range 
figure;
plot(w, 20*log10(abs(HwnT)));
title('Magnitude Response in the range (0,w_s/2) — Kaiser Window'); xlabel('Angular Frequency (rad/s) ');
ylabel('Magnitude (dB)'); grid on;
axis([0, 2100, -140, 15]);

%Magnitude response of the passband 
figure; 
plot(w, 20*log10(abs(HwnT)));
title('Magnitude Response of the Passband — Kaiser Window'); 
xlim([wc1 wc2]);
grid on;
xlabel('Angular Frequency (rad/s) '); 
ylabel('Magnitude (dB)');

%Magnitude response of the passband ripple 
figure; 
plot(w, 20*log10(abs(HwnT)));
title('Magnitude Response of the Passband ripples - Kaiser Window'); 
xlim([wc1 wc2]);
axis([wc1, wc2, -0.05, 0.05]); 
grid on;
xlabel('Angular Frequency (rad/s) '); 
ylabel('Magnitude (dB)');




%Rectangular window - Time Domain 
figure; 
stem((0:N-1), h_nT);
grid on;
title('Impulse response - Rectangular Window');
ylabel('Amplitude');
xlabel('n')

%Rectangular window - Frequency Domain 
w_Rec=freqz(h_nT);
figure; plot(w,20*log10(w_Rec)); 
grid on;
axis([0, 2100, -140, 15]);
title('Magnitude Response in the range (0,w s/2) - Rectangular window '); 
ylabel('Magnitude (dB)');
xlabel('Angular Frequency (rad/s) ')

%Comparing the filters	in Frequency Domain using Kaiser & Rectangular windows 
figure;
plot(w, 20*log10(abs(HwnT))); 
hold on;
plot(w,20*log10(w_Rec)); 
axis([0, 2100, -140, 15]);
title('Comparing the filters in Frequency Domain using Kaiser & Rectangular windows'); 
legend('Kaiser Window','Rectangular Window');
ylabel('Magnitude (dB)'); 
xlabel('Angular Frequency (rad/s) ')

figure;
plot(w, 20*log10(abs(HwnT))); 
hold on, 
plot(w,20*log10(w_Rec)); 
axis([wc1, wc2, -1, 1]);
title('Comparing passband Ripples '); 
legend('Kaiser Window','Rectangular Window'); 
ylabel('Magnitude (dB)');
xlabel('Angular Frequency (rad/s) ')




%Input signal specifications 
samples = 600;
w1 = wc1/2 ;
w2 = (wc1+wc2)/2 ;
w3 = (ws/2 +wc2)/2 ;
n_dis = 0:1:samples;
x_nT = sin(w1*n_dis*T)+sin(w2*n_dis*T)+sin(w3*n_dis*T);

fft_len = length(x_nT)+length(hwnT)-1; % length for fft in x dimension 
x_fft = fft(x_nT,fft_len);
hw_fft = fft(hwnT,fft_len');
y_fft = hw_fft.*x_fft; % get the filter output in frequency domain 
y_ifft = ifft(y_fft,fft_len); % take Inverse Fast Fourier Transform
y_delay = y_ifft(floor(N/2)+1:length(y_ifft)-floor(N/2));% delay the output by shifting

% Input signal Frequency Domain Representation 
figure;
subplot(2,1,1); 
fft_len = 1023;
x_fft = fft(x_nT,fft_len);

% construct the input signal in frequency domain using fft values
x1_fft = [abs(x_fft(fft_len/2+1:fft_len)),abs(x_fft(1)),abs(x_fft(2:fft_len/2+1))]; 
w_len = ws*linspace(0,1,fft_len)-ws/2;
plot(w_len,x1_fft);
title('Input signal Frequency Domain Representation'); 
xlabel('Frequency (rad/s) ')
ylabel('Magnitude')

% Input signal Time Domain Representation 
subplot(2,1,2)
stem(n_dis,x_nT) 
axis([0, 75, -3, 3])
title('Input signal Time Domain Representation'); 
xlabel('n')
ylabel('Amplitude') 
hold on
n_cont = 0:0.001:samples;
x_t = sin(w1*n_cont*T)+sin(w2*n_cont*T)+sin(w3*n_cont*T); 
plot(n_cont,x_t,"r--")
legend('Input signal','sin(575t)+sin(1400t)+sin(1875t) ');

% Output signal Frequency Domain Representation 
figure
subplot(2,1,1) 
fft_len = 1023;
w_len = ws*linspace(0,1,fft_len)-ws/2; 
y_delay_fft = fft(y_delay,fft_len);

% construct the output signal in frequency domain using fliter output fft values 
y1_delay_fft =[abs(y_delay_fft(fft_len/2+1:fft_len)),abs(y_delay_fft(1)),abs(y_delay_fft(2:fft_len/2+1))]; plot(w_len,y1_delay_fft);
title('Output signal Frequency Domain Representation'); xlabel('Frequency (rad/s) ')
ylabel('Magnitude')

subplot(2,1,2) 
ideal_out=zeros(1,fft_len); 
for i = 1:fft_len
    ideal_out(i)=sin(w2*i*T);
end
ideal_fft= fft(ideal_out,fft_len);
ideal1_fft =[abs(ideal_fft(fft_len/2+1:fft_len)),abs(ideal_fft(1)),abs(ideal_fft(2:fft_len/2+1))];
plot(w_len,ideal1_fft);
title('Expected signal Frequency Domain Representation');
axis([-2500,2500, 0, 550])
xlabel('Frequency (rad/s) ') 
ylabel('Magnitude')

%Output signal Time Domain Representation 
figure;
subplot(2,1,1) 
stem((0:1:718),y_ifft) 
axis([0, 200, -1.5, 1.5])
xlabel('n') 
ylabel('Amplitude')
title("Output signal Time Domain Representation");

subplot(2,1,2) 
n1=1:1:fft_len; 
stem(n1,sin(w2*n1*T));
title('Expected signal Time Domain Representation'); 
xlabel('n');
ylabel('Amplitude');
axis ( [0,200,-1.5,1.5]);

function alpha = get_Alpha(Aa)
    if(Aa<=21)
        alpha=0; 
    elseif(Aa>21 && Aa<=50)
        alpha=0.5842*((Aa-21)^0.4)+0.07886*(Aa-21);
    else
        alpha=0.1101*(Aa-8.7);
    end
 end
 

function D = get_D(Aa)
    if(Aa<=21)
        D=0.9222;
    else
        D=(Aa-7.95)/14.36;
    end
end
 

function N = get_N(D,Bt,ws) 
    N=ceil(((ws*D)/Bt)+1);
    if(mod(N,2) == 0) 
        N=N+1;
    end
end

function I0_x = get_Bessel(x)
        k=1;
        value = (1/factorial(k)*(x/2)^k)^2; 
        I0_x = 1 + value;
        while (value > 10^-6)
            k=k+1;
            value=(1/factorial(k)*(x/2)^k)^2;
            I0_x = I0_x + value;
        end
end



