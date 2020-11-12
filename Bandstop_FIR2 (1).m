f_samp = 260e3;

%Band Edge speifications
fp1 = 61.4e3;
fs1 = 65.4e3;
fs2 = 85.4e3;
fp2 = 89.4e3;


Wc1 = fp1*2*pi/f_samp;
Wc2  = fp2*2*pi/f_samp;
Ws1 = fs1*2*pi/f_samp;
Ws2  = fs2*2*pi/f_samp;
%Kaiser paramters
A = -20*log10(0.15);
if(A < 21)
    beta = 0;
elseif(A <51)
    beta = 0.5842*(A-21)^0.4 + 0.07886*(A-21);
else
    beta = 0.1102*(A-8.7);
end

Wn = [(fs1+fp1)/2 (fs2+fp2)/2]*2/f_samp;        %average value of the two paramters
N_min = 1 + ceil((A-8) / (2.285*(Ws1-Wc1)));       %empirical formula for N_min

%Window length for Kaiser Window
n=N_min+13;
disp(n);

%Ideal bandstop impulse response of length "n"

bs_ideal =  ideal_lp(pi,n) -ideal_lp((Ws2+Wc2)/2,n) + ideal_lp((Ws1+Wc1)/2,n);

%Kaiser Window of length "n" with shape paramter beta calculated above
kaiser_win = (kaiser(n,beta))';

FIR_BandStop = bs_ideal .* kaiser_win;
%ssfvtool(FIR_BandStop);         %frequency response

%magnitude response
[H,f] = freqz(FIR_BandStop,1,1024, f_samp);
plot(f,abs(H))
plot(f,abs(H));
yline(0.85,':');
yline(0.15,':');
xline(fp1,':');
xline(fs1,':');
xline(fs2,':');
xline(fp2,':');
grid

function hd = ideal_lp(wc,M)
    alpha = (M-1)/2;
    n = (0:1:(M-1));
    m = n - alpha + eps;
    hd = sin(wc*m) ./ (pi*m);
end