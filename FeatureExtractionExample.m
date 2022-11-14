% In this script, an example of the feature extraction scheme used in the work is presented. 
% For its execution it is necessary to have the Matlab Audio Toolbox.
load BreathCycleExample.mat %For this example we use a respiratory cyle containing a wheeze, the variable name is T, its sample rate is 4 KHz
%% section 3.1
% For the Fourier Transform the fft function is used.
x=(T-mean(T))./std(T);
X=fft(x);
%% section 3.2 and 3.3
% 


% For the feature extraction scheme 3 types of descriptors are implemented,
% Spectral features, Mel frequency cepstral coefficients and Chroma features

% First we use the spectrogram function, it receives the time signal and allows 
% to choose the size of the windows to be used and the number of overlapping 
% samples between windows.

hs=113;% Overlapped samples
L=125;% Window size
x=(T-mean(T))./std(T);
[s1, fs1]=spectrogram(x,hanning(L),hs,L,4000);

THETA=zeros(size(s1,2),36); %Matrix in wich the features will be stored

[Mfcc]=mfcc(s1,4000,'LogEnergy','ignore'); % The the Mel Frequency Cepstral Coefficients we use the function mfcc
Chroma=ChromaVector(abs(s1).^2).'; % The function used to extract the chroma vector is at the end of this script
THETA(:,1:13)=Mfcc;
THETA(:,14:25)=Chroma;

% The spectral descriptors are described in https://la.mathworks.com/help/audio/ug/spectral-descriptors.html#mw_rtc_SpectralDescriptorsExample_7E5C532B
% and in the book  "An Introduction to Audio Content Analysis Applicationsin Signal Processing and Music Informatics" by Alexander Lerch.
% These characteristics were extracted in 3 frequency bands defined as: Low (60-300Hz), Mid (300-800Hz) and High (800-1800Hz). 

Fi=300;% Fi and FF are the frequency intervals to be worked on, in this example we use those of the mid-band.
Ff=800;
Bi=round(Fi*L/4000);
Bf=round(Ff*L/4000);
s=s1(Bi:Bf,:);
fs=fs1(Bi:Bf);
% The functions used for the spectral descriptors receive the power spectrum, 
% for which s1 is used, which is the variable containing the spectrogram calculated 
% with the desired windowing scheme and the indicated frequency interval.
THETA(:,26)=spectralEntropy(abs(s).^2,fs);
[FSkewness,FSpread,FCentroid]=spectralSkewness(abs(s).^2,fs);
THETA(:,27)=FSkewness;
THETA(:,28)=FCentroid;
THETA(:,29)=FSpread;
THETA(:,30)=spectralSlope(abs(s).^2,fs);
THETA(:,31)=spectralCrest(abs(s).^2,fs);
THETA(:,32)=spectralFlux(abs(s).^2,fs);
THETA(:,33)=spectralKurtosis(abs(s).^2,fs);
THETA(:,34)=spectralFlatness(abs(s).^2,fs);
THETA(:,35)=spectralRolloffPoint(abs(s).^2,fs);
THETA(:,36)=spectralDecrease(abs(s).^2,fs);


% The Descriptive statistics og the section 3.3 are extracted using the
% functions mean and std.
Mean=mean(THETA);
StandardDeviation=std(THETA);
CoefVariance=StandardDeviation./Mean;


% For the Fractional state space reprensentation Grunwald–Letnikov definition 
% of fractional order derivative, the function for its implementation is
% prensented at the end of the script as FracDer.

% In this example we calculate the 0.7 derivative of one of the feauture
% signals
h=12/4000; %Time diferential
alpha=0.7; % Fractional order, in the work this takes values between -1 and 2
n=100; %Size of the buffer for the derivative calculation
D=FracDer(THETA(:,19),h,alpha,n);












function [Chromagram] = ChromaVector(s)
    BW=0.5;
    FrecPrin=[61.41,69.3,73.42,77.78,82.41,87.31,92.50,98,103.83,110,116.54,123.47];
    LowEdge=2^(-BW/12).*FrecPrin;
    HighEdge=2^(BW/12).*FrecPrin;
    Num=size(s,2);
    Chromagram=zeros(12,Num);
    for j=1:Num   
        X=s(:,j); 
        for gg=1:length(FrecPrin)
            for oct=1:5
            BinLow=floor((length(X)./4000)*LowEdge(gg)*2^(oct-1));
            BinHigh=floor((length(X)./4000)*HighEdge(gg)*2^(oct-1));
            if BinHigh>=2000
                BinHigh=1999;
            end
            Chromagram(gg,j)=(mean(X(1+(BinLow:BinHigh))))+Chromagram(gg,j);
            end            
        end
        Chromagram(:,j)=(Chromagram(:,j))./sum(Chromagram(:,j));   
    end
end


function [D] = FracDer(x,h,alpha,n)
h_a = h^alpha;
N=length(x);
D=zeros(100,1);
for i = 1 : N
    sigma = 0;
    cm=1;
    for m = 0 : 1 : n
    if i - m < 1
        break;
    end 
        sigma = sigma + cm* x(i - m);
        cm=cm*(1-(1+alpha)/(m+1));
    end
    D(i) = sigma / h_a;
end

end
