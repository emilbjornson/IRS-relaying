%%This Matlab script generates Figure 4 in the paper:
%
%Emil Björnson, Özgecan Özdogan, Erik G. Larsson, “Intelligent Reflecting
%Surface vs. Decode-and-Forward: How Large Surfaces Are Needed to Beat
%Relaying?,” IEEE Wireless Communications Letters, To appear
%
%Download article: https://arxiv.org/pdf/1906.03949
%
%This is version 1.0 (Last edited: 2019-10-28)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%paper as described above.


close all;
clear;


%% Set simulation parameters

%Carrier frequency (in GHz)
fc = 3;

%Bandwidth
B = 10e6;

%Noise figure (in dB)
noiseFiguredB = 10;

%Compute the noise power in dBm
sigma2dBm = -174 + 10*log10(B) + noiseFiguredB; 
sigma2 = db2pow(sigma2dBm);

%Define the antenna gains at the source, relay/IRS, and destination. The
%numbers are in linear scale
antennaGainS = db2pow(5);
antennaGainR = db2pow(5);
antennaGainD = db2pow(0);

%Define the channel gain functions based on the 3GPP Urban Micro in 
%"Further advancements for E-UTRA physical layer aspects (Release 9)." 
%3GPP TS 36.814, Mar. 2010. Note that x is measured in m and that the
%antenna gains are included later in the code
pathloss_3GPP_LOS = @(x) db2pow(-28-20*log10(fc)-22*log10(x));
pathloss_3GPP_NLOS = @(x) db2pow(-22.7-26*log10(fc)-36.7*log10(x));


%Set the amplitude reflection coefficient
alpha = 1;

%Define the range of different number of reflection elements in the IRS
Nrange = [25 50 100 150];

%Set the range of rate values
Rbar = [4 6];


%Define distances in simulation setup

d_SR = 80; %Distance between the source and IRS/relay
dv = 10; %Minimum distance between destination and the IRS/relay

%Define the range of d1 values in the simulation setup
d1range = 40:100;


%Prepare to save simulation results
powerFractionRelay = zeros(length(d1range),length(Rbar));
Nmin = zeros(length(Rbar),1);


%% Go through all rate values
for ind = 1:length(Rbar)
    
    %Prepare to save transmit powers
    P_IRS = zeros(length(d1range),length(Nrange));
    P_DF = zeros(length(d1range),1);
    P_SISO = zeros(length(d1range),1);
    
    %Compute required SINR values
    SINR = 2^(Rbar(ind))-1; %SISO and IRS
    SINR_DF = 2^(2*Rbar(ind))-1; %DF relaying
    
    
    %Go through all values of d1
    for k = 1:length(d1range)
        
        %Extract value of d1
        d1 = d1range(k);
        
        %Compute distance between the source and destination
        d_SD = sqrt(d1^2+dv^2);
        
        %Compute distance between the IRS/relay and destination
        d_RD = sqrt((d1-d_SR)^2+dv^2);

        
        %Compute the channel gains using the 3GPP models and antenna gains
        betaSR = pathloss_3GPP_LOS(d_SR)*antennaGainS*antennaGainR;
        betaRD = pathloss_3GPP_LOS(d_RD)*antennaGainR*antennaGainD;
        betaSD = pathloss_3GPP_NLOS(d_SD)*antennaGainS*antennaGainD;
        
        
        %Compute the transmit power in mW in the SISO case, using Eq. (11)
        P_SISO(k) = SINR*sigma2/betaSD;
        
        
        %Compute the transmit power in mW in the IRS case, using Eq. (12)
        P_IRS(k,:) = SINR*sigma2./(sqrt(betaSD) + Nrange*alpha*sqrt(betaSR*betaRD)).^2;
        
        
        %Compute the transmit power in mW in the DF relaying case, using Eq. (14)
        if betaSR>=betaSD
            P_DF(k) = SINR_DF*sigma2*(betaSR+betaRD-betaSD)/(2*betaRD*betaSR);
            powerFractionRelay(k,ind) = 2*(betaSR-betaSD)/(betaSR+betaRD-betaSD);
        else
            P_DF(k) = SINR_DF*sigma2/betaSD;
            powerFractionRelay(k,ind) = 0;
        end
        

        %Compute the number of reflecting elements needed to get a lower
        %transmit power with the IRS than with DF relaying
        if d1 == d_SR
            
            rho = P_DF(k)/sigma2;
            
            Nmin(ind) = sqrt((sqrt(1+rho*2*betaRD*betaSR/(betaSR+betaRD-betaSD))-1)/(rho*betaSR*betaRD))-sqrt(betaSD/(betaSR*betaRD));
            
        end
        
    end
    
    
    
    %Plot simulation results
    figure;
    hold on; box on;
    plot(d1range,10*log10(P_SISO),'k--','LineWidth',2);
    plot(d1range,10*log10(P_IRS(:,1)),'r-','LineWidth',2);
    plot(d1range,10*log10(P_DF),'b-.','LineWidth',2);
    for n = 2:length(Nrange)
        plot(d1range,10*log10(P_IRS(:,n)),'r-','LineWidth',2);
    end
    xlabel('Distance $d_1$ [m]','Interpreter','Latex');
    ylabel('Transmit power [dBm]','Interpreter','Latex');
    legend('SISO','IRS','DF relay','Location','NorthWest');
    set(gca,'fontsize',18);
    xlim([40 100]);
    
end
