%%This Matlab script generates Figure 5 in the paper:
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

%Define the channel gain functions based on the 3GPP Urban Micro in
%"Further advancements for E-UTRA physical layer aspects (Release 9)."
%3GPP TS 36.814, Mar. 2010. Note that x is measured in m and that the
%antenna gains are included later in the code
pathloss_3GPP_LOS = @(x) db2pow(-28-20*log10(fc)-22*log10(x));
pathloss_3GPP_NLOS = @(x) db2pow(-22.7-26*log10(fc)-36.7*log10(x));

%Define the antenna gains at the source, relay/IRS, and destination. The
%numbers are in linear scale
antennaGainS = db2pow(5);
antennaGainR = db2pow(5);
antennaGainD = db2pow(0);

%Set the amplitude reflection coefficient
alpha = 1;

%Set the range of rate values
Rbar = [0.01 0.1:0.1:10];

%Set parameters related to circuit power consumption
Ps = 100; %Power dissipation in the transceiver hardware of the source
Pd = 100; %Power dissipation in the transceiver hardware of the destination
Pe = 5;   %Power dissipation per element in the IRS (mW)
Pr = 100; %Power dissipation in the transceiver hardware of the relay
nu = 0.5; %Efficiency of the power amplifier at the source


%Define distances in simulation setup

d_SR = 80; %Distance between the source and IRS/relay
dv = 10; %Minimum distance between destination and the IRS/relay

%Define the range of d1 values in the simulation setup
d1 = 70;

%Copmute distance between the source and destination
d_SD = sqrt(d1^2+dv^2);

%Compute distance between the IRS/relay and destination
d_RD = sqrt((d1-d_SR)^2+dv^2);


%Compute the channel gains using the 3GPP models and antenna gains
betaSR = pathloss_3GPP_LOS(d_SR)*antennaGainS*antennaGainR;
betaRD = pathloss_3GPP_LOS(d_RD)*antennaGainR*antennaGainD;
betaSD = pathloss_3GPP_NLOS(d_SD)*antennaGainS*antennaGainD;


%Prepare to save simulation results
EE_SISO = zeros(length(Rbar),1);
EE_IRS = zeros(length(Rbar),1);
EE_DF = zeros(length(Rbar),1);
Nopt = zeros(length(Rbar),1);


%% Go through all rate values
for ind = 1:length(Rbar)
    
    %Compute required SINR values
    SINR = 2^(Rbar(ind))-1; %SISO and IRS
    SINR_DF = 2^(2*Rbar(ind))-1; %DF relaying
    
    
    %Compute the transmit power in the SISO case, using Eq. (17)
    P_SISO = SINR*sigma2/betaSD;
    
    %Compute the energy efficiency in the SISO case
    %(the factor 1000 is used to convert mW to W)
    EE_SISO(ind) = 1000*B*Rbar(ind)/(P_SISO/nu + Ps + Pd);
    
    
    %Compute the transmit power in the DF relaying case, using Eq. (19)
    P_DF = SINR_DF*sigma2*(betaSR+betaRD-betaSD)/(2*betaRD*betaSR);
    
    %Compute the energy efficiency in the DF relaying case
    %(the factor 1000 is used to convert mW to W)
    EE_DF(ind) = 1000*B*Rbar(ind)/(P_DF/nu + Ps/2 + Pd + Pr);
    
    
    %Compute the power-minimizing number of reflecting elements
    Nopt(ind) = (2*SINR*sigma2/(alpha^2*betaSR*betaRD*Pe))^(1/3) - sqrt(betaSD/(betaSR*betaRD))/alpha;
    
    if Nopt(ind)<0
        Nopt(ind) = 0;
    end
    
    %Compute the transmit power in the IRS case, using Eq. (18)
    P_IRS = SINR*sigma2./(sqrt(betaSD) + Nopt(ind)*alpha*sqrt(betaSR*betaRD)).^2;
    
    %Compute the energy efficiency in the IRS case
    %(the factor 1000 is used to convert mW to W)
    EE_IRS(ind) = 1000*B*Rbar(ind)/(P_IRS/nu + Ps + Pd + Nopt(ind)*Pe);
    
end



%Plot simulation results
figure;
hold on; box on;
plot(Rbar,EE_DF/1e6,'b-.','LineWidth',2);
plot(Rbar,EE_IRS/1e6,'r-','LineWidth',2);
plot(Rbar,EE_SISO/1e6,'k--','LineWidth',2);
xlabel('Achievable rate [bit/s/Hz]','Interpreter','Latex');
ylabel('Energy efficiency [Mbit/Joule]','Interpreter','Latex');
legend('DF relay','IRS','SISO','Location','NorthWest');
set(gca,'fontsize',18);
