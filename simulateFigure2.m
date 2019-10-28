%%This Matlab script generates Figure 2 in the paper:
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

%% Set parameter values

%Carrier frequency (in GHz)
fc = 3;

%Distances in m
d = 10:0.1:100;

%Define the antenna gains at the transmitter and receiver
antennaGainTdBi = db2pow(5);
antennaGainRdBi = db2pow(5);

%Compute channel gains based on the 3GPP Urban Micro in "Further
%advancements for E-UTRA physical layer aspects (Release 9)."
%3GPP TS 36.814, Mar. 2010. Note that the antenna gains are included.
beta_3GPP_LOS = db2pow(antennaGainTdBi+antennaGainRdBi-28-20*log10(fc)-22*log10(d));
beta_3GPP_NLOS = db2pow(antennaGainTdBi+antennaGainRdBi-22.7-26*log10(fc)-36.7*log10(d));



%% Plot simulation results

figure;
hold on; box on; grid on;
plot(d,10*log10(beta_3GPP_LOS),'k-.','LineWidth',2);
plot(d,10*log10(beta_3GPP_NLOS),'r--','LineWidth',2);
xlabel('Distance $d$ [m]','Interpreter','Latex');
ylabel('Channel gain $\beta(d)$ [dB]','Interpreter','Latex');
legend({'UMi-LOS','UMi-NLOS'},'Interpreter','Latex');
set(gca,'fontsize',18);
ylim([-110 -50]);
