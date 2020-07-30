%Q1
%% Reading s4p reflection coeff and mutual coupling of array
Nport = 4;
FileName = 'N:\MASTERS\Quarter 3\Antenna Systems\Matlab\data\Patch4C S-matrix.s4p';
[S4,Nport4,Frequencies]=readsp2(FileName,Nport);

%% Reading SGH S11 parameter

Nport = 1;
FileName = 'N:\MASTERS\Quarter 3\Antenna Systems\Matlab\data\SGH S11.s1p';
[SGH,Nport1,FrequenciesSGH]=readsp2(FileName,Nport);

%% Calculating impedance and plotting

%Impedances of the patch
%Assume Z0 as 50 Ohm
Z0 = 50;
%Calculating impedances
Zl = zeros([4, size(Frequencies, 2)]);

for ind = 1:4
    Zl(ind,:) = Z0.*(1 + squeeze(S4(:,ind,ind))) ...
        ./(1 - squeeze(S4(:,ind,ind)));
    %Plot reflection coeff
    name = ['S', num2str(ind), num2str(ind)];
    nameZ = ['Real(Zin',num2str(ind), ')'];
    nameZi = ['Imag(Zin',num2str(ind), ')'];
    figure(1);
    plot(Frequencies(1,:)./10^9, mag2db(squeeze(abs(S4(:,ind,ind)))), 'DisplayName', name, 'LineWidth', 1.5); hold on;
    figure(2);
    plot(Frequencies(1,:)./10^9, ((real(Zl(ind,:)))), 'DisplayName', nameZ, 'LineWidth', 1.5); hold on;
    figure(3);
    plot(Frequencies(1,:)./10^9, ((imag(Zl(ind,:)))), 'DisplayName', nameZi, 'LineWidth', 1.5); hold on;
end

%Beautifyin the figures
figure(1);
hold off;
legend show;
title('Return Loss of Antenna Array Elements');
xlabel('Frequencies (in GHz)');
ylabel('Reflection Coefficients (in dB)');

figure(2);
hold off;
legend show;
title('Real Part of Input Impedance of Antenna Array Elements');
xlabel('Frequencies (in GHz)');
ylabel('Input Impedance (in Ohm)');
%ylim([-100, 500])

figure(3);
hold off;
legend show;
title('Imaginary Part of Input Impedance of Antenna Array Elements');
xlabel('Frequencies (in GHz)');
ylabel('Input Impedance (in Ohm)');
%ylim([-300, 500])

%Impedances of the SGH
ZlSGH = Z0.*(1 + (SGH)./(1 - (SGH)));

%Plotting S11, realZ and imagZ for SGH
figure(4);
plot(FrequenciesSGH(1,:)./10^9, mag2db(squeeze(abs(SGH))),'LineWidth', 1.5);
%title('Reflection Coefficient of SGH');
title('Return Loss of SGH');
xlabel('Frequencies (in GHz)');
ylabel('Reflection Coefficient (in dB)');

%Real part of horn antennas
figure(5);
plot(FrequenciesSGH(1,:)./10^9, ((real(ZlSGH))), 'LineWidth', 1.5);
title('Real Part of Input Impedance of SGH');
xlabel('Frequencies (in GHz)');
ylabel('Input Impedance (in Ohm)');

%Imaginary part of the horn antennas
figure(6);
plot(FrequenciesSGH(1,:)./10^9, ((imag(ZlSGH))), 'LineWidth', 1.5);
title('Imag Part of Input Impedance of SGH');
xlabel('Frequencies (in GHz)');
ylabel('Input Impedance (in Ohm)');

%% Coupling between the elements

%S4 = smoothdata(S4);
for ind = 1:4
    for indj = 1:4
        if(ind == indj)
            continue;
        else
            figure(ind+6);
            name = ['S', num2str(ind), num2str(indj)]; 
            plot(Frequencies(1,:)./10^9, pow2db(squeeze(abs(S4(:,ind,indj)))), 'DisplayName', name, 'LineWidth', 1.5); hold on; 
        end
    end
end

%Beautifying the figures
figure(7);
hold off;
legend show;
title('Mutual Coupling Between Elements S11 to S14');
xlabel('Frequency (in GHz)');
ylabel('Sij in dB');

figure(8);
hold off;
legend show;
title('Mutual Coupling Between Elements S21 to S24');
xlabel('Frequency (in GHz)');
ylabel('Sij in dB');

figure(9);
hold off;
legend show;
title('Mutual Coupling Between Elements S31 to S34');
xlabel('Frequency (in GHz)');
ylabel('Sij in dB');

figure(10);
hold off;
legend show;
title('Mutual Coupling Between Elements S41 to S44');
xlabel('Frequency (in GHz)');
ylabel('Sij in dB');


%% Scanning angle and active reflection coefficient

%Resonant frequency
freq = 5e9;

%Scanning angles
drad = pi/180;
th = eps:15*drad:60*drad;
thd = 0:15:60;

%Calculate initial phase shift
c = 3e8;
lam = c/freq;
d = 0.06; %Distance between arraty elements is 6 cm
beta = -2.*pi.*d./lam.*sin(th);
%beta = [0 0 0 0 0];
%betaMesh = meshgrid(beta, beta);

%Define active reflectio coefficient matrix
%S = zeros([4, size(Frequencies, 2)]);
%Alternate S-parameter matrix to calculate active impedance
S4a = zeros(size(S4));
for ind = 1:size(S4, 2)
    S4a(:,ind,:) = [S4(:,ind,ind) squeeze(S4(:,ind,1:ind-1))...
        squeeze(S4(:,ind,ind+1:size(S4,3)))];
end

%Caculating the active impedance
ZinA = zeros([size(Frequencies, 2) size(beta, 2) size(S4a, 2)]);

%ReqMat is the active S11 matrix
reqMat = zeros([size(Frequencies, 2) size(beta, 2) size(S4a, 2)]);

for ind = 1:size(beta, 2)
    N = meshgrid(0:size(S4a,3)-1, 0:size(S4a,3)-1);
    Nb = 1j.*N.*beta(ind);
    for indj = 1:size(S4a,1)
        reqMat1 = sum(exp(Nb).*squeeze(S4a(indj,:,:)), 2);
        reqMat(indj, ind, :) = sum(exp(Nb).*squeeze(S4a(indj,:,:)), 2);
        ZinA(indj, ind, :) = Z0.*(1+reqMat1)./(1-reqMat1);
    end
end 

%Plotting
for ind = 1:4
    figure(ind);
    plot(Frequencies(1,:)./10^9, mag2db((abs(reqMat(:,1,ind)))), 'LineWidth', 1.5, 'DisplayName', "Active-0 deg"); hold on;
    plot(Frequencies(1,:)./10^9, mag2db((abs(reqMat(:,2,ind)))), 'LineWidth', 1.5, 'DisplayName', "Active-15 deg"); hold on;
    plot(Frequencies(1,:)./10^9, mag2db((abs(reqMat(:,3,ind)))), 'LineWidth', 1.5, 'DisplayName', "Active-30 deg"); hold on;
    plot(Frequencies(1,:)./10^9, mag2db((abs(reqMat(:,4,ind)))), 'LineWidth', 1.5, 'DisplayName', "Active-45 deg"); hold on;
    plot(Frequencies(1,:)./10^9, mag2db((abs(reqMat(:,5,ind)))), 'LineWidth', 1.5, 'DisplayName', "Active-60 deg"); hold on;
    plot(Frequencies(1,:)./10^9, mag2db((abs(S4(:,ind,ind)))), 'LineWidth', 1.5, 'DisplayName', "Passive");
    
    title(['Active and Passive reflection coefficients, S', num2str(ind), num2str(ind)]);
    xlabel('Frequency (in GHz)');
    ylabel('Sii (Active and Passive) in dB');
    legend show;
end

%% Q2 Gain and raditaion pattern

%Gain of the patch vs. frequency Combined
%Reading S21 of the SGH
Nport = 2;
FileName = 'N:\MASTERS\Quarter 3\Antenna Systems\Matlab\data\CBand SGH transfer - R apertures 3.27m.s2p';
[S21,Nport2,Frequencies21]=readsp2(FileName,Nport);

%Frequency center (Assume 5GHz)
%freq = 5e9;

%Power calculation PSGH
%idx = find(Frequencies21 == 5e9);
PSGH = mag2db(abs(S21(:,2,1)));

%Gain calculation Using the manual
G = 5.0195 + 4.2742.*(Frequencies21./10^9) - 0.3062.*(Frequencies21./10^9).^2;

%Gain in dB
Gdb = pow2db(G);

%Gain of the patch at 0 deg
%load data
figure(1)
for ind = 1:4
    load(['N:\MASTERS\Quarter 3\Antenna Systems\Matlab\data\Patch', num2str(ind),'_HCO.mat']);
    PAUT = mag2db(abs(CpxData((positions == thd(1)),:)));
    GAUT = Gdb - PSGH' + PAUT;
    name =['Gpatch', num2str(1)];
    plot(frequencies(:,1)./10^9, GAUT(1,:), 'LineWidth', 1.5, 'DisplayName', name); hold on;
end

load(['N:\MASTERS\Quarter 3\Antenna Systems\Matlab\data\Combined_HCO.mat']);
PAUT = mag2db(abs(CpxData((positions == thd(1)),:)));
GAUT = Gdb - PSGH' + PAUT;

%Gain of the SGH
% figure(12);
% plot(Frequencies21(1,:), Gdb(1,:));
% title('Gain (in dB) Vs. Frequency Plot of SGH');
% xlabel('Frequency (in GHz)');
% ylabel('Gain in (dB)');

%Combined gain of the patches
figure(1);
hold off;
legend show;
title('Gain (in dB) Vs. Frequency Plot of Patch, Scan Angle = 0 deg');
xlabel('Frequency (in GHz)');
ylabel('Gain in (dB)');

%% Individual gains of the patches

for ind = 1:4
    FileName = strcat('N:\MASTERS\Quarter 3\Antenna Systems\Matlab\data\Patch',num2str(ind),'_HCO.mat');
    load(FileName);
    PAUT = mag2db(abs(CpxData((positions == thd(1)),:)));
    GAUT_ind = Gdb - PSGH' + PAUT;
    figure(13+ind);
    plot(frequencies(:,1), GAUT_ind(1,:));
    title(strcat('Gain (in dB) Patch', num2str(2), 'Vs. Frequency Plot, \theta = 0'));
    xlabel('Frequency (in GHz)');
    ylabel('Gain in (dB)');    
end

%% Comparing with the theory

%Gain of single element
%load('N:\MASTERS\Quarter 3\Antenna Systems\Matlab\data\PatchC_SingleElmnt_HCO_Y-140.mat');
load('N:\MASTERS\Quarter 3\Antenna Systems\Matlab\data\Patch3_HCO.mat');
PAUT = mag2db(abs(CpxData((positions == thd(1)),:)));

GAUT_single = Gdb - PSGH' + PAUT;

k = 2.*pi.*frequencies./c;
AF = zeros([size(th, 2) size(frequencies, 1)]);
Narray = 4;

%AF Calculations
for ind = 1:size(beta, 2)
    Psi=k.*d.*sin(th(ind))+beta(ind);
    AF(ind, :) =(1/Narray).*sin((Narray/2).*Psi)./sin((1/2).*Psi);    
end

%At position 0; Theoretical gain
Gt = GAUT_single + mag2db(AF(1, :));

%Plotting the comparison
figure(18);
plot(frequencies(:,1)./(10^9), GAUT(1,:), 'LineWidth', 1.5); hold on;
plot(frequencies(:,1)./(10^9), Gt(1,:), 'LineWidth', 1.5);
title('Theoretical and Practical Gain value comparison, Scan Angle = 0 deg');
xlabel('Frequency (in GHz)');
ylabel('Gain in (dB)');
legend('Practical', 'Theoretical');

%% Embedded pattern

%Load all the files
for ind = 1:4
    FileName = strcat('N:\MASTERS\Quarter 3\Antenna Systems\Matlab\data\Patch',num2str(ind),'_HCO.mat');
    temp_var = strcat('patch_',num2str(ind));
    Pattern.(temp_var) = load(FileName);
end

%For different beta, Array pattern calculation
S21array = zeros([size(Pattern.patch_1.CpxData) size(beta, 2)]);
% 
% for ind = 1:size(beta, 2)
%     S21array(:,ind) = (Pattern.patch_1.CpxData((Pattern.patch_1.positions == thd(ind)), :) ...
%         + Pattern.patch_2.CpxData((Pattern.patch_2.positions == thd(ind)), :).*exp(1j*beta(ind)) ...
%         + Pattern.patch_3.CpxData((Pattern.patch_3.positions == thd(ind)), :).*exp(2j*beta(ind)) ...
%         + Pattern.patch_4.CpxData((Pattern.patch_4.positions == thd(ind)), :).*exp(3j*beta(ind)))';
%     figure(18+ind);
%     plot(Pattern.patch_1.frequencies(:,1)./(10^9), mag2db(abs(S21array(:,ind))), 'LineWidth', 1.5);
%     %title('Array Pattern of Patch' + num2str(ind));
%     xlabel('Frequencies (in GHz)');
%     ylabel('Array Pattern (in dB)');
% end
%load('N:\MASTERS\Quarter 3\Antenna Systems\Matlab\data\PatchC_SingleElmnt_HCO_Y-140.mat');
load('N:\MASTERS\Quarter 3\Antenna Systems\Matlab\data\Combined_HCO.mat');
%Plotting with respect to theta
for ind = 1:size(beta, 2)
    S21array(:,:,ind) = (Pattern.patch_1.CpxData ...
        + Pattern.patch_2.CpxData.*exp(1j*beta(ind)) ...
        + Pattern.patch_3.CpxData.*exp(2j*beta(ind)) ...
        + Pattern.patch_4.CpxData.*exp(3j*beta(ind)));
    AF_Req = AFcode(freq, th(ind), 0.06);
    
    figure(18+ind);
    plot(Pattern.patch_1.positions(:,1), mag2db(abs(S21array(:,(Pattern.patch_2.frequencies == freq),ind))), 'LineWidth', 1.5); hold on;
    plot(positions(:,1), mag2db((abs(CpxData(:,101)))),'LineWidth',1.5);
    %plot(ASingle.positions, mag2db((abs(test(:,1)))),'LineWidth',1.5);
    title(['Array Pattern Theory and Practical, scan angle',num2str(thd(ind))]);
    xlabel('\theta (in deg)');
    ylabel('Array Pattern (in dB)');
    legend('Computed', 'Combined');
    ylim([-60, -10]);
end

%% Plotting embedded patterns

figure(24);
%freq = 5.8e9;
for ind = 1:4
    temp_var = strcat('patch_',num2str(ind));
    %figure(23+ind);
    name = ['Patch_', num2str(ind)];
    plot(Pattern.patch_1.positions(:,1), mag2db(abs(Pattern.(temp_var).CpxData(:,(Pattern.patch_2.frequencies == freq)))), 'DisplayName', name, 'LineWidth', 1.5); hold on;
    %title('Array Pattern of Patch' + num2str(ind));
end
singleElem = load('N:\MASTERS\Quarter 3\Antenna Systems\Matlab\data\PatchC_SingleElmnt_HCO_Y-140.mat');
name = ['Single Element'];
plot(singleElem.positions(:,1), mag2db(abs(singleElem.CpxData(:,(singleElem.frequencies == freq)))), 'DisplayName', name, 'LineWidth', 1.5);
title('Emdedded pattern of each element and single element, freq = 5GHz');
xlabel('\theta (deg)');
ylabel('Array Pattern (in dB)');
hold off;
legend show;

%% AF calculation and comparison with combined pattern

%Theroy array pattern
%APth = zeros([size(singleElem.CpxData, 1)]);

%k0 = 2.*pi.*freq./c;
%AFreq = zeros([size(size(singleElem.CpxData, 1))]);
%beta_dash = -2.*pi.*d.*freq./c.*sin(singleElem.positions);
%Psi_dash  = k0.*d.*sin(singleElem.positions);
%AF_req = (1/Narray).*sin((Narray/2).*Psi_dash)./sin((1/2).*Psi_dash);
%Narray = 4;

Pattern_Req = zeros([size(singleElem.CpxData) size(AF, 1)]);

for ind = 1:size(AF, 1)
    for indj = 1:size(Pattern_Req, 1)
        Pattern_Req(indj,:,ind) = (AF(ind, :).*singleElem.CpxData(indj,:));
    end
    disp(AF(ind,1));
end

plot(singleElem.positions(:,1), mag2db(abs(Pattern_Req(:,singleElem.frequencies == freq,5))), ...
    'LineWidth', 1.5); hold on;
plot(singleElem.positions(:,1), mag2db(abs(singleElem.CpxData(:,singleElem.frequencies == freq))), ...
    'LineWidth', 1.5); 
title('Array Pattern Theory vs. Practical, \theta = 0, freq = 5 GHz');
xlabel('\theta');
ylabel('Array pattern (in dB)');
legend('1', '2');
