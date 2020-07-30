%% Array Factor

ASingle = load('N:\MASTERS\Quarter 3\Antenna Systems\Matlab\data\PatchC_SingleElmnt_HCO_Y-140.mat');
k01 = 2*pi*5e9/3e8;
th_2 = 0*pi/180; %scanning angle
beta_n = -k01*d*sin(th_2); %Phase shift

dth = pi/180; % resolution 

pos = -pi/2:dth:pi/2;

N = 4; %Number of elements

%Plotting Array factor:

AF = zeros(size(ASingle.positions,1),size(ASingle.positions,2));

for m = 1:N
    AF = AF + 1*exp(1j*(m-1)*(beta_n+ (k01*d*sind(ASingle.positions)))); % slides 
end

AF_norm = abs(AF/max(max(AF)));
test = ASingle.CpxData(:,101).*abs(AF(:,1));

% Gain of array:
test2 =  ASingle.CpxData(91,:).*abs(AF(:,1));
% 
% figure()
% plot(ASingle.frequencies(:,1)*10^-9,pow2db(abs(test2(91,:))));
% 
% 
% figure()
% plot(ASingle.positions, pow2db((abs(AF))),'LineWidth',2);
% grid on;
% title('Array Factor');
% xlabel('Elevation Angle \theta [degrees]');
% ylabel('Array Factor [dB]');
% ylim([-40,10]);

% figure()
% plot(ASingle.positions, pow2db((abs(test(:,1)))),'LineWidth',2);
% hold on;
% 
% %plot(Anet.positions(91:271,1),pow2db(S21_AUT(91:271,1)),'LineWidth',2);
% grid on;
% title('Array Factor mulitplied with single element pattern');
% xlabel('Elevation Angle \theta [degrees]');
% ylabel('Radiation Pattern [dB]');
% legend('ArrayFactor X SingleElement','Combined Radiation');
% ylim([-40,-10]);