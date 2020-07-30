%% Array factor calculation

function [AF] = AFcode(freq, scan_angle, d)
k01 = 2*pi*freq/3e8;
th_2 = scan_angle; %scanning angle
beta_n = k01*d*sin(th_2); %Phase shift

dth = pi/180; % resolution 

pos = -pi/2:dth:pi/2;

N = 4; %Number of elements

%Plotting Array factor:

AF = zeros(size(pos,1),size(pos,2));

for m = 1:N
    AF = AF + 1*exp(1j*(m-1)*(beta_n+ (k01*d*sin(pos)))); % slides 
end

% AF_norm = abs(AF/max(max(AF)));
% %test = ASingle.CpxData(:,101).*abs(AF(1,:));
% 
% figure()
% plot(pos*180/pi, 20*log10((abs(AF))),'LineWidth',2);
% grid on;
% title('Array Factor');
% xlabel('Elevation Angle \theta [degrees]');
% ylabel('Array Factor [dB]');
end 