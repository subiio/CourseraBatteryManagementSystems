load E2model
load("E2_DYN_35_P25.mat")
T = 25;
time = DYNData.script1.time(:); deltat = time(2) - time(1);
time = time - time(1);
current = DYNData.script1.current(:);
voltage = DYNData.script1.voltage(:);
soc = DYNData.script1.soc(:);
sochat = zeros(size(soc));
socbound = zeros(size(soc));

SigmaX0 = diag([1e-3 1e-3 1e-2]);
SigmaV = 2e-1;
SigmaW = 1e1;

ekfData = initEKF(voltage(1),T,SigmaX0,SigmaV,SigmaW,model);

hwait = waitbar(0, 'Computing...');
for k = 1:length(voltage)
    vk = voltage(k);
    ik = current(k);
    Tk = T;

    [sochat(k),socbound(k),ekfData] = iterEKF(vk,ik,Tk,deltat,ekfData);
    if mod(k,1000)==0, waitbar(k/length(current),hwait); end;
end
close(hwait);

% figure;
% plot(time/60, 100*sochat, time/60, 100*soc);
% hold on
% plot([time/60; NaN(1,37660); time/60],...
%     [100*(sochat + socbound); NaN(1,37660); 100*(sochat - socbound)]);
% title('SOC estimation using EKF'); xlabel('Time (min)'); ylabel('SOC  (%)');
% legend('Estimate','Truth','Bounds'); grid on
% fprintf('RMS SOC estimation error = %g%%\n',sqrt(mean((100*(soc-sochat)).^2)));
% 
% figure; plot(time/60,100*(soc-sochat)); holdon
% plot([time/60; NaN(37660,1); time/60],[100*socbound; NaN(1,37660); -100*socbound]);
% title('SOC estimation errors using EKF');
% xlabel('Time(min)'); ylabel('SOC error (%)'); ylim([-4 4]);
% legend('Estimation error','Bounds'); grid on
% 
% ind = find(abs(soc-sochat)>socbound);
% fprintf('Percent of time error outside bounds = %g%%\n', ...
%    length(ind)/length(soc)*100);

% 첫 번째 figure
figure;
plot(time/60, 100*sochat, time/60, 100*soc);
hold on;
plot([time/60; NaN; time/60], ...
     [100*(sochat + socbound);  NaN; 100*(sochat - socbound)]);
title('SOC estimation using EKF');
xlabel('Time (min)');
ylabel('SOC  (%)');
legend('Estimate','Truth','Bounds');
grid on;
fprintf('RMS SOC estimation error = %g%%\n', sqrt(mean((100*(soc-sochat)).^2)));

% 두 번째 figure
figure;
plot(time/60, 100*(soc-sochat));
hold on;
plot([time/60; NaN; time/60], ...
     [100*socbound; NaN; -100*socbound]);
title('SOC estimation errors using EKF');
xlabel('Time (min)');
ylabel('SOC error (%)');
ylim([-4 4]);
legend('Estimation error','Bounds');
grid on;

% 인덱스 찾기 및 출력
ind = find(abs(soc-sochat) > socbound);
fprintf('Percent of time error outside bounds = %g%%\n', length(ind) / length(soc) * 100);