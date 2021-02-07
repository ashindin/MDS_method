% Программа для чтения и обработки исходного массива данных. Часть 1.
% - гетеродинирование
% - децимация 
% - фильтрация и постобработка
% - построение разностной фазы
% - сравнение с аналитическими данными

%---------------------------------------------------------------------------------------------------------------------------------
%----------------------------------- Чтение полных массивов данных, гетеродинирование и децимация --------------------------------

Data_ID = fopen('../03_FDTD_simulation/E_out_col100.bin');
Data = fread(Data_ID,'double');
fclose(Data_ID);

Data_ID = fopen('../03_FDTD_simulation/cos.bin');
cos = fread(Data_ID,'double');
fclose(Data_ID);

Data_ID = fopen('../03_FDTD_simulation/sin.bin');
sin = fread(Data_ID,'double');
fclose(Data_ID);

%-------------------------------------------------- Параметры исходного сигнала --------------------------------------------------
% Исходный сигнал (данные) представляет собой склейку из 600 идентичных кадров, содержащих излучаемый импульс и
% отраженный. 

% Частота (сэмплинг): 125 МГц

% Размер кадра в отсчетах
Frame_size = 625000;

% Количество отсчетов в сигнале
N_init = length(Data);

% Исходный сэмплинг
sample_rate_init = 125000000;

% Количество кадров в данных
Frame_number = N_init/625000;

% Время симуляции кадра (сек)
time_simul = Frame_size/sample_rate_init;

% Временная шкала в секундах для одного кадра
time_init = zeros(Frame_size,1);

for n = 1:Frame_size
    time_init(n,1) = (n-1)*(10^(-6))/125;   
end

% Строим массив исходных данных, соответствующий первому кадру
Data_1 = Data(1:625000,1);
figure;
plot(time_init,Data_1);
%plot(Data)


% Частотная шкала в МГц
freq_init = zeros(Frame_size,1);

for n = 1:Frame_size
    freq_init(n,1) = ((n-1)*(1/(time_simul)))/1000000;
end


% Строим спектр исходных данных
Spectrum_init = fft(Data_1);
% figure;
% plot(freq_init,abs(Spectrum_init));


% Период следования импульсов 
tau_impulse = 0.1;

% Длительность импульса в секундах
tau = 0.0001;




%----------------------------------------------------- Real to Complex -----------------------------------------------------------

% Домножаем на комплексную экспоненту на частоте 4.5 МГц

Data = Data.*(cos-1j*sin);

clear cos sin;

% Частота комплексной экспоненты
f = 4500000;



% Строим полученные массив данных после преобразования
%figure;
%plot(time_init,real(Data_conv));

% Строим спектр результата преобразования, соответствующий первому кадру
Data_conv_1 = Data(1:625000,1);
Spectrum_conv_1 = fft(Data_conv_1);
figure;
plot(freq_init,abs(Spectrum_conv_1));

% Преобразуем одномерный массив данных в двумерный. Один столбец - один кадр.
%Data = reshape(Data,625000,Frame_number);



%-------------------------------------------------------- Децимация --------------------------------------------------------------
% Понижаем частоту в 125 раз (до 1 МГц)

decim_rate = 125;
sample_rate_decim = 1000000;

% Децимацию придется делать по кускам

Data_1 = Data(1:625000*201,1);
Data_2 = Data(625000*201+1:625000*401,1);
Data_3 = Data(625000*401+1:625000*601,1);

clear Data;

Data_1 = decimate(Data_1,5);
Data_1 = decimate(Data_1,5);
Data_1 = decimate(Data_1,5);

Data_2 = decimate(Data_2,5);
Data_2 = decimate(Data_2,5);
Data_2 = decimate(Data_2,5);

Data_3 = decimate(Data_3,5);
Data_3 = decimate(Data_3,5);
Data_3 = decimate(Data_3,5);

% Data = decimate(Data,5);
% Data = decimate(Data,5);
% Data = decimate(Data,5);

Data_decim = [Data_1;Data_2;Data_3];

clear Data_1 Data_2 Data_3;

dlmwrite('Data_decimated.txt',Data_decim);

% Количество отсчетов после децимации
N_decim = N_init/decim_rate;

% Размер кадра в отсчетах после децимации
Frame_size_decim = Frame_size/decim_rate;


% Временная шкала в секунда для одного кадра
time_decim = zeros(Frame_size_decim,1);

for n = 1:Frame_size_decim
    time_decim(n,1) = (n-1)*(10^(-6));   
end

% Строим полученный массив исходных данных после децимации, соответствующий первому кадру
Data_decim_1 = Data_decim(1:Frame_size_decim,1);
figure;
plot(time_decim,real(Data_decim_1));

% Частотная шкала в МГц
freq_decim = zeros(Frame_size_decim,1);

for n = 1:Frame_size_decim
    freq_decim(n,1) = ((n-1)*(1/(time_simul)))/sample_rate_decim ;
end


% Строим спектр результата децимации
Spectrum_decim = fft(Data_decim_1);
figure;
plot(freq_decim,abs(Spectrum_decim));

%---------------------------------------------------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------------------------------------------------

%---------------------------------------------------------------------------------------------------------------------------------
%--------------------------------------------- Фильтрация и постобработка (FFT и Bessel) -----------------------------------------

% Cэмплинг после децимации
sample_rate_decim = 1000000;

% Количество кадров в данных (соответствует количеству отраженных импульсов)
Frame_number = 601;

% Преобразуем одномерный массив данных в двумерный. Один столбец - один кадр.
Data_decim = reshape(Data_decim,5000,Frame_number);

% Зададим параметры окна, содержащего только отраженный импульс, для обработки.

% Размер окна обработки в секундах
tau_proc = 0.004;

% Размер окна обработки в отсчетах
tau_proc_samp = 0.004*sample_rate_decim;
N_proc = 4001; % тип integer

% Период импульсов в отсчетах
tau_impulse_samp = 0.1*sample_rate_decim;

% impulse_1 = Data_decim(1000:1000+tau_proc_samp,1);
% impulse_2 = Data_decim(1000 + tau_impulse_samp:1000 + tau_impulse_samp+tau_proc_samp,1);
% impulse_3 = Data_decim(1000 + 2*tau_impulse_samp:1000 + 2*tau_impulse_samp+tau_proc_samp,1);



% Создаем массив, содержащий только данные, попадающие в окно обработки каждого кадра
impulse = zeros(N_proc,Frame_number);

for i = 1:Frame_number
 impulse(:,i) = Data_decim(500:500+tau_proc_samp,i);
end

% Cтроим сигнал на выбранных участках
% figure;
% plot(real(impulse_1));
% figure;
% plot(real(impulse_2));
% figure;
% plot(real(impulse_3));


% Находим спектр сигнала в каждом окне обработки
spectrum = zeros(N_proc,Frame_number);

for i = 1:Frame_number
   spectrum(:,i) = fftshift(fft(impulse(:,i))); 
end



% for n = 1:N_proc
%     if  abs(spectrum_1(n,1))== 0  
%         spectrum_1(n,1) = 0;  
%     end
%     if  abs(spectrum_2(n,1))== 0  
%         spectrum_2(n,1) = 0;  
%     end
%     if  abs(spectrum_3(n,1))== 0  
%         spectrum_3(n,1) = 0;  
%     end
% end


% Рассчитываем массив разностных фаз (фазовый набег между следующим импульсом и предыдущим)

phase_dif_fft = zeros(N_proc,Frame_number);

phase_dif_fft(:,1) = unwrap(unwrap(angle(spectrum(:,1))) - unwrap(angle(spectrum(:,1))));

for i = 2:601
   phase_dif_fft(:,i) = unwrap(unwrap(angle(spectrum(:,i))) - unwrap(angle(spectrum(:,i-1))));
end



%dif_1 = unwrap(unwrap(angle(spectrum_2)) - unwrap(angle(spectrum_1)));
%dif_2 = unwrap(unwrap(angle(spectrum_3)) - unwrap(angle(spectrum_2)));


% Частотная шкала в кГц со смещением (отрицательными частотами) на интервале (окне) обработки
freq_decim_shifted = zeros(N_proc,1);

for n = 1:N_proc
    freq_decim_shifted(n,1) = ((n-1)*(1/(0.004)))/1000;
end

freq_decim_shifted = freq_decim_shifted - freq_decim_shifted(int32(((N_proc-1)/2)+1),1);


% Строим разностную фазу
% figure;
% plot(freq_decim_shifted,abs(spectrum(:,2)));
% figure;
% plot(freq_decim_shifted,phase_dif(:,2),'black');
% hold on
% plot(freq_decim_shifted,phase_dif(:,3),'red');
% 
% hold off




%--------------------------------------------------- Старый алгоритм на Lab View -------------------------------------------------

freq_lv = zeros(N_proc,1);

for n = 1:N_proc
    freq_lv(n,1) = 2*pi*((n-1)*(1/N_proc));
end

% Сетка времени в секундах 
time_lv = zeros(N_proc,1);

for n = 1:N_proc
    time_lv(n,1) = (n);
end



% Симулируем работу гетеродина LV

spectrum_phase_lv = zeros(N_proc,Frame_number);

% spectrum_lv_r_1 = zeros(N_proc,1);
% spectrum_phase_lv_r_1 = zeros(N_proc,1);
% 
% spectrum_lv_r_2 = zeros(N_proc,1);
% spectrum_phase_lv_r_2 = zeros(N_proc,1);
% 
% spectrum_lv_r_3 = zeros(N_proc,1);
% spectrum_phase_lv_r_3 = zeros(N_proc,1);

% Зададим параметры НЧ фильтра

%------------------------------------------------------- FIR-фильтр --------------------------------------------------------------
 %pb_freq = 2*pi*5/(10^6);
% sb_freq = 2*pi*15/(10^6);
% 
% 
 %lpFilt = designfilt('lowpassfir','PassbandFrequency',pb_freq, ...
        %  'StopbandFrequency',sb_freq,'PassbandRipple',0.5, ...
         % 'StopbandAttenuation',65,'DesignMethod','kaiserwin');
%fvtool(lpFilt)

%------------------------------------------------------- IIR-фильтр --------------------------------------------------------------
%  lpFilt = designfilt('lowpassiir','FilterOrder',2, ...
%           'PassbandFrequency',50,'PassbandRipple',0.2, ...
%           'SampleRate',1e6);
%       
%       lpFilt2 = designfilt('lowpassiir','FilterOrder',2, ...
%           'PassbandFrequency',10,'PassbandRipple',0.2, ...
%           'SampleRate',1e6);
      
      
%----------------------------------------------------------- Bessel --------------------------------------------------------------
 Wo = 2*pi*25;
 Fs = 10^6;
 [z,p,k] = besself(5,Wo);
 [zd,pd,kd] = bilinear(z,p,k,Fs);    % Analog to digital mapping
 sos = zp2sos(zd,pd,kd);      




% Находим фазу по старому методу
for i = 1:Frame_number
    disp(i);
    for n = 1:N_proc
     % spectrum_lv_r_1(n,1) = sum(impulse_1.*exp(-1j*freq_lv(n,1)*time_lv));
     % spectrum_lv_r_2(n,1) = sum(impulse_2.*exp(-1j*freq_lv(n,1)*time_lv));
     % spectrum_lv_r_3(n,1) = sum(impulse_3.*exp(-1j*freq_lv(n,1)*time_lv));
   
    % spectrum_lv_r_1(n,1) = sum(impulse_1.*exp(-1j*freq_lv(n,1)*time_lv));
    
    % spectrum_phase_lv_r_1(n,1) = sum(unwrap(angle(filter(lpFilt,impulse_1.*exp(-1j*freq_lv(n,1)*time_lv)))))/N_proc;
    % spectrum_phase_lv_r_2(n,1) = sum(unwrap(angle(filter(lpFilt,impulse_2.*exp(-1j*freq_lv(n,1)*time_lv)))))/N_proc;
    % spectrum_phase_lv_r_3(n,1) = sum(unwrap(angle(filter(lpFilt,impulse_3.*exp(-1j*freq_lv(n,1)*time_lv)))))/N_proc;
    
     %spectrum_phase_lv_r_1(n,1) = sum(unwrap(angle(filtfilt(lpFilt,impulse_1.*exp(-1j*freq_lv(n,1)*time_lv)))))/N_proc;
    % spectrum_phase_lv_r_2(n,1) = sum(unwrap(angle(filtfilt(lpFilt,impulse_2.*exp(-1j*freq_lv(n,1)*time_lv)))))/N_proc;
    % spectrum_phase_lv_r_3(n,1) = sum(unwrap(angle(filtfilt(lpFilt,impulse_3.*exp(-1j*freq_lv(n,1)*time_lv)))))/N_proc;
    
     % spectrum_phase_lv_r_1(n,1) = sum((angle      (filtfilt(lpFilt2,(filtfilt(lpFilt,impulse_1.*exp(-1j*freq_lv(n,1)*time_lv)))))     ))/N_proc;
     % spectrum_phase_lv_r_2(n,1) = sum((angle      (filtfilt(lpFilt2,(filtfilt(lpFilt,impulse_2.*exp(-1j*freq_lv(n,1)*time_lv)))))     ))/N_proc;
     % spectrum_phase_lv_r_3(n,1) = sum((angle      (filtfilt(lpFilt2,(filtfilt(lpFilt,impulse_3.*exp(-1j*freq_lv(n,1)*time_lv)))))     ))/N_proc;
   
   
   %-------------------------------------------------------------- For Bessel ----------------------------------------------------
    
    spectrum_phase_lv(n,i) = sum(unwrap(angle      (sosfilt(sos,flip(sosfilt(sos,impulse(:,i).*exp(-1j*freq_lv(n,1)*time_lv)))))     ))/N_proc;
  
    
   
    end
end

% for n = 1:N_proc
%    % spectrum_lv_r_1(n,1) = sum(impulse_1.*exp(-1j*freq_lv(n,1)*time_lv));
%    % spectrum_lv_r_2(n,1) = sum(impulse_2.*exp(-1j*freq_lv(n,1)*time_lv));
%    % spectrum_lv_r_3(n,1) = sum(impulse_3.*exp(-1j*freq_lv(n,1)*time_lv));
%    
%   % spectrum_lv_r_1(n,1) = sum(impulse_1.*exp(-1j*freq_lv(n,1)*time_lv));
%     
%    % spectrum_phase_lv_r_1(n,1) = sum(unwrap(angle(filter(lpFilt,impulse_1.*exp(-1j*freq_lv(n,1)*time_lv)))))/N_proc;
%    % spectrum_phase_lv_r_2(n,1) = sum(unwrap(angle(filter(lpFilt,impulse_2.*exp(-1j*freq_lv(n,1)*time_lv)))))/N_proc;
%    % spectrum_phase_lv_r_3(n,1) = sum(unwrap(angle(filter(lpFilt,impulse_3.*exp(-1j*freq_lv(n,1)*time_lv)))))/N_proc;
%     
%     %spectrum_phase_lv_r_1(n,1) = sum(unwrap(angle(filtfilt(lpFilt,impulse_1.*exp(-1j*freq_lv(n,1)*time_lv)))))/N_proc;
%    % spectrum_phase_lv_r_2(n,1) = sum(unwrap(angle(filtfilt(lpFilt,impulse_2.*exp(-1j*freq_lv(n,1)*time_lv)))))/N_proc;
%    % spectrum_phase_lv_r_3(n,1) = sum(unwrap(angle(filtfilt(lpFilt,impulse_3.*exp(-1j*freq_lv(n,1)*time_lv)))))/N_proc;
%     
%    % spectrum_phase_lv_r_1(n,1) = sum((angle      (filtfilt(lpFilt2,(filtfilt(lpFilt,impulse_1.*exp(-1j*freq_lv(n,1)*time_lv)))))     ))/N_proc;
%     %spectrum_phase_lv_r_2(n,1) = sum((angle      (filtfilt(lpFilt2,(filtfilt(lpFilt,impulse_2.*exp(-1j*freq_lv(n,1)*time_lv)))))     ))/N_proc;
%    % spectrum_phase_lv_r_3(n,1) = sum((angle      (filtfilt(lpFilt2,(filtfilt(lpFilt,impulse_3.*exp(-1j*freq_lv(n,1)*time_lv)))))     ))/N_proc;
%    
%    
%    %-------------------------------------------------------------- For Bessel ----------------------------------------------------
%     
% %     spectrum_phase_lv_r_1(n,1) = sum(unwrap(angle      (sosfilt(sos,flip(sosfilt(sos,impulse(:,1).*exp(-1j*freq_lv(n,1)*time_lv)))))     ))/N_proc;
% %     spectrum_phase_lv_r_2(n,1) = sum(unwrap(angle      (sosfilt(sos,flip(sosfilt(sos,impulse(:,2).*exp(-1j*freq_lv(n,1)*time_lv)))))     ))/N_proc;
% %     spectrum_phase_lv_r_3(n,1) = sum(unwrap(angle      (sosfilt(sos,flip(sosfilt(sos,impulse(:,3).*exp(-1j*freq_lv(n,1)*time_lv)))))     ))/N_proc;
%     
%     spectrum_phase_lv(n,1) = sum(unwrap(angle      (sosfilt(sos,flip(sosfilt(sos,impulse(:,1).*exp(-1j*freq_lv(n,1)*time_lv)))))     ))/N_proc;
%     spectrum_phase_lv(n,2) = sum(unwrap(angle      (sosfilt(sos,flip(sosfilt(sos,impulse(:,2).*exp(-1j*freq_lv(n,1)*time_lv)))))     ))/N_proc;
%     spectrum_phase_lv(n,3) = sum(unwrap(angle      (sosfilt(sos,flip(sosfilt(sos,impulse(:,3).*exp(-1j*freq_lv(n,1)*time_lv)))))     ))/N_proc;
%    
% end




% Находим разностную фазу

%dif_lv_1 = angle(spectrum_lv_r_2) - angle(spectrum_lv_r_1);
%dif_lv_2 = angle(spectrum_lv_r_3) - angle(spectrum_lv_r_1);

%dif_lv_1 = unwrap(unwrap(unwrap(unwrap(fftshift(spectrum_phase_lv_r_2)) - unwrap(fftshift(spectrum_phase_lv_r_1)))));
%dif_lv_2 = unwrap(unwrap(unwrap(unwrap(fftshift(spectrum_phase_lv_r_3)) - unwrap(fftshift(spectrum_phase_lv_r_2)))));

%dif_lv_1 = (((unwrap(spectrum_phase_lv_r_2) - unwrap(spectrum_phase_lv_r_1))));
%dif_lv_2 = (((unwrap(spectrum_phase_lv_r_3) - unwrap(spectrum_phase_lv_r_1))));

phase_dif_lv = zeros(N_proc,Frame_number);

phase_dif_lv(:,1) = unwrap(unwrap(unwrap(unwrap(fftshift(spectrum_phase_lv(:,1))) - unwrap(fftshift(spectrum_phase_lv(:,1))))));

for i = 2:601
   phase_dif_lv(:,i) = unwrap(unwrap(unwrap(unwrap(fftshift(spectrum_phase_lv(:,i))) - unwrap(fftshift(spectrum_phase_lv(:,i-1))))));
end



% Строим разностную фазу
% figure;
% plot(freq_decim_shifted,phase_dif_lv(:,2),'black');
% 
% hold on
% plot(freq_decim_shifted,phase_dif_lv(:,3),'red');


%------------------------------------------------------- Data Converting ---------------------------------------------------------

% Полученные данные имеют шаг по частоте 250 Гц. Попробуем через децимацию сделать 1кГц

% Преобразование массивов данных

phase_dif_fft_1kHz_step = zeros(1001,Frame_number);
phase_dif_lv_1kHz_step = zeros(1001,Frame_number);

for i = 1:Frame_number
    phase_dif_fft_1kHz_step(:,i) = decimate(phase_dif_fft(:,i),4);
    phase_dif_lv_1kHz_step(:,i) = decimate(phase_dif_lv(:,i),4);
end

% dif_lv_1_1kHz = decimate(dif_lv_1,4);
% dif_lv_2_1kHz = decimate(dif_lv_2,4);

freq_decim2_shifted = zeros(1001,1);

for n = 1:1001
    freq_decim2_shifted(n,1) = n;
end

freq_decim2_shifted = freq_decim2_shifted - freq_decim2_shifted(501,1);


% figure;
% plot(freq_decim2_shifted,phase_dif_lv_1kHz_step(:,2),'black');
% 
% hold on
% plot(freq_decim_shifted,phase_dif_lv(:,2),'red');


%---------------------------------------------------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------------------------------------------------

%---------------------------------------------------------------------------------------------------------------------------------
%---------------------------- Построение разностной фазы, сравнение с аналитическими данным --------------------------------------


% Делаем корректировку

Phase_dif_fft = phase_dif_fft_1kHz_step;
Phase_dif_lv = phase_dif_lv_1kHz_step; 

Phase_dif_fft = Phase_dif_fft - Phase_dif_fft(300,:);
Phase_dif_lv = Phase_dif_lv - Phase_dif_lv(300,:);

% Востановленная фаза почему-то получается с обратным знаком
Phase_dif_fft = -Phase_dif_fft;
Phase_dif_lv = -Phase_dif_lv;

dlmwrite('Phase_dif_fft.txt',Phase_dif_fft);
dlmwrite('Phase_dif_lv.txt',Phase_dif_lv);

% Сравниваем с аналитической фазой

%Phase_dif_anal = dlmread('Phase_dif_anal.txt');

Data_ID = fopen('../02_calculation_of_phase_variations/DPhi_1m.bin');
Phase_dif_anal = fread(Data_ID,'double');
fclose(Data_ID);
Phase_dif_anal=reshape(Phase_dif_anal, 601, 1001)';

time = zeros(601,1);

for n = 1:601
    time(n,1) = 0.1*(n-1);
end

%figure
%plot(freq_decim2_shifted,Phase_dif_anal(:,60))


figure;
imagesc(time,freq_decim2_shifted,Phase_dif_anal);
caxis([-pi; 2*pi]);
colorbar();



freq_decim2_shifted = zeros(1001,1);

for n = 1:1001
    freq_decim2_shifted(n,1) = n;
end

freq_decim2_shifted = freq_decim2_shifted - freq_decim2_shifted(501,1);



%Phase_dif_anal = Phase_dif_anal.';
%Phase_dif_lv = Phase_dif_lv;


% figure;
% plot(freq_decim2_shifted,-Phase_dif_lv(:,2),'black');
% 
% hold on
% plot(freq_decim2_shifted,Phase_dif_anal_2,'red');
% 
% figure;
% plot(freq_decim2_shifted,-Phase_dif_lv(:,3),'black');
% 
% hold on
% plot(freq_decim2_shifted,Phase_dif_anal_3,'red');

% figure;
% imagesc(time,freq_decim2_shifted,-Phase_dif_lv);
% caxis([-pi, 2*pi]);
% colorbar();
% 
% figure;
% imagesc(time,freq_decim2_shifted,Phase_dif_anal);
% caxis([-pi, 2*pi]);
% colorbar();

figure;
imagesc(time,freq_decim2_shifted,abs(Phase_dif_anal-Phase_dif_fft));
caxis([0, 0.9]);
colorbar();

figure;
imagesc(time,freq_decim2_shifted,Phase_dif_anal);
%caxis([-2*pi, 2*pi]);
colorbar();

figure;
imagesc(time,freq_decim2_shifted,Phase_dif_lv);
%caxis([-2*pi, 2*pi]);
colorbar();


figure;
imagesc(time,freq_decim2_shifted,abs(Phase_dif_anal-Phase_dif_lv));
%imagesc(time,freq_decim2_shifted,Phase_dif_anal-Phase_dif_lv);
%caxis([0, 0.9]);
colorbar();

std_fft = std(std(Phase_dif_anal-Phase_dif_fft))
mean_fft = mean(mean(Phase_dif_anal-Phase_dif_fft))

std_lv = std(std(Phase_dif_anal-Phase_dif_lv))
mean_lv = mean(mean(Phase_dif_anal-Phase_dif_lv))

% figure;
% mesh(time,freq_decim2_shifted,Phase_dif_lv);





