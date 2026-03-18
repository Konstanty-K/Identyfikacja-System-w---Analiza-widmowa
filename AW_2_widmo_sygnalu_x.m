%% Sekcja 2: Widmo amplitudowe |X(k)|

Fs = 1/Tp;              
Xk = fft(x);            

f = (0:N/2-1)*(Fs/N);   

A = abs(Xk)/N;          
A = A(1:N/2);           
A(2:end) = 2*A(2:end);  

figure;
stem(f, A);
xlabel('Czestotliwosc [Hz]');
ylabel('Amplituda');
title('Widmo amplitudowe x');
grid on;


xlabel('Częstotliwość [Hz]');
ylabel('|X(k)|'); 
title('Widmo amplitudowe x');
grid on; 
grid minor; 
xlim([0 50]); grid on;