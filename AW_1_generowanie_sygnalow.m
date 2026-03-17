%% Sekcja 1: Generowanie sygnałów (Tp=0.001, N=2000)
rng(77,'twister'); % stałe ziarno generatora liczb losowych

Tp = 0.001;
N=2000
n = 0:N-1;
sigma = 0.8;
x = sin(2*pi*5*n*Tp) + 0.5*sin(2*pi*10*n*Tp) + 0.25*sin(2*pi*30*n*Tp);
e = sigma * randn(1,N);
H = tf([0.1], [1, -0.9], Tp);
v = lsim(H, e, n*Tp);

figure; 
subplot(3,1,1)
plot(n*Tp, [x(:)]);  % Transpozycja (:) na kolumny  Nx3

subplot(3,1,2)
plot(n*Tp, [e(:)]); 

subplot(3,1,3)
plot(n*Tp, [v(:)]); 
legend('x', 'e (biały)', 'v (kolorowy)'); 
title('Sygnały w dziedzinie czasu dyskretnego');
xlabel('Czas [s]');
ylabel('Amplituda');
grid on;
grid minor;