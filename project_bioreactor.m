%% project_bioreactor %% 
clear 
close all
clc
load('dati1_bioreactor.mat')

N = length(u);
t = 1:N;
[Xest_a, D_a, actMSE_a , estMSE_a] = pt_A(alpha, Kp, R, Sin, Tc, u, X, y, N);
[Xest_b, D_b, actMSE_b , estMSE_b] = pt_B(Kp, R, Sin, Tc, u, X, y, N);

%% Punto A: Comparazione errori di stima delle variabili di stato rispetto al relativo intervallo di confidenza
figure(1)
for i=1:2
    subplot(2,2,i),plot(t,X(:,i)-Xest_a(i,:)','r',t,3*sqrt(D_a(i,:)),'b--',t,-3*sqrt(D_a(i,:)),'b--'), grid on
    title(['Punto C - caso A (variabile di stato x' int2str(i) ')'])
    legend(['Stima errore variabile di stato x(' int2str(i) ')' ],[ 'Intervallo di confidenza'])
    ylabel('Stima errore'), xlabel('Tempo') , xlim([0 N]), 
    if i == 1 
        ylim([-0.8 0.8])
    else
        ylim([-100 100])
    end
    grid on 
end
% Punto B: Comparazione errori di stima delle varibili di stato rispetto al relativo intervallo di confidenza (varibili x1 e x2)
for k=1:2
    subplot(2,2,2+k),plot(t,X(:,k)-Xest_b(k,:)','r',t,3*sqrt(D_b(k,:)),'b--',t,-3*sqrt(D_b(k,:)),'b--'), grid on
    title(['Punto C - caso B (variabile di stato x' int2str(k) ')'])
    legend(['Stima errore variabile di stato x(' int2str(k) ')'],[ 'Intervallo di confidenza'])
    ylabel('Stima errore'), xlabel('Tempo') , xlim([0 N]), 
    if k == 1 
        ylim([-0.8 0.8])
    else
        ylim([-100 100])
    end
    grid on
end
% variabile x3 
figure(2)
plot(t,alpha-Xest_b(3,:)','r',t,3*sqrt(D_b(3,:)),'b--',t,-3*sqrt(D_b(3,:)),'b--'), grid on
title('Punto C - caso B (variabile di stato x3)')
legend(['Stima di errore della variabile di stato x(3)' ],[ 'Intervallo di confidenza'])
ylabel('Stima errore'), xlabel('Tempo') , xlim([0 N]), grid on 

figure(7)
plot(t,alpha,'r',t,Xest_b(3,:), 'g')

%% Plot delle osservazioni y(t):
figure(3)
subplot(3,1,1), plot(t,y(:), 'LineWidth', 3)
title('Osservazioni di y(t)'), xlim([0 N]), ylim([0 60])                    %
% Plot delle stime delle varibili di stato x1 e x2 nel punto A e delle variabili di stato x1, x2 e x3 nel punto B
subplot(3,1,2),plot(t, Xest_a(1,t), 'g',t, Xest_a(2,t), 'r')
title('Stima delle variabili di stato (punto A)')
legend('Stima andamento di x1' , 'Stima andamento di x2')
xlim([0 N])
subplot(3,1,3),plot(t, Xest_b(1,t), 'b',t, Xest_b(2,t), 'm', t, Xest_b(3,t), 'k')
title('Stima delle variabili di stato (punto B)')
legend('Stima andamento di x1' , 'Stima andamento di x2', 'Stima andamento di x3')
xlim([0 N])

%% Differenza tra valore reale e valore stimato degli stati x1 e x2 nel punto A
figure(4)
for i=1:2
    subplot(2,2,i),plot(t,X(:,i),'g',t, Xest_a(i,:)','r')
    title(['Valori reali e stimati della variabile x(' int2str(i) ') - punto A']) 
    legend(['Valore reale di x(' int2str(i) ')'],['Valore stimato di x(' int2str(i) ')'])
    xlabel('Tempo'), xlim([0 N])
end
% Differenza tra valore reale e valore stimato degli stati x1 e x2 nel punto B
for k=1:2
    subplot(2,2,2+k),plot(t,X(:,k),'g',t, Xest_b(k,:)','r')
    title(['Valori reali e stimati della variabile x(' int2str(k) ') - punto B']) 
    legend(['Valore reale di x(' int2str(k) ')'],['Valore stimato di x(' int2str(k) ')'])
    xlabel('Tempo'), xlim([0 N])
end
 
%% Comparazione tra errore reale (verde) ed errore predetto dal filtro (rosso) nel punto A
figure(5)
subplot(2,1,1), plot(t,actMSE_a,'g',t,estMSE_a,'r')
title('Comparazione tra errore reale ed errore predetto dal filtro (punto A)')
legend('Errore reale' , 'Errore predetto dal filtro'), xlim([0 N])
% Comparazione tra errore reale (verde) ed errore predetto dal filtro (rosso) nel punto B
subplot(2,1,2), plot(t,actMSE_b,'g',t,estMSE_b,'r')
title('Comparazione tra errore reale ed errore predetto dal filtro (punto B)') 
legend('Errore reale' , 'Errore predetto dal filtro'), xlim([0 N])

%% Andamento P_i_i(k|k) nel punto A e nel punto B
figure(6)
subplot(2,1,1), plot(t,D_a(1,t),'g',t, D_a(2,t),'r')
title('Andamento di P_i_i(k|k) nel punto A') 
legend('P_1_1(k|k)','P_2_2(k|k) '),  xlim([0 N])
subplot(2,1,2),plot(t, D_b(1,t),'b',t, D_b(2,t),'m',t, D_b(3,t),'k')
title('Andamento di P_i_i(k|k) nel punto B') 
legend( 'P_1_1(k|k)' , 'P_2_2(k|k)' , 'P_3_3(k|k)'),  xlim([0 N])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%