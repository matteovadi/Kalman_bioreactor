
function [Xest_a, D_a, actMSE_a , estMSE_a] = pt_A(alpha, Kp, R, Sin, Tc, u, X, y, N)

% INIZIALIZZAZIONE 
x0 = [1, 0];
x = x0;
xp = zeros(2,1);
P = [0.1 0; 
     0 1111.11];
Xest_a = zeros(2,N);
Xest_a(:,1) = x;
D_a = zeros(2,N);
D_a(:,1) = diag(P);

actMSE_a = zeros(N,1);
estMSE_a = zeros(N,1);
actMSE_a(1) = norm(X(1,:) - Xest_a(:,1)');    % errore reale 
estMSE_a(1) = sqrt(trace(P));            % errore predetto dal filtro
 
% RECURSION:
 for k = 1:N-1
     % Prediction
     xp(1) = x(1) - Tc * u(k) * x(1) + Tc * alpha(k) * x(1) * x(2);
     xp(2) = x(2) + Tc * u(k) * (Sin-x(2)) - Tc * Kp * alpha(k) * x(1) * x(2);
     F = [1 - Tc * u(k) + Tc * alpha(k) * x(2), Tc * alpha(k) * x(1);
          - Tc * Kp * alpha(k) * x(2), 1 - Tc * u(k) - Tc * Kp * alpha(k) * x(1)];
     Pp = F*P*F';
     
     % Correction
     H = [1 0];
     K = Pp * H'* inv(H * Pp * H'+ R);  
     xc = xp + K * (y(k+1) - xp(1));      
     Pc = Pp * (eye(2) - H' * K');
     x = xc;
     P = Pc;
     
     Xest_a(:,k+1) = x;
     D_a(:,k+1) = diag(P); 
     
     actMSE_a(k+1) = norm(X(k,:)- Xest_a(:,k)');   
     estMSE_a(k+1) = sqrt(trace(P));  
  end
end