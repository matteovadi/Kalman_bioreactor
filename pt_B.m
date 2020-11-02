
function [Xest_b, D_b, actMSE_b , estMSE_b] = pt_B(Kp, R, Sin, Tc, u, X, y, N)

% INIZIALIZZAZIONE 
x0 = [1, 53, 0.5];
x = x0;
xp = zeros(3,1);
P = [0.1 0 0;
     0 45 0;
     0 0 0.027];
Xest_b = zeros(3,N);
Xest_b(:,1) = x;
D_b = zeros(3,N);
D_b(:,1) = diag(P);
G = zeros(3);
G(3,3) = 1;           
sigmaquadrow = 1e-6;
Q = zeros(3);
Q(3,3) = sigmaquadrow;

actMSE_b = zeros(N,1);
estMSE_b = zeros(N,1);
actMSE_b(1) = norm(X(1,:)-Xest_b(1:2,1)');    % errore reale 
estMSE_b(1) = sqrt(trace(P(1:2,1:2)));            % errore predetto dal filtro
 
% RECURSION:
 for k = 1:N-1
     % Prediction
     xp(1) = x(1) - Tc * u(k) * x(1) + Tc * x(3) * x(1) * x(2);
     xp(2) = x(2) + Tc * u(k) * (Sin-x(2)) - Tc * Kp * x(3) * x(1) * x(2);
     xp(3) = x(3);
     F = [1 - Tc * u(k) + Tc * x(3) * x(2), Tc * x(3) * x(1), Tc * x(1) * x(2);
          - Tc * Kp * x(3) * x(2), 1 - Tc * u(k) - Tc * Kp * x(3) * x(1), - Tc * Kp * x(1) * x(2);
          0, 0, 1]; 
     Pp = F*P*F' + G*Q*G';
     
     % Correction
     H = [1 0 0];
     K = Pp * H'* inv(H * Pp * H'+ R);  
     xc = xp + K * (y(k+1) - xp(1));      
     Pc = Pp * (eye(3)- H' * K');
     x = xc;
     P = Pc;
     
     Xest_b(:,k+1) = x;
     D_b(:,k+1) = diag(P); 
     
     actMSE_b(k+1) = norm(X(k,:) - Xest_b(1:2,k)');   
     estMSE_b(k+1) = sqrt(trace(P(1:2,1:2)));  
  end
end