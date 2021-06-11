%dot_q_py = load("dotq.txt") ; 
V = load("V_inertial_record.txt") ; 
%%V_prev = load("V.txt") ; 
q_py = load("q_part_record.txt") ; 
t = 0:0.001:9.999 ; 

%figure(1) ; plot(omegab(:,1), dot_q_sim(:,2), omegab(:,1), dot_q_sim(:,3), omegab(:,1), dot_q_sim(:,4)) ; 
%hold on ; 
%figure(1) ; plot(omegab(:,1), dot_q_py(:, 11), 'g--', omegab(:,1), dot_q_py(:, 12), 'k--', omegab(:,1), dot_q_py(:, 13), 'c--'); 
%hold off ; 


figure(2) ; plot(omegab(:,1), q_part(:,2), omegab(:,1), q_part(:,3), omegab(:,1), q_part(:,4)) ;
hold on ; 
figure(2) ; plot(omegab(:,1), q_py(:,8),  'g--', omegab(:,1), q_py(:,9), 'k--', omegab(:,1), q_py(:,10), 'c--') ;
hold off; 

n = size(V_prev); n = n(1)/8 ; 

%V = zeros(n, 6, 8) ; 
% for i = 1:1:n
%     % j es el n�mero de cuerpo, e i*h es el tiempo
%     for j = 1:1:8
%     V(i, :, j) = V_prev((8*(i-1)+j) , :) ; 
%     end
% end

figure(4) ;  hold on; 
plot(omegab(:,1), V(:, 4), 'b') ; % Velocidad angular del primer cuerpo (base)
plot(omegab(:,1), V(:, 5), 'c') ; 
plot(omegab(:,1), V(:, 6), 'm') ; 
plot(omegab(:,1), omegab(:,2), 'g--',omegab(:,1), omegab(:,3), 'k--', omegab(:,1), omegab(:,4), '--' ,'Color', [0, 0, 0] );
hold off ; 

% 
% figure(4) ; hold on ; 
% for i = 4:1:6
%     plot(omegab(:,1), V(:, i, 7)) ; % Velocidad angular del septimo cuerpo 
% end
% plot(omegab(:,1), omega7(:,2), 'g--',omegab(:,1), omega7(:,3), 'k--', omegab(:,1), omega7(:,4),  'c--' );
% hold off ; 
% 
% figure(5) ; hold on ; 
% for i = 1:1:3
%     plot(omegab(:,1), V(:, i, 7)) ; % Velocidad angular del septimo cuerpo 
% end
% plot(omegab(:,1), v7(:,2), 'g--',omegab(:,1), v7(:,3), 'k--', omegab(:,1), v7(:,4),  'c--' );
% hold off ; 

