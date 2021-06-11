dot_q_save = load("dotq.txt") ; 
q_py = load("q.txt") ; 
%V_0 = load("V_0.txt") ; 
V_prev = load("V.txt") ; 

figure(1) ; plot(omega0(:,1), q_sim(:,2), omega0(:,1), q_sim(:,3), omega0(:,1), q_sim(:,4)) ; 
hold on ; 
figure(1) ; plot(omega0(:,1), q_py(:, 1), 'g--' , omega0(:,1), q_py(:, 2), 'k--', omega0(:,1), q_py(:, 3), 'm--'); 






n = size(V_prev); n = n(1)/3 ; 

V = zeros(n, 6, 3) ; 
for i = 1:1:n
    % j es el número de cuerpo, e i*h es el tiempo
    for j = 1:1:3
    V(i, :, j) = V_prev((3*(i-1)+j) , :) ; 
    end
end

figure(3) ; hold off; hold on ; 
for i = 4:1:6
    plot(omega0(:,1), V(:, i, 3)) ; % Velocidad angular del quinto cuerpo 
end
