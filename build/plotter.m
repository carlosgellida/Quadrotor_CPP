A = load("q_part_record.txt"); % Contiene la solución con cuaternios, y suma simple
%B = load("q_part_record2.txt"); % Contiene la solución con Euler-angle y matrix exponencial
t = 0:0.001:9.999 ;
figure(1) ; 
hold on ; 
plot(t, A(:, 8), 'b') ;  
%plot(t, B(:, 8), 'g--') ;  
plot(q_part(:, 1), q_part(:,2), 'r') ; 

figure(2) ; 
hold on ; 
plot(t, A(:, 9), 'b') ; 
%plot(t, B(:, 9), 'g--') ; 
plot(q_part(:, 1), q_part(:, 3), 'r') ; 

figure(3) ; 
hold on ; 
plot(t, A(:, 10), 'b'); 
%plot(t, B(:, 10), 'g--') ; 
plot(q_part(:, 1), q_part(:, 4), 'r') ; 

