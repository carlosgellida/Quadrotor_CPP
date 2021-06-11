A = load("q_part_record.txt"); % Contiene la solución con cuaternios, y suma simple
B = load("q.txt"); % Contiene la solución utilizando el programa hecho en Python
t = 0:0.001:9.999 ;
t2 = 0:0.001:10 ; 
figure(1) ; 
hold on ; 
plot(t, A(:, 8), 'b') ;  
plot(t2, B(:, 11), 'g--') ;  
plot(q_part(:, 1), q_part(:,2), 'r') ; 

figure(2) ; 
hold on ; 
plot(t, A(:, 9), 'b') ; 
plot(t2, B(:, 12), 'g--') ; 
plot(q_part(:, 1), q_part(:, 3), 'r') ; 

figure(3) ; 
hold on ; 
plot(t, A(:, 10), 'b'); 
plot(t2, B(:, 13), 'g--') ; 
plot(q_part(:, 1), q_part(:, 4), 'r') ; 

