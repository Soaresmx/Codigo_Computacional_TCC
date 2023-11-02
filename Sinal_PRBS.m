%Autor: Lucas Raduy Gomes de Camargo
% Código adptado para a planta de vazão e pH. 

P=7;                            %Sequencia PRBS
L=(2^P)-1;                      %tamanho do vetor necessario para criar o PRBS
O=10;

x=randi([0,1],P,O);
x=[x;zeros(L-P,O)];             %Deixa o vetor do tamanho necessario para a gereação do sinal
xorOut=0;

for i=1:L
    x=circshift(x,1);           %Desloca o vetor
    for j=1:O
        xorOut=xor(x((7),j),x((6),j));    %Aplica a operação xor nos bits corretos para um PRBS15
        
        x(1,j)=xorOut;                %substitui a primeira posição do vetor com o valor calculado pela operação XOR
    end
end
prbs=x;
x_norm=(prbs).*(10);
figure
stairs(prbs(:,1))
ylim([-0.2 1.1])
title('Sinal PRBS')                                 
t = linspace(5,1.5e3,127)';                % Adptação no vetor tempo para garantir a estabilização do sinal da planta quando submetida a uma mudança de amplitude do sinal.
t = (t) + 5*rand(length(prbs),1);
x_duplo2=[t,prbs];                         % Adaptação no formato em que os dados serão salvos. 
csvwrite('sigControl.data',x_duplo2);

fprintf('Média: %f \n',mean(x));
fprintf('Variância: %f \n',var(x));