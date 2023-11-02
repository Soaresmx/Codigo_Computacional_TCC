

%%
% Trabalho de TCC 02 - Aluno: Daniel Souza Soares.
% Orientador: Victor Frencl; Coorientador: Alexandre Tuoto
% Código para a Identificação de Sinais e Sistemas.

% Este código utiliza-se do banco de dados gerado pelo script "GetData.py".

% Fechando todos os arquivos abertos e limpando a memória do sistema

clear, close all; clc;

%%
% Carregando os arquivos do Excel, estes arquivos foram gerados durante
% ensaios com a planta via script "GetData.py".

% Carregando o arquivo referencia para a identificação de sistemas.
x_1 = load('output1.csv');

% Carregando arquivos auxiliares que serão utilizados para a validação da identificação de sistmas.
x_2 = load('output2.csv');
x_3 = load('output3.csv');
x_4 = load('output4.csv');

% Esta informação é util para coleta do intervalo de tempo total de
% amostras. Dessa forma, serão criados quatro vetores independentes para
% este procedimento. 

timeStamp_1 = x_1(2:end,1);
timeStamp_2 = x_2(2:end,1);
timeStamp_3 = x_3(2:end,1);
timeStamp_4 = x_4(2:end,1);

% Os vetores abaixo serão utilizados para calcular o tempo de Amostragem dos experimentos de modo
% individual, descartando a 1 linha. 

sample_times_1 = x_1(2:end,2);
sample_times_2 = x_2(2:end,2);
sample_times_3 = x_3(2:end,2);
sample_times_4 = x_4(2:end,2);

% Média do tempo de amostragem de cada ensaio, isto é cada sinal de entrada e saída amostrado. 
Ts_1 = mean(sample_times_1);
Ts_2 = mean(sample_times_2);
Ts_3 = mean(sample_times_3);
Ts_4 = mean(sample_times_4);

% Tempo médio de amostragem de todos os sinais amostrados. 

Ts = (Ts_1+Ts_2+Ts_3+Ts_4)/4;

% Criando os vetores de entrada e saída de cada arquivo separadamente

% Entrada e saída referente ao primeiro ensaio.   
input_1 = x_1(2:end,5);
output_1 = x_1(2:end,4);

% Entrada e saída referente ao segundo ensaio.
input_2 = x_2(2:end,5);
output_2 = x_2(2:end,4);

% Entrada e saída referente ao terceiro ensaio.
input_3 = x_3(2:end,5);
output_3 = x_3(2:end,4);

% Entrada e saída referente ao quarto ensaio.
input_4 = x_4(2:end,5);
output_4 = x_4(2:end,4);

%% Técnica para modelagem de um filtro digital tipo Média Móvel

% Definindo o tamanho da janela do filtro, este valor é de escolha do
% usuário.
no = 100; 

% Criação do filtro de Média Móvel para o primeiro ensaio.
filtro1 = ones(1, no) / no;
output_1 = conv(output_1, filtro1, 'same');

% Criação do filtro de Média Móvel para o segundo ensaio.
filtro2 = ones(1, no) / no;
output_2 = conv(output_2, filtro2, 'same');

% Criação do filtro de Média Móvel para o terceiro ensaio.
filtro3 = ones(1, no) / no;
output_3 = conv(output_3, filtro3, 'same');

% Criação do filtro de Média Móvel para o quarto ensaio.
filtro4 = ones(1, no) / no;
output_4 = conv(output_4, filtro4, 'same');

%% Nesta seção serão realizados a construção dos sinais, utilizando-se dos recursos disponível
% pela biblioteca do MATLAB. 

% Construção dos modelos que serão utilizados para a identificação de sistemas,
% utilizando os dados de saída, entrada e tempo de amostragem de cada experimento. 

data_y_val  = iddata(output_1,input_1,Ts_1); 
data_y_mod1 = iddata(output_2,input_2,Ts_2);
data_y_mod2 = iddata(output_3,input_3,Ts_3);
data_y_mod3 = iddata(output_4,input_4,Ts_4);

%%  Seção para a visualização gráfica de cada sinal. 

figure (2)

subplot(4,1,1)
hold on;
title(['Sinal para identificação c/ filtro média móvel utilizando-se uma janela de ' num2str(no) ' amostras']);
plot(timeStamp_1,input_1)
plot(timeStamp_1,output_1)
legend('Sin. entrada','Sin. saida');
hold off;

subplot(4,1,2)
hold on;
title (['Modelo aux 1 c/ filtro média móvel utilizando-se uma janela de ' num2str(no) ' amostras']);
plot(timeStamp_2,input_2)
plot(timeStamp_2,output_2)
legend('Sin. entrada','Sin. saida');
hold off;

subplot(4,1,3)
hold on;
title (['Modelo aux 2 c/ filtro média móvel utilizando-se uma janela de ' num2str(no) ' amostras']);
plot(timeStamp_3,input_3)
plot(timeStamp_3,output_3)
legend('Sin. entrada','Sin. saida');
hold off;

subplot(4,1,4)
hold on;
title (['Modelo aux 3 c/ filtro média móvel utilizando-se uma janela de ' num2str(no) ' amostras']);
plot(timeStamp_4,input_4)
plot(timeStamp_4,output_4)
legend('Sin. entrada','Sin. saida');
hold off;



%% Plot das Taxas de Amostragem: Ts, em formato de Histograma.

clc

figure (2)

title("Tempos de amostragens")
subplot(4,1,1);
histogram(sample_times_1);
title('Experimento 1');

subplot(4,1,2)
histogram(sample_times_2);
title('Experimento 2');

subplot(4,1,3)
histogram(sample_times_2);
title('Experimento 3');

subplot(4,1,4)
histogram(sample_times_2);
title('Experimento 4');

% Backup com os dados do experimento 
%save experimento 
%% Cálculo do modelo de identificação de sistemas. 

% Esta Seção adaptou o script utilizado no TCC do aluno: Lucas Raduy. 

% Para a obtenção de um modelo utilizando-se da técnica via FT do sinal é realizado o cálculo de 
% combinações de pólos e zeros. Com isso, o modelo escolhido será o que possuir o menor MSE.
% Nesse sentido, para este TCC, foram testados modelos com até 3 polos e 2
% zeros, pois aumentando a possibilidade do sistemas ser representado em
% polinômios de ordem elevada, há um aumento da probabilidade do sistema
% conter pólos e zeros imagináveis, o que compromete a estabilidade do
% sistema. 

% De acordo com a teoria de controle, o número máximo de zeros = qtd. polos - 1.

% O usuário deverá definir a quantidade máxima de polos que o modelo poderá
% possuir. 
numero_polos= 3;

numero_zeros= numero_polos-1;
modo= 'MSE';
np=(1:numero_polos)';
nz=(0:numero_zeros)';
k=1;

% Gerador da combinacao de polos e zeros

for i=1:numero_polos
    for j=0:i-1
        npnz(k,1)= i;
        npnz(k,2)= j;
        k=k+1;
    end
end

% Serve para o usuário acompanhar a evolução da contagem

disp(length(npnz))

% Método para o Cálculo dos modelos:

for i=1:length(npnz)

    % Calcula os modelos:
    modelo_y1_TF= tfest(data_y_val,npnz(i,1),npnz(i,2));
    modelo_y2_TF= tfest(data_y_mod1,npnz(i,1),npnz(i,2));
    modelo_y3_TF= tfest(data_y_mod2,npnz(i,1),npnz(i,2));
    modelo_y4_TF= tfest(data_y_mod3,npnz(i,1),npnz(i,2));

    % Testa os modelos:
    [y1_val, y1_NRMSE(i), y1_x0] = compare(data_y_val,modelo_y1_TF);
    [y2_val, y2_NRMSE(i), y2_x0] = compare(data_y_mod1,modelo_y2_TF);
    [y3_val, y3_NRMSE(i), y3_x0] = compare(data_y_mod2,modelo_y3_TF);
    [y4_val, y4_NRMSE(i), y4_x0] = compare(data_y_mod3,modelo_y4_TF);

    % Calcula o MSE, referente a cada modelo:
    fit_y1(i,1)=npnz(i,1);
    fit_y1(i,2)=npnz(i,2);
    fit_y1(i,3)=goodnessOfFit(y1_val.y, data_y_val.y,modo)

    fit_y2(i,1)=npnz(i,1);
    fit_y2(i,2)=npnz(i,2);
    fit_y2(i,3)=goodnessOfFit(y2_val.y,data_y_mod1.y,modo);

    fit_y3(i,1)=npnz(i,1);
    fit_y3(i,2)=npnz(i,2);
    fit_y3(i,3)=goodnessOfFit(y3_val.y,data_y_mod2.y,modo);

    fit_y4(i,1)=npnz(i,1);
    fit_y4(i,2)=npnz(i,2);
    fit_y4(i,3)=goodnessOfFit(y4_val.y,data_y_mod3.y,modo);

    % Núemro de iterações a serem realizadas.
    disp(i)
end

% Nesta etapa, os resultados obtidos serão salvos como uma variável a
% título de Backup. 
%save resultados
%% Varedura e escolha da melhor FT. 

% Esta Seção adaptou o script utilizado no TCC do aluno: Lucas Raduy. 

% Checagem da melhor combinação de polos e zeros da saída, experimento 1
[MSE,identf_1] = min(fit_y1(:,3));
bestNP_1=fit_y1(identf_1,1);
bestNZ_1=fit_y1(identf_1,2);

% bestNP_1=3;
% bestNZ_1=2;  %FPE: 98.63, MSE: 98.62 

%Exibe os melhores dados
disp(fit_y1(identf_1,:));

% Constroe a FT do sinal em análise

best_Y1= tfest(data_y_val,bestNP_1,bestNZ_1)

best_Y=tf(best_Y1.Numerator,best_Y1.Denominator);

best_Y.OutputName=['vazao s/ filtro'];
best_Y.InputName=['inversor'];
best_Y.Name='Modelo';

% Sintonia de controladores PID sugerida pelo MATLAB baseado na FT encontrada. 
step(best_Y,0:1000)
[C_sys_1, info_1] = pidtune(best_Y,'PID')

%% Técnica para diminuição da FT, o que deve melhorar a estabilidade do sistema.  

pid1 = reduce(best_Y,(bestNZ_1));
[num_pid1,den_pid1] = ss2tf(pid1.A,pid1.B,pid1.C,pid1.D);
pid1 = tf(num_pid1,den_pid1)
step(pid1,-10:0.5:1000)
[C_sys_11, info_11] = pidtune(pid1,'PI')

%% Checagem da melhor combinação de polos e zeros do exper 2

[MSE2,identf_2] = min(fit_y2(:,3));
bestNP_2=fit_y2(identf_2,1);
bestNZ_2=fit_y2(identf_2,2);

% bestNP_1=3; %
% bestNZ_1=2; % FPE: 0.7959, MSE: 0.7953  

%Exibe os melhores dados
disp(fit_y2(identf_2,:));

% Constroe a FT do sinal em análise

best_Y2= tfest(data_y_mod1,bestNP_2,bestNZ_2)

best_Y22=tf(best_Y2.Numerator,best_Y2.Denominator);


best_Y22.OutputName=['vazao'];
best_Y22.InputName=['inversor'];
best_Y22.Name='Modelo';

% Sintonia de controladores PID sugerida pelo MATLAB baseado na FT encontrada. 
step(best_Y22,0:500)
[C_sys_2, info_2] = pidtune(best_Y22,'PID')

%% Técnica para diminuição da FT, o que deve melhorar a estabilidade do sistema.  

pid2 = reduce(best_Y22,(bestNZ_2));
[num_pid2,den_pid2] = ss2tf(pid2.A,pid2.B,pid2.C,pid2.D);
pid2 = tf(num_pid2,den_pid2)
step(pid2,-10:0.5:500)
[C_sys_22, info_22] = pidtune(pid2,'PID')

%% Checagem da melhor combinação de polos e zeros do exper 3

[MSE3,identf_3] = min(fit_y3(:,3));
bestNP_3=fit_y3(identf_3,1);
bestNZ_3=fit_y3(identf_3,2);

% bestNP_1=3;
% bestNZ_1=2;  % FPE: 0.6508, MSE: 0.6504  

%Exibe os melhores dados
disp(fit_y3(identf_3,:));

% Constroe a FT do sinal em análise

best_Y3= tfest(data_y_mod2,bestNP_3,bestNZ_3)

best_Y33=tf(best_Y3.Numerator,best_Y3.Denominator);

best_Y33.OutputName=['vazao'];
best_Y33.InputName=['inversor'];
best_Y33.Name='Modelo';

% Sintonia de controladores PID sugerida pelo MATLAB baseado na FT encontrada. 
step(best_Y33,0:500)
[C_sys_3, info_3] = pidtune(best_Y33,'PID')

%% Técnica para diminuição da FT, o que deve melhorar a estabilidade do sistema.  

pid3 = reduce(best_Y33,(bestNZ_3));
[num_pid3,den_pid3] = ss2tf(pid3.A,pid3.B,pid3.C,pid3.D)
pid3 = tf(num_pid3,den_pid3)
step(pid3,-10:0.5:500)
[C_sys_33, info_33] = pidtune(pid3,'PID')

%% Checagem da melhor combinação de polos e zeros do exper 4

[MSE4,identf_4] = min(fit_y4(:,3));
bestNP_4=fit_y4(identf_4,1);
bestNZ_4=fit_y4(identf_4,2);

% bestNP_1=3;
% bestNZ_1=2;  % FPE: 2.221, MSE: 2.22

%Exibe os melhores dados
disp(fit_y4(identf_4,:));

% Constroe a FT do sinal em análise

best_Y4= tfest(data_y_mod3,bestNP_4,bestNZ_4)

best_Y44=tf(best_Y4.Numerator,best_Y4.Denominator);

%best_Y(1,2)=0;

best_Y44.OutputName=['vazao'];
best_Y44.InputName=['inversor'];
best_Y44.Name='Modelo';

% Sintonia de controladores PID sugerida pelo MATLAB baseado na FT encontrada. 
step(best_Y44,0:500)
[C_sys_4, info_4] = pidtune(best_Y44,'PID')

%% Técnica para diminuição da FT, o que deve melhorar a estabilidade do sistema.  

pid4 = reduce(best_Y44,(bestNZ_4));
[num_pid4,den_pid4] = ss2tf(pid4.A,pid4.B,pid4.C,pid4.D)
pid4 = tf(num_pid4,den_pid4)
step(pid4,-10:0.5:500)
[C_sys_44, info_44] = pidtune(pid4,'PID')

%save resultados_sintonias
%% Exibição de todos os parâmetros encontrados de sintonia de PID

%load ("resultados_sintonias.mat")
C_sys_1, %info_1
C_sys_2, %info_2
C_sys_3, %info_3
C_sys_4, %info_4

%% Exibição de todos os parâmetros encontrados de sintonia de PID c/ redução da FT.

%load ("resultados_sintonias.mat")
C_sys_11, %info_1
C_sys_22, %info_2
C_sys_33, %info_3
C_sys_44, %info_4
%save result_finais

%% Obtenção de Polos e Zeros de cada sinal 

%pzmap(best_Y)
%pzmap(C_sys_11)


%pzmap(best_Y22)
%pzmap(C_sys_22)
% 
% pzmap(best_Y33)
% pzmap(C_sys_33)
% 
% pzmap(best_Y44)
% pzmap(C_sys_44)
% 
% 


%% Compara os resultados, usando a curva do exp 1 como ref

figure (3)

subplot(4,1,1)
%hold on;
compare(data_y_val, best_Y,'b')
title('Comparação do modelo com o sinal de saída do processo.');
xlabel('Tempo');
%hold off;

subplot(4,1,2)
%hold on;
compare(data_y_mod1,best_Y,'b')
title('');
xlabel('Tempo');
%hold off;

subplot(4,1,3)
%hold on;
compare(data_y_mod2,best_Y,'b')
title('');
xlabel('Tempo');
%hold off;


subplot(4,1,4)
compare(data_y_mod3,best_Y,'b')
title('');
xlabel('Tempo');
%hold off;

%%

figure (4)

subplot(4,1,1)
%hold on;
compare(data_y_mod1,best_Y22,'b')
title('A identf. desta curva será a Ref. p/ os modls.');
xlabel('Tempo');
%hold off;

subplot(4,1,2)
%hold on;
compare(data_y_val,best_Y22,'b')
title('');
xlabel('Tempo');
%hold off;

subplot(4,1,3)
%hold on;
compare(data_y_mod2,best_Y22,'b')
title('');
xlabel('Tempo');
%hold off;

subplot(4,1,4)
compare(data_y_mod3,best_Y22,'b')
title('');
xlabel('Tempo');
%hold off;

%%
figure (5)

subplot(4,1,1)
%hold on;
compare(data_y_mod2,best_Y33,'b')
title('A identf. desta curva será a Ref. p/ os modls.');
xlabel('Tempo');
%hold off;

subplot(4,1,2)
%hold on;
compare(data_y_val,best_Y33,'b')
title('');
xlabel('Tempo');
%hold off;

subplot(4,1,3)
%hold on;
compare(data_y_mod3,best_Y33,'b')
title('');
xlabel('Tempo');
%hold off;

subplot(4,1,4)
compare(data_y_mod1,best_Y33,'b')
title('');
xlabel('Tempo');
%hold off;

%%

figure (6)

subplot(4,1,1)
%hold on;
compare(data_y_mod3,best_Y44,'b')
title('A identf. desta curva será a Ref. p/ os modls.');
xlabel('Tempo');
%hold off;

subplot(4,1,2)
%hold on;
compare(data_y_val,best_Y44,'b')
title('');
xlabel('Tempo');
%hold off;

subplot(4,1,3)
%hold on;
compare(data_y_mod2,best_Y44,'b')
title('');
xlabel('Tempo');
%hold off;

subplot(4,1,4)
compare(data_y_mod1,best_Y44,'b')
title('');
xlabel('Tempo');
%hold off;

%% Save_Indent

clear OUT_1

Ts_1; input_1 ; output_1;

OUT_1(:,1) = best_Y1.Numerator;
OUT_1((end+1):length(best_Y1.Denominator),1) = zeros;
OUT_1(:,2) = best_Y1.Denominator;

writematrix(OUT_1, "FT_planta.csv");

%%
clear OUT_2

Ts_1; input_1 ; output_1;

OUT_2(:,1) = round(pid1.Numerator{1,1},4);
OUT_2((end+1):length(pid1.Denominator),1) = zeros;
OUT_2(:,2) = round(pid1.Denominator{1,1},4);

writematrix(OUT_2, "FT_simplif_ident.csv");

%% Para o TCC 02 - Foi previsto que a sintonia de PID seria feita em MA via Resposta ao Degrau

% Váriaveis utilizadas
clear OUT;

Ts_1; input_1 ; output_1;

OUT(:,1) = best_Y1.Numerator;
OUT((end+1):length(best_Y1.Denominator),1) = zeros;
OUT(:,2) = best_Y1.Denominator;

writematrix(OUT, "FT_resp_Degrau.csv");


clear OUT2;

Ts_1; input_1 ; output_1;

OUT2(:,1) = best_Y1.Numerator;
OUT2((end+1):length(best_Y1.Denominator),1) = zeros;
OUT2(:,2) = best_Y1.Denominator;

writematrix(OUT2, "FT_ident_Planta.csv");



%% Identificação Utilizando do Critério PRBS

clear OUT

OUT(:,1) = input_1;
OUT(:,2) = output_1;
OUT(:,3) = Ts;

writematrix(OUT,"Sinal_PRBS.csv");

SIGNAL = load("Sinal_PRBS.csv");  % Carrega o vetor de entrada da probe 1

% Pré-definições do modelo

u = SIGNAL(:,1);
y = SIGNAL(:,2);

% Aplique o filtro ao sinal
%y = lowpass(y1,0.00001);

N = length(u);          % Leitura do tamanho do vetor de entrada

Ts = 0.01;

t = Ts:Ts:Ts*N;              % Vetor de tempo

% % % % % % ENTRADAS DO MODELO % % % % % % % % % %

nRI = 3;                 % Quantidade de Regressores na Entrada
nRO = 1;                 % Quantidade de Regressores na Saida
nRR = 5;                 % Quantidade de Regressores de Resíduo
p_max = -1;           % Número de Iterações Máximas
RMSE_min = 0.05;         % MSE mínimo aceitável

% % % % % % % Cálculo da Matriz PSI % % % % % %

% Maior regressor
if(nRI>nRO)
    if(nRI>nRR)
        M = nRI + 1;
    else
        M = nRR + 1;
    end
else
    if(nRO>nRR)
        M = nR0 + 1;
    else
        M = nRR + 1;
    end
end

% Cálculo do vetor Y
Y = y(M:(N));

% Condicoes Iniciais de PSI (primeiro regressor de saída)
PSI = y(M-1:N-1);

% Adiciona colunas no vetor conforme a quantidade de regressores na saida
for i=2:nRO
    PSI = [PSI y((M-i):(N-i))];
end

% Adiciona colunas no vetor conforme a quantidade de regressores na entrada
for i=1:nRI
    PSI = [PSI u((M-i):(N-i))];
end

% % % % % Cálculo dos Parâmetros e Saída Inicial % % % % % % % % % % % % %

% Cálculo dos parâmetros do modelo
O = pinv(PSI) * Y;

% Cálculo da Saída
YH = PSI * O;

% Cálculo do Resíduo
RE = [zeros(M-1, 1); (Y-YH)];

% Cálculo da Matriz do Resíduo
PSI_RE = RE((M-1):(N-1));
for i=2:nRR
    PSI_RE = [PSI_RE RE((M-i):(N-i))];
end

% Cálculado do MSE E RMSE
MSE = 0;
for i=1:(N-M+1)
    MSE = MSE +  (Y(i)-YH(i))*(Y(i)-YH(i));
end
RMSE = sqrt(MSE/N);

% Iterador
p=0;

while (RMSE > RMSE_min) && (p <= p_max)

    % Cálculo de PSI extendida
    PSI_EXT = [PSI PSI_RE];

    % Cálculo dos parâmetros do modelo
    O = pinv(PSI_EXT) * Y;

    % Cálculo da Saída
    YH = PSI_EXT * O;

    % Cálculo do Resíduo
    RE = [zeros(M-1, 1); (Y-YH)];

    % Cálculo da Matriz do Resíduo
    PSI_RE = RE((M-1):(N-1));
    for i=2:nRR
        PSI_RE = [PSI_RE RE((M-i):(N-i))];
    end

    % Cálculado do MSE E RMSE
    MSE = 0;
    for i=1:(N-M+1)
        MSE = MSE +  (Y(i)-YH(i))*(Y(i)-YH(i));
    end
    RMSE = sqrt(MSE/N);

    % Iteração
    p = p+1;

end

% Atribuindo as Condições Iniciais
YH = [y(1:M-1);YH];

fprintf("Menor RMSE encontrado: %.6f\nNumero de Iterações: %d\n\n", RMSE, p);

% % % % % % % % % % % % % GRAFICOS % % % % % % % % % % % % % % % % % % % %


% Primeiro Grafico (Entrada u(i))
subplot(2,2,1);
plot(t,u)                                           % Plot da tensão de entrada
title("Entrada");                                   % Titulo do grafico
xlabel("Amostras");                                 % Titulo eixo x
ylabel("Amplitude");                                % Titulo eixo y
grid on                                             % Grid habilitado

% Segundo Grafico (Saída y(i))
subplot(2,2,2);
plot(t,y,'r','LineWidth',1)                                             % Plot da saida real
title("Saida");                                     % Titulo do grafico
xlabel("Amostras");                                 % Titulo eixo x
ylabel("Amplitude");                                % Titulo eixo y
grid on                                             % Grid habilitado

% Terceiro Grafico (Saída yh(i))
subplot(2,2,3);
plot(t, YH, 'gr','LineWidth',1)                                            % Plot da saída do modelo
title("Saída Modelo ARMAX");                        % Titulo do grafico
xlabel("Amostras");                                 % Titulo eixo x
ylabel("Amplitude");                                % Titulo eixo y
grid on                                             % Grid habilitado

% Quarto Grafico (Sobreposição das Saídas)
subplot(2,2,4);
plot(t,y,'r','LineWidth',2);                       % Plot da saída real
hold on                                             % Habilita a Sobreposição
plot(t, YH, 'gr','LineWidth',0.5);                 % Plot da saída do modelo
title("Sobreposição das Saídas Real e do Modelo");  % Titulo do grafico
xlabel("Amostras");                                 % Titulo eixo x
ylabel("Amplitude");                                % Titulo eixo y
grid on                                             % Grid habilitado

DAT = iddata(YH,u,Ts);
SYS = tf(DAT.y',DAT.u');