% Cooking a sausage link with various methods
% By David Rey
% ME3333 Spring 2020

clear all; close all; clc;

%% Sausage Parameters

%Note tA,tB,tC are the time to cook values for tasks A,B,C

Ti = 5+273;                     %initial temp in Kelvin
Tf = 75+273;                    %final temp in Kelvin

dia = 0.025;                    %sausage dia
k = 0.3;                        %thermal conductivity
den = 1050;                     %density
cap = 3200;                     %specific heat capacitance
alpha = k/(den*cap);            %thermal diffusivity

h = [200 22 22 22 95];          %heat transfer coefficients
Tsur = [100 175 200 230 80];    %surround temps in Kelvin
Tsur = Tsur + 273;

L = dia/2;                      %characteristic length, more conservative estimate

bi = [0 0 0 0 0];               %initializations
tA = [0 0 0 0 0];
tB = [0 0 0 0 0];
tC = [0 0 0 0 0];

Fo_B = [0 0 0 0 0];

Y = [0 0 0 0 0];
Guess = [0 0 0 0 0];
Diff = [0 0 0 0 0];

ReErrA = [0 0 0 0 0];
ReErrB = [0 0 0 0 0];

TC = [0 0 0 0 0];

%% Task A: assuming lumped capacitance; determine bi number and time to cook for each method

for m = 1:5                     %for loop that calculates Bi numbers & time to cook for each method 
bi(m) = h(m)*L/k;
tA(m) = -(den*L*cap*log((Tf-Tsur(m))/(Ti-Tsur(m))))/(2*h(m));
end

trA = round(tA,2);              %time to cook rounded for visiblity on bar graph

figure(1);
bar(trA,0.5,'facecolor',[0.9290, 0.6940, 0.1250]);
text(1:length(trA),tA,num2str(trA'),'vert','bottom','horiz','center'); 
box off
grid on
ylabel('Cooking Time (seconds)')
xlabel('Cooking Method')
ax = gca;
ax.XTick = [1 2 3 4 5]; 
ax.XTickLabels = {'Boiling \newline Bi = 8.33','Oven @ 175C \newline Bi = 0.91','Oven @ 200C \newline Bi = 0.91','Oven @ 230C \newline Bi = 0.91','Sous Vide \newline Bi = 3.95'};
ax.XTickLabelRotation = 45;
title ('Task A: Lumped Capacitance Assumption')


%% Task B: assuming as semi-infite solid; determine Fo number and time to cook for each method

for m = 1:5                             %for loop that perform transcendental eqn and calculates time to cook & Fo numbers for each method
    Y(m) = (Tf-Ti)/(Tsur(m)-Ti);        %desired output
    
    for ts = 10:1:78000                 %step throu time steps of 10s from 10s to 10,000s;for each step calc a guess to take the difference of the desired output - guess
        Guess(m) = erfc(L/(2*(alpha*ts)^0.5)) - exp((h(m)*L/k)+(h(m)^2*alpha*ts)/(k^2))*(erfc(L/(2*(alpha*ts)^0.5)+(h(m)*(alpha*ts)^0.5/k)));
        
        Diff(m) = Y(m) - Guess(m);
        
        if Diff <= 0.0001
            tB(m) = ts;                 %if the guess is within 0.0001 of accuracy then set the value of the guess as the time in takes to cook (tB)
            Fo_B(m) = (alpha*tB(m)/L^2);%use determined tB to get the Fo number
            
            break
        else
            tB(m) = ts;
            Fo_B(m) = (alpha*tB(m)/L^2);%if the guess is never within accuracy, take value at upper limit as the time in takes to cook (tB); for sous vide case
        end
    end
end

trB = round(tB,2);
figure(2);
bar(trB,0.5,'facecolor',[0.6350, 0.0780, 0.1840]);
text(1:length(trB),tB,num2str(trB'),'vert','bottom','horiz','center'); 
box off
grid on
ylabel('Cooking Time (seconds)')
xlabel('Cooking Method')
ax = gca;
ax.XTick = [1 2 3 4 5]; 
ax.XTickLabels = {'Boiling \newline Fo = 5.54','Oven @ 175C \newline Fo = 2.66','Oven @ 200C \newline Fo = 2.05','Oven @ 230C \newline Fo = 1.62','Sous Vide \newline Fo = 44.57'};
ax.XTickLabelRotation = 45;
title ('Task B: Semi-infinite Solid Assumption')

%% Task C: using finite difference scheme to determine time to cook for each method
N = 20;                             %number of nodes
Fo_C = 0.1;                         %for stability, want Fo < or = to 0.5; chose value this value for smaller time step

dr = (L)/(N-1);                     %spacial step
dt = Fo_C*dr^2/alpha;               %time step

for i = 1:5                         %for each of the 5 conditions...
    
    T_old = zeros(N,1);             %array that stores old temperature values, initialize arrays with zeros
    T_new = T_old;                  %array that stores new temperature values, initialize arrays with zeros           
    T_old(:,1) = Ti;                %set temperature values to initial value
    
    for p = 1:1000000               %step through 1,000,000 seconds
    
        for m = 2:N-1               %interior nodes
            T_new(m) = Fo_C*(1-dr/(2*(m-1)*dr))*T_old(m-1)+(1-2*Fo_C)*T_old(m)+Fo_C*(1+dr/(2*(m-1)*dr))*T_old(m+1);
        end
        
        T_new(1) = T_new(2);        %boundary node, center of sausage 
        T_new(N) = (h(i)*Tsur(i)+(k/dr)*T_new(N-1))/(h(i)+k/dr);   %boundary node, sausage surface
    
        if T_new(1) >= Tf
            tC(i) = (p-1)*dt;      %if center temp reach 75C or more, determine the time to reach and break out of the loop
            TC(i) = T_new(N);
            break
        end
        
        T_old = T_new;              %new values become the old
        
    end
end

trC = round(tC,2);
figure(3);
bar(trC,0.5,'facecolor',[0.3010, 0.7450, 0.9330]);
text(1:length(trC),trC,num2str(trC'),'vert','bottom','horiz','center'); 
box off
grid on
ylabel('Cooking Time (seconds)')
xlabel('Cooking Method')
ax = gca;
ax.XTick = [1 2 3 4 5]; 
ax.XTickLabels = {'Boiling','Oven @ 175C','Oven @ 200C','Oven @ 230C','Sous Vide'};
ax.XTickLabelRotation = 45;
title ('Task C: Using Finite-Difference Approach')

%% Relative Error Graphs
data = [tA/60;tB/60;tC/60];

for m = 1:5
ReErrA(m) = abs((tA(m)-tC(m))/tC(m))*100;   %get relative error
ReErrB(m) = abs((tB(m)-tC(m))/tC(m))*100;   %get relative error
end

ReErrA = round(ReErrA,2);
ReErrB = round(ReErrB,2);

figure(4);                      
bar(ReErrA,0.5,'facecolor',[0.9290, 0.6940, 0.1250]);

text(1:length(ReErrA),ReErrA,num2str(ReErrA'),'vert','bottom','horiz','center'); 
box off
grid off
ylabel('Relative Error (%)')
xlabel('Cooking Method')
ax = gca;
ax.XTick = [1 2 3 4 5]; 
ax.XTickLabels = {'Boiling \newline Bi = 8.33','Oven @ 175C \newline Bi = 0.91','Oven @ 200C \newline Bi = 0.91','Oven @ 230C \newline Bi = 0.91','Sous Vide \newline Bi = 3.95'};
ax.XTickLabelRotation = 45;
title ('Relative Error of Lumped Capacitance Approach \newline ')

figure(5);
bar(ReErrB,0.5,'facecolor',[0.6350, 0.0780, 0.1840]);
text(1:length(ReErrB),ReErrB,num2str(ReErrB'),'vert','bottom','horiz','center'); 
box off
grid off
ylabel('Relative Error (%)')
xlabel('Cooking Method')
ax = gca;
ax.XTick = [1 2 3 4 5]; 
ax.XTickLabels = {'Boiling \newline Fo = 5.54','Oven @ 175C \newline Fo = 2.66','Oven @ 200C \newline Fo = 2.05','Oven @ 230C \newline Fo = 1.62','Sous Vide \newline Fo = 44.57'};
ax.XTickLabelRotation = 45;
title ('Relative Error of Semi-Infinitely Long Solid Approach \newline ')

TC = TC - 273;
TC = round(TC,2);
figure (6)
bar(TC,0.5,'facecolor',[0.3010, 0.7450, 0.9330]);
text(1:length(TC),TC,num2str(TC'),'vert','bottom','horiz','center'); 
box off
grid off
ylabel('Surface Temperature (C)')
xlabel('Cooking Method')
ax = gca;
ax.XTick = [1 2 3 4 5]; 
ax.XTickLabels = {'Boiling','Oven @ 175C','Oven @ 200C','Oven @ 230C','Sous Vide'};
ax.XTickLabelRotation = 45;
title ('Sausage Surface Temperature When Fully Cooked')
