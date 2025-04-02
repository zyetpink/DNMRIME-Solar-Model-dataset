%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Open Source Code for the Paper: 
% "RIME Optimization with Dynamic Multi-dimensional Random Mechanism and Nelder-Mead Simplex for Photovoltaic Parameter"
% This code is released under an open-source license for academic and research purposes.
% Please refer to the paper for detailed methodology and experiments.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  

%%
% fobj = @YourCostFunction    
% dim = number of your variables   
% MaxFEs = Maximum numbef of fitness evaluations
% lb=[lb1,lb2,...,lbn] where lbn is the lower bound of variable n  
% ub=[ub1,ub2,...,ubn] where ubn is the upper bound of variable n  
% If all the variables have equal lower bound you can just
% define lb and ub as two single number numbers

%%
clear all 
clc
rng('default')
Function_name='F18'; % Name of the test function that can be from F1 to F23



% Load details of the selected benchmark function
[lb,ub,dim,fobj]=Get_Functions_details(Function_name);
MaxFEs=1000*dim;
N = 50; %popultaion size

[Best_pos,Convergence_curve] = DNMRIME(N,MaxFEs,lb,ub,dim,fobj);
figure('Position',[269   240   660   290])
%Draw search space
subplot(1,2,1);
func_plot(Function_name);
title('Parameter space')
xlabel('x_1');
ylabel('x_2');
zlabel([Function_name,'( x_1 , x_2 )'])

%Draw objective space
subplot(1,2,2);
plot(Convergence_curve,'Color','r')
title('Objective space')
xlabel('Iteration');
ylabel('Best score obtained so far');

axis tight
grid on
box on
legend('DNMRIME')

display(['The best solution obtained by ,DNMRIME is : ', num2str(Best_pos)]);


