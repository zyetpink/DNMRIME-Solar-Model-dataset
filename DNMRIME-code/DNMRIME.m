%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Open Source Code for the Paper: 
% "RIME Optimization with Dynamic Multi-dimensional Random Mechanism and Nelder-Mead Simplex for Photovoltaic Parameter"
% This code is released under an open-source license for academic and research purposes.
% Please refer to the paper for detailed methodology and experiments.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [Best_rime,Convergence_curve] = DNMRIME(N, MaxFEs, lb, ub, dim, fobj)
disp('DNMRIME is now tackling your problem')
%%%%%%%%%%%%%%%%%%%%%%%----Initialization----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Best_rime = zeros(1, dim); % Initialize the best solution
Best_rime_rate = inf; % Initialize best fitness value as infinity (use -inf for maximization problems)
Rimepop = initialization(N, dim, ub, lb);% % Call the `initialization` function to initialize the population
Rime_rates = zeros(1, N); % Initialize the fitness value
newRime_rates = zeros(1, N);% Prepare an array to store updated fitness values of the population
Lb = lb .* ones(1, dim); % lower boundary
Ub = ub .* ones(1, dim); % upper boundary
FEs = 0;
it= 1;
Convergence_curve = []; % Store convergence data


W = 5; % Soft-rime parameters
NMsPro=0.1;% Probability of executing the Nelder-Mead simplex method

%维度跳跃机制的比例
initialJumpRate = 0.3;% The parameter of Lambda
decayRate = 0.05; % The parameter of Delta 
jumpMagnitude = 0.1; 





% Calculate the fitness value of the initial position
for i = 1:N
    Rime_rates(1, i) = fobj(Rimepop(i, :)); % Calculate the fitness value for each search agent
    FEs = FEs + 1;
    % Make greedy selections
    if Rime_rates(1, i) < Best_rime_rate
        Best_rime_rate = Rime_rates(1, i);
        Best_rime = Rimepop(i, :);
    end
end

% Main loop
while FEs < MaxFEs
    RimeFactor = (rand - 0.5) * 2 * cos((pi * FEs / (MaxFEs / 10))) * (1 - round(FEs * W / MaxFEs) / W); % Parameters of Eq.(17),(18),(19)
    E = sqrt(FEs / MaxFEs); % Eq.(20),the parameter of Epsilon
    newRimepop = Rimepop; % Recording new populations
    normalized_rime_rates = normr(Rime_rates);  % Normalize
    for i = 1:N
        for j = 1:dim
            % Soft-rime search strategy
            r1 = rand();
            if r1 < E
                newRimepop(i, j) = Best_rime(1, j) + RimeFactor * ((Ub(j) - Lb(j)) * rand + Lb(j)); % Eq.(3)
            end
            % Hard-rime puncture mechanism
            r2 = rand();
            if r2 < normalized_rime_rates(i)
                newRimepop(i, j) = Best_rime(1, j); % Eq.(7)
            end
        end
    end
    for i = 1:N
        % Boundary absorption
        Flag4ub = newRimepop(i, :) > ub;
        Flag4lb = newRimepop(i, :) < lb;
        newRimepop(i, :) = (newRimepop(i, :) .* (~(Flag4ub + Flag4lb))) + ub .* Flag4ub + lb .* Flag4lb;
        newRime_rates(1, i) = fobj(newRimepop(i, :));
        FEs = FEs + 1;
        % Positive greedy selection mechanism
        if newRime_rates(1, i) < Rime_rates(1, i)
            Rime_rates(1, i) = newRime_rates(1, i);
            Rimepop(i, :) = newRimepop(i, :);
            if newRime_rates(1, i) < Best_rime_rate
               Best_rime_rate = Rime_rates(1, i);
               Best_rime = Rimepop(i, :);
            end
        end
    end
     
     %%%%%%%%%%%%%%%----Dynamic multi-dimensional random mechanism----%%%%%%%%%%%%%%%%%
     jumpRate = initialJumpRate * exp(-decayRate * FEs); % Eq.(34),the parameter of Epsilon
     dimJump = round(dim * jumpRate); %Eq.(34),the parameter of k
     jumpIndices = randperm(dim, dimJump); 

     b=1-FEs/MaxFEs;% Eq.(31),the parameter of Eta

     for i = 1:N
        for j = jumpIndices% Dimension selection

            jumpSize =(Ub(j) - Lb(j)) * jumpMagnitude * (2 * rand - 1);
            newRimepop(i, j) = b*Levy(1)*(Rimepop(i, j) + jumpSize);        
            newRimepop(i, j) = max(min(newRimepop(i, j),    Ub(j)), Lb(j)); % Ensure the value is within the boundaries
        end
        r=rand;
        r1=(2*pi)*r;% Eq.(36)
        for j=1:dim % in j-th dimension 
            %newRimepop(i,j)= Rimepop(i,j)*abs(sin(r1)) - 0.75*sin(r1)*abs(-0.75*Best_rime(j)-0.75*Rimepop(i,j));
            newRimepop(i,j) = Rimepop(i,j) * abs(sin(r1)) - (9/16) * sin(r1) * abs(-Best_rime(j) - Rimepop(i,j));% Eq.(35)
        end
       
        newRime_rates(1, i) = fobj(newRimepop(i, :));
        FEs = FEs + 1;
        if newRime_rates(1, i) < Rime_rates(1, i)
            Rime_rates(1, i) = newRime_rates(1, i);
            Rimepop(i, :) = newRimepop(i, :);
        
            if newRime_rates(1, i) < Best_rime_rate
                   Best_rime_rate = newRime_rates(1, i);
                   Best_rime =  newRimepop(i, :);
            end
        end
      
     end %end DMRM

     %%%%%%%%%%%%%%%----Dynamic multi-dimensional random mechanism----%%%%%%%%%%%%%%%%%
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%----NMs-----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     for i=1:N
         r3=rand();
         if r3<NMsPro
            options = optimset('MaxFunEvals', floor(MaxFEs * 0.1));
            [x, fval, ~, output]  = fminsearchbnd(fobj,Rimepop(i, :),Lb,Ub,options);
                if fval < Rime_rates(1, i)
                    Rime_rates(1, i) = fval;
                    Rimepop(i, :) = x;
         
                    if fval < Best_rime_rate
                        Best_rime_rate = fval;
                        Best_rime = x;
                    end
                end
            FEs = FEs + output.funcCount;
         end
     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%----NMs-----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    Convergence_curve(it) = Best_rime_rate;
    it=it+1;
end
end

%Lévy flight 
function L=Levy(d)
beta=1.5;
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
u=randn(1,d)*sigma;
v=randn(1,d);
step=u./abs(v).^(1/beta);
L=0.01*step;
end


