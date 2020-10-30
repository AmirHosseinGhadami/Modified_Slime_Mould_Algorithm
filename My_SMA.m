function [Destination_fitness,bestPositions,Convergence_curve]=My_SMA(N,Max_iter,lb,ub,dim,fobj)
disp('SMA is now tackling your problem')

% initialize position
bestPositions=zeros(1,dim);
Destination_fitness=inf;%change this to -inf for maximization problems
AllFitness = inf*ones(N,1);%record the fitness of all slime mold
weight = ones(N,dim);%fitness weight of each slime mold
%Initialize the set of random solutions
X=initialization(N,dim,ub,lb);
Convergence_curve=zeros(1,Max_iter);
it=1;  %Number of iterations
lb=ones(1,dim).*lb; % lower boundary 
ub=ones(1,dim).*ub; % upper boundary
z=0.03; % parameter

%-----------------------------------changes---------------------------------------------------
%init the archive 
archive.vb = [];
archive.positions= [];
archive.fitness = [];
%---------------------------------------------------------------------------------------------

% Main loop
while  it <= Max_iter
    
    

    %sort the fitness
    for i=1:N
        % Check if solutions go outside the search space and bring them back
        Flag4ub=X(i,:)>ub;
        Flag4lb=X(i,:)<lb;
        X(i,:)=(X(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        AllFitness(i) = fobj(X(i,:));
    end
   

    
    
    [SmellOrder,SmellIndex] = sort(AllFitness);  %Eq.(2.6)
    worstFitness = SmellOrder(N);
    bestFitness  = SmellOrder(1);

    S=bestFitness - worstFitness+eps;  % plus eps to avoid denominator zero

    %calculate the fitness weight of each slime mold
    for i=1:N
        for j=1:dim
            if i<=(N/2)  %Eq.(2.5)
                weight(SmellIndex(i),j) = 1+rand()*log10((bestFitness-SmellOrder(i))/(S)+1);
            else
                weight(SmellIndex(i),j) = 1-rand()*log10((bestFitness-SmellOrder(i))/(S)+1);
            end
        end
    end
     %update the best fitness value and best position
    if bestFitness < Destination_fitness
        bestPositions=X(SmellIndex(1),:);
        Destination_fitness = bestFitness;
    end
  %-----------------------------------changes_start------------------------------
  %archive.positions = [archive.positions;X];
  %add Allfitness to the archive
  archive.fitness = [archive.fitness ; AllFitness];
  
  
  
  
   

    if  ~isempty(archive.vb)
        %archive.fitness = archive.fitness(1:20,:); 
        %best_index = best_index(1:20,:);
        %archive.vb = archive.vb(best_index,:);
      
        %best_vb  = archive.vb(best_index,:);
        %twenty_good_vb = best_vb(1:50,:) ;
        
        
        best_vb_not_zero_index =  archive.vb~=0;
        
        best_vb_not_zero =archive.vb(best_vb_not_zero_index);
            
    end


  
  
    using_vb = zeros(N,dim);
    to_add_vb=zeros(N,dim);
    for i=1:N
        for j=1:dim
            if ~isempty(archive.vb)&& ~isempty(best_vb_not_zero)
                random = randi([1,length(best_vb_not_zero)]);
                using_vb(i,j) = best_vb_not_zero(random);

            else
                a = atanh(-(it/Max_iter)+1);
                using_vb(i,j) = unifrnd(-a,a);
            
            end
        end
    end

  %-----------------------------------changes_end------------------------------- 
  %------------------------------------------------------------------------------
    b = 1-it/Max_iter;
    for i=1:N 
        if rand<z     %Eq.(2.7)
            X(i,:) = (ub-lb)*rand+lb;
        else
            p  = tanh(abs(AllFitness(i)-Destination_fitness));  %Eq.(2.2)
            %-------------------------------changes_start-----------------------------
            %vb = unifrnd(-a,a,1,dim);  %Eq.(2.3) %deleted
            
            %--------------------------------changes_end------------------------------
            %-------------------------------------------------------------------------
            vc = unifrnd(-b,b,1,dim);
            for j=1:dim
                r = rand();
                A = randi([1,N]);  % two positions randomly selected from population
                B = randi([1,N]);
                
                if r<p    %Eq.(2.1)
                    %---------------------------changes_start-------------------------------------------
                    %instead of vb(j) ,vb(i,j)
                    X(i,j) = bestPositions(j)+ using_vb(i,j)*(weight(i,j)*X(A,j)-X(B,j));
                    
                    to_add_vb(i,j)= using_vb(i,j);
                    %---------------------------changes_end--------------------------------------------
                    %----------------------------------------------------------------------------------
                else
                    X(i,j) = vc(j)*X(i,j);
                end
            end
        end
    end
    %-------------------------------------------changes_start-----------------------------------------
    [archive.fitness,best_index] = sort(archive.fitness);
    archive.vb =[archive.vb ; to_add_vb];
    archive.fitness = archive.fitness(1:30,:); 
    best_index = best_index(1:30,:);
    archive.vb = archive.vb(best_index,:);
    
    
    %-------------------------------------------changes_end-------------------------------------------
    Convergence_curve(it)=Destination_fitness;
    it=it+1;
end

end
