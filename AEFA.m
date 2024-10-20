
function [Fbest,Lbest,BestValues,MeanValues,Div,ExpExt,FirstP_D1,t1]=AEFA(func_num,N,max_it,FCheck,tag,Rpower,D)
%V:   Velocity.
%a:   Acceleration.
%Q:   Charge  
%D:   Dimension of the test function.
%N:   Number of charged particles.
%X:   Position of particles.
%R:   Distance between charged particle in search space.
%lb:  lower bound of the variables
%ub:  upper bound of the variables
%Rnorm: Euclidean Norm 
data=[];Div=[];ExpExt=[];
Rnorm=2; 
% Dimension and lower and upper bounds of the variables
% [lb,ub,D]=benchmark_range(func_num);  
lb=-100;ub=100;
%------------------------------------------------------------------------------------
%random initialization of charge population.
%X=initialization(D,N,ub,lb); 
X=rand(N,D).*(ub-lb)+lb;
% save('initial_100D.csv','X','-append','-ascii')
% X=load('initial_100D.csv');
%create the best so far chart and average fitnesses chart.
BestValues=[];MeanValues=[];
U=zeros;
I=zeros;
V=zeros(N,D);
change=0;
vel=zeros(N,D);
%-------------------------------------------------------------------------------------
tic;m=1;
for iteration=1:max_it
%Evaluation of fitness values of charged particles. 
 for i=1:N 
    %calculation of objective function for charged particle 'i'
%     fitness(i)=benchmark(X(i,:),func_num,D);
fitness(i)=cec17_func(X(i,:)',func_num);
 end
 asd=size(fitness);
 %fitness =benchmark(X,func_num,D);
 
    if tag==1
    [best, best_X]=min(fitness); %minimization.
    else
    [best, best_X]=max(fitness); %maximization.
    end        
    if iteration==1
       Fbest=best;Lbest=X(best_X,:);
    end
    if tag==1
      if best<Fbest  %minimization.
       Fbest=best;Lbest=X(best_X,:);
      end
    else 
      if best>Fbest  %maximization
       Fbest=best;Lbest=X(best_X,:);
      end
    end
    
BestValues=[BestValues Fbest];
MeanValues=[MeanValues mean(fitness)];

%-----------------------------------------------------------------------------------
% Charge 
Fmax=max(fitness); Fmin=min(fitness); Fmean=mean(fitness); 

if Fmax==Fmin
   M=ones(N,1);
   Q=ones(N,1);
else
    
   if tag==1 %for minimization
      best=Fmin;worst=Fmax; 
      
   else %for maximization
       
      best=Fmax;worst=Fmin; 
   end
  
   Q=exp((fitness-worst)./(best-worst)); 

end

Q=Q./sum(Q); 
%----------------------------------------------------------------------------------

fper=2; %In the last iteration, only 2-6 percent of charges apply force to the others.

%----------------------------------------------------------------------------------
 %%%%total electric force calculation
 if FCheck==1
     cbest=fper+(1-iteration/max_it)*(100-fper); 
     cbest=round(N*cbest/100);
 else
     cbest=N; 
 end 
    [Qs s]=sort(Q,'descend');
 for i=1:N
     E(i,:)=zeros(1,D);
     for ii=1:cbest
         j=s(ii);
         if j~=i
            R=norm(X(i,:)-X(j,:),Rnorm); %Euclidian distanse.
         for k=1:D 
             E(i,k)=E(i,k)+ rand*(Q(j))*((X(j,k)-X(i,k))/(R^Rpower+eps));
              
         end
         end
     end
 end
%---------------------------------------------------------------------------------- 
%Calculation of Coulomb constant
%----------------------------------------------------------------------------------
alpha=30;
K0=500;
   K2(iteration)=K0*exp(-alpha*iteration/max_it);
%---------------------------------------------------------------------------------- 
%%%Calculation of accelaration.
a=K2(iteration)*E; 
%---------------------------------------------------------------------------------- 
%----------------------------------------------------------------------------------
%Charge movement
%----------------------------------------------------------------------------------
V=rand(N,D).*V+a;
X=X+V;
X=max(X,lb);X=min(X,ub);   % Check the bounds of the variables
FirstP_D1(iteration) = X(1);
%----------------------------------------------------------------------------------
%plot charged particles 
%mask it if you do not need to plot them
%----------------------------------------------------------------------------------
% swarm(1:N,1,1)=X(:,1);
% swarm(1:N,1,2)=X(:,2);
% clf    
%     plot(swarm(:, 1, 1), swarm(:, 1, 2), 'X')  % drawing swarm movements
%     hold on;
%     plot(swarm(best_X,1,1),swarm(best_X,1,2),'*r') % drawning of best charged particle
%     axis([lb ub lb ub]);
%      title(['\fontsize{12}\bf Iteration:',num2str(iteration)]);
% pause(.2)
%---------------------------------------------------------------------------------

%% Diversity Calculation
% Div1=[];
% for i=1:N
%     B=abs((sum(norm(mean(X(i,:))-X(i,:)))^(1/2))/N);
%    Div1=[Div1,B];
% end
%  Div=[Div,sum(Div1)];
%  %% Exploration and Exploitation analysis
% for j=1:D
%     for i=1:N
%     A=abs(sum((sum((mean(X(:,j))-X(i,j))))/N)/D);
%     end
%     
% end
% ExpExt=[ExpExt,A];
ExpExt=zeros;Div=zeros;

end
% plot(data)
t1=toc;
end
function D = Diversity(x, ava_Xd)       
% This program is used to compute the diversity measure value of the population
% Reference£ºOlorunda O , Engelbrecht A P . Measuring Exploration/Exploitation in Particle Swarms using Swarm Diversity[C], 2008.

[N,D] = size(x);             
ava_Xd = repmat(ava_Xd,N,1);    
Distance = sum([x - ava_Xd].^2,2).^(1/2);
D = mean(Distance);
end