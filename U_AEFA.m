function [Fbest,Lbest,BestValues,MeanValues]=U_AEFA(func_num,N,max_it,FCheck,tag,Rpower,D)
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
Rnorm=2; A=[];
lb=-100;ub=100;
X=rand(N,D).*(ub-lb)+lb;
%create the best so far chart and average fitnesses chart.
BestValues=[];MeanValues=[];
V=zeros(N,D);
change=0;
%-------------------------------------------------------------------------------------
for iteration=1:max_it
    %Evaluation of fitness values of charged particles.
    for i=1:N
        fitness(i)=cec17_func(X(i,:)',func_num);
    end
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
    %% DVR
    if iteration>=2
        A=[A;Lbest-g_opt];
        m=10;
        if size(A,1)>m
            A=A(m:end,:);
        end
        qq=randperm(size(A,1));
        if rand<0.5
            kbest=Lbest+rand*A(qq(1),:);
        else
            kbest=Lbest-rand*A(qq(1),:);
        end
        fitness_diff=cec17_func(kbest',func_num);
        if fitness_diff<=best
            best=fitness_diff;
            Lbest=kbest;
            X(best_X,:)=kbest;
        end
    end
    BestValues=[BestValues Fbest];
    MeanValues=[MeanValues mean(fitness)];
    %-----------------------------------------------------------------------------------
    %% Charge
    Fmax=max(fitness);Fmin=min(fitness); Fmean=mean(fitness);
    if Fmax==Fmin
        M=ones(N,1);
        Q=ones(N,1);
    else
        if tag==1 %for minimization
            best1=Fmin;worst=Fmax;
        else %for maximization
            best1=Fmax;worst=Fmin;
        end
        Q=exp((fitness-worst)./(best1-worst));
    end
    Q=Q./sum(Q);
    %----------------------------------------------------------------------------------
    fper=2; %In the last iteration, only 2-6 percent of charges apply force to the others.
    %----------------------------------------------------------------------------------
    %% total electric force calculation
    if FCheck==1
        cbest=fper+(1-iteration/max_it)*(100-fper);
        cbest=round(N*cbest/100);
    else
        cbest=N;
    end
    [~, s]=sort(Q,'descend');
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
    %% Calculation of Coulomb constant
    K0=500;L1=300;beta=3;
    K(iteration)=K0./(1+exp(beta*(iteration-max_it/2)/L1));  % Hierarchical Coulomb's constant
    %----------------------------------------------------------------------------------
    %% Calculation of accelaration.
    a=K(iteration)*E;
    w1=1-iteration^3/max_it^3;
    w2=iteration^3/max_it^3;
    X=X(s,:);a=a(s,:);V=V(s,:);
    for ii=1:cbest
        j=s(ii);
        g_opt=Lbest;
        V(j,:)=rand.*V(j,:)+w1*a(j,:)+w2*(g_opt-X(j,:));
        X(j,:)=X(j,:)+V(j,:);
    end
    for i=cbest+1:N
        V(i,:)=rand.*V(i,:)+a(i,:);
        X(i,:)=X(i,:)+V(i,:);
    end
    %% Non-uniform mutation
    if (abs(Fbest-best)<.01*Fbest)% || abs(gbest_rank-rank(1))>gbest_rank)
        change=change+1;
        change=mod(change,50);
    else
        change=0;
    end
    [X]=mutation(X,Q,iteration,max_it,change,(ub-lb),4,5,Lbest);
    %% Check the bounds of the variables
    X=max(X,lb);X=min(X,ub);
    
end
end