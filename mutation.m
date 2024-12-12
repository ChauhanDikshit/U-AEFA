%% For non-uniform mutation please cite the refer paper: Kar, Devroop, et al. "Fuzzy mutation embedded hybrids of gravitational search and Particle Swarm Optimization methods for engineering design problems."
%% Engineering Applications of Artificial Intelligence 95 (2020): 103847.
function [population,centroid]=mutation(population,charge,count,iter,change,range,alpha,beta,gbest)
[~,s]=sort(charge);
for i=1:size(population,1)
    for j = 1:size(population,1)
       distance(1,j) = norm(population(i,:)-population(j,:));
    end
    distance(find(distance==0))=[];
    [Novel_val, Novel_IX] = sort(distance,'ascend');
    dist = sum(Novel_val(Novel_IX(1:5)));
    pm=1/(dist+1);
    pc=.5+.5*(tanh((change/alpha)-beta));
    for j=1:size(population,2)
        if ((rand < pm) || (pc > .98))

            if abs(population(i,j)) <= 0.00001
                delta=.5*range*rand;
            else
                delta=min(.5*range*((iter-count)/iter)^20,abs(population(i,j)));
            end
            if rand <.5
                population(i,:)=population(i,:)-rand*delta;
            else
                population(i,:)=population(i,:)+rand*delta;
            end
        end
    end
end
