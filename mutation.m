% function [population,centroid]=mutation(population,charge,count,iter,change,range,alpha,beta,rho,phi)
function [population,centroid]=mutation(population,charge,count,iter,change,range,alpha,beta,gbest)
% % mass=massCalculation(rank);
% mass=charge+.00005;
% centroid=zeros(1,size(population,2));
% for i=1:size(population,2)
%     centroid(1,i)=(mass*population(:,i))/sum(mass);
% end

[~,s]=sort(charge);
for i=1:size(population,1)
    for j = 1:size(population,1)
       distance(1,j) = norm(population(i,:)-population(j,:));
    end
    distance(find(distance==0))=[];
    [Novel_val, Novel_IX] = sort(distance,'ascend');
    dist = sum(Novel_val(Novel_IX(1:5)));
%     k=randperm(size(population,1));
%     if k~=i
%     dist=norm(population(k,:)-population(i,:));
%     else
%      dist=norm(population(s(i),:)-population(i,:));  
%     end
    pm=1/(dist+1);
    pc=.5+.5*(tanh((change/alpha)-beta));
    %pm=min(pc,pm);
    pm = 0.8*pm + 0.2*pc ;
%     pm = rho*pm + phi*pc ;
    for j=1:size(population,2)
        if ((rand < pm) || (pc > .98))
            %fprintf('l');
            if abs(population(i,j)) <= 0.00001
                delta=.5*range*rand;%ranje(j)
            else
                delta=min(.5*range*((iter-count)/iter)^20,abs(population(i,j)));
            end
            %delta=max(delta,population(i,j)/10); %bad results
            %delta=.5*range*((1-(count/iter))^2);
            if rand <.5
                population(i,:)=population(i,:)-rand*delta;
            else
                population(i,:)=population(i,:)+rand*delta;
            end
        end
    end
    
%     end
end