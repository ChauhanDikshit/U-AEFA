% "U-AEFA: Online and offline learning-based unified artificial electric
% field algorithm for real parameter optimization." Knowledge-Based Systems (2024)
% Dikshit Chauhan, Anupam Trivedi, Anupam Yadav
% anupuam@gmail.com
clear all;
clc;
for i=1
    rng('default');
    rng(1);
    N=100;D=30;
    max_FE=10000*D;
    max_it=round(max_FE/N);
    FCheck=1; R=1;
    tag=1; % 1: minimization, 0: maximization
    rand ('state', sum(100*clock))
    func_num=i
    [Fbest_U_AEFA,Lbest1,BestValues1,MeanValues1]=U_AEFA(func_num,N,max_it,FCheck,tag,R,D);Fbest_U_AEFA
    [Fbest_aefa,Lbest_aefa,BestValues_aefa,MeanValues_aefa]=AEFA(func_num,N,max_it,FCheck,tag,R,D);Fbest_aefa
end


