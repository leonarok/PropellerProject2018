function [] = percentdone_print(num,num_fin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------PRINTS THE PERCENTAGE COMPLETED TOWARDS CONVERGENCE GOAL---------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

frac=num_fin/num;
clc
fprintf('%6.2f%% completed\n',frac*100)


end

