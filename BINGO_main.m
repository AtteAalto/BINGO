%% === BINGO ===
addpath('./BINGO_files/')

%initialization
[data,state,parameters]=BINGO_init(data);

% MCMC Burn-in
[~,chain,~,state,stats]=BINGO(data,state,parameters);
disp_stats(' BURN-IN COMPLETE',stats,chain,parameters.its)

% Actual sampling
parameters.its=10000;
[Plink,chain,xstore,state,stats]=BINGO(data,state,parameters);
disp_stats(' SAMPLING COMPLETE',stats,chain,parameters.its)

%The end result:
confidence_matrix=Plink/chain;


%% Collect more samples

%Store old information
chain_old=chain;
Plink_old=Plink;
xstore_old=xstore;

%Run new iterations
parameters.its=10000;
[Plink,chain,xstore,state,stats]=BINGO(data,state,parameters);

%Combine old and new results
xstore=chain_old/(chain+chain_old)*xstore_old+chain/(chain+chain_old)*xstore;
Plink=Plink_old+Plink;
chain=chain_old+chain;
disp_stats(' SAMPLING COMPLETE',stats,chain,parameters.its)

%The end result:
confidence_matrix=Plink/chain;
