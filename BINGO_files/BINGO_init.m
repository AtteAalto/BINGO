function [data,state,parameters]=BINGO_init(data)

%% Set the parameters of the sampler.

%Step length of the Crank-Nicolson sampler (should be >.05)
parameters.etraj = .1;

%Other step length parameters:
parameters.egamma = .1;
parameters.ea = .005;
parameters.eb = .005;
parameters.ebeta = .125;
parameters.er = .0001;
parameters.eq = .0002;

%link_pr = p / (1-p) where p is the prior probability for the existence of
%a link. This controls the sparsity level of the network. Default is 1/n if
%it is not given.
%parameters.link_pr = 1/10;

%Heuristic temperature variable to speed up the topology sampling. Value 1
%is default option if it is not set, and it corresponds to exactly correct
%sampling.
parameters.Theur = 1;

%Number of iterations (in the burn-in)
parameters.its = 3000;

%The number of steps in the trajectory between measurements. A good default
%is 4, but the sampling gets faster if it is reduced.
nstep = 4;

%Number of pseudoinputs. Tradeoff between accuracy/speed.
nr_pi = 50;

%The cell array Tsam in the data-structure contains the sampling times of
%the measurements. If the user doesn't specify these, this code creates 
%them, assuming constant sampling frequency (1 sample / time unit). 
if isfield(data,'Tsam')
    if size(data.Tsam,2) < 1.5
        if size(data.Tsam{1},2) < 1.5
            Tsam = {};
            for j = 1:size(data.ts,2)
                Tsam = {Tsam{1:size(Tsam,2)}, data.Tsam{1}(1,1)*(0:size(data.ts{j},2)-1)};
            end
            data.Tsam = Tsam;
        end
    else
        for j = 1:size(data.ts,2)
            if size(data.Tsam{j},2) < 1.5
                data.Tsam{j} = data.Tsam{j}(1,1)*(0:size(data.ts{j},2)-1);
            end
        end
    end
else
    Tsam = {};
    for j = 1:size(data.ts,2)
        Tsam = {Tsam{1:size(Tsam,2)}, 0:size(data.ts{j},2)-1};
    end
    data.Tsam = Tsam;
end

%Check the dimension of the time series
n = size(data.ts{1},1);
for j = 2:size(data.ts,2)
    if size(data.ts{j},1) ~= n
        error('Time series dimensions do not agree!')
    end
end

%Perform linear interpolation if there is missing measurements. Note that
%the code also takes into account the missing measurements in the sampling.
if isfield(data,'missing')
    data.ts = missing_data_interpolation(data);
end


%Create times and indices for plotting the trajectory estimate.
indlast = 0;
for jser = 1:size(data.ts,2)
    data.plot_index{jser} = (1:nstep*(size(data.ts{jser},2)-1)+1) + indlast;
    indlast = nstep*(size(data.ts{jser},2)-1) + 1 + indlast;   
    d = (data.Tsam{jser}(2:end)-data.Tsam{jser}(1:end-1))/nstep; 
    data.fine_times{jser} = [0, cumsum(reshape(ones(nstep,1)*d, 1, []))];
end


%% IMPORTANT! 
% Scaling data. Here each dimension (gene) is scaled so that the difference 
% between the maximum and minimum value of each gene is one. This does not 
% make sense if the data contains genes that are not expressed properly. 
% Therefore the data should be carefully pre-processed, if this scaling is 
% used!

maxs = zeros(size(data.ts{1},1),1);
mins = 1e5*ones(size(data.ts{1},1),1);
for jser = 1:size(data.ts,2)
    maxs = max(maxs,max(data.ts{jser},[],2));
    mins = min(mins,min(data.ts{jser},[],2));
end

for jser = 1:size(data.ts,2)
    data.ts{jser} = data.ts{jser}./(maxs-mins);
end

maxs = maxs./(maxs-mins);
mins = maxs-1;

%% Initialize the state of the sampler. Do not touch this part!

n_in = 0;
if isfield(data,'input')
    n_in = size(data.input{1},1);
end

%Initial connectivity guess
state.S = rand(n,n+n_in)>.9;

%Initial state for hyperparameters beta
state.bets = abs(randn(n,n+n_in)*.5) + 1e-3;

%Initial pseudoinputs. 
state.psi = [(maxs-mins).*rand(n,nr_pi)+mins.*(ones(1,nr_pi));rand(n_in,nr_pi)];

%Initial trajectory x is obtained by linear interpolation of the time series data.
Ser = zeros(4,size(data.ts,2));
Ser(1,1) = 1;
Ser(2,1) = size(data.ts{1},2);
Ser(3:4,1) = [1; nstep*(Ser(2,1)-1)+1];
for jser = 2:size(data.ts,2)
    Ser(1:2,jser) = [Ser(2,jser-1)+1;Ser(2,jser-1)+size(data.ts{jser},2)];
    Ser(3:4,jser) = [Ser(4,jser-1)+1; Ser(4,jser-1)+nstep*(Ser(2,jser)-Ser(1,jser))+1];
end

clear('xs')
for jser = 1:size(Ser,2)
    for jj = 1:(Ser(2,jser)-Ser(1,jser))
        xs(:,Ser(3,jser)+(jj-1)*nstep:Ser(3,jser)+jj*nstep-1) = data.ts{jser}(:,jj)*(1-(0:nstep-1)/nstep)+data.ts{jser}(:,jj+1)*(0:nstep-1)/nstep;
    end
    xs(:,Ser(4,jser)) = data.ts{jser}(:,end);
end
state.xs = xs;


nry = 0;
nry_q = 0;
Ttot = 0;
for l = 1:size(Ser,2)  
    nry = nry + sum((data.ts{l}(:,2:end)-data.ts{l}(:,1:end-1)).^2./(data.Tsam{l}(2:end)-data.Tsam{l}(1:end-1)),2);
    nry_q = nry_q + sum((data.ts{l}(:,2:end)-data.ts{l}(:,1:end-1)).^2,2);
    Ttot = Ttot + data.Tsam{l}(end)-data.Tsam{l}(1);    
end 
nry = nry/Ttot;
nry_q = nry_q/Ttot;


%Initialise other hyperparameters
state.r = .0006*.1*ones(n,1);
state.gamma = nry;        
state.q = nry_q/20;
state.ma = .1*ones(n,1);
state.mb = .05*ones(n,1);

%Initialise cost function values
state.P = -1e8*ones(n,1);
state.J = 1e8*ones(n,1);


