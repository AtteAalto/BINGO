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
parameters.nstep = nstep;

%Number of pseudoinputs. Tradeoff between accuracy/speed.
nr_pi = 50;

%The cell array Tsam in the data-structure contains the sampling times of
%the measurements. If the user doesn't specify these, this code creates 
%them, assuming constant sampling frequency (1 sample / time unit). 
if isfield(data,'Tsam')
    Tsam = data.Tsam;
    if size(data.Tsam,2) < 1.5
        if size(data.Tsam,1) < 1.5
            if size(data.Tsam{1},2) < 1.5
                Tsam = {};
                for j1 = 1:size(data.ts,1)
                    for j2 = 1:size(data.ts,2)
                        Tsam{j1,j2} = data.Tsam{1}(1,1)*(0:size(data.ts{j1,j2},2)-1);
                    end
                end
                data.Tsam = Tsam;
            end
        else
            for j1 = 1:size(data.ts,1)
                for j2 = 1:size(data.ts,2)
                    Tsam{j1,j2} = data.Tsam{j1,1}(1,1)*(0:size(data.ts{j1,j2},2)-1);
                end
            end
            data.Tsam = Tsam;
        end       
    else
        for j1 = 1:size(data.ts,2)
            for j2 = 1:size(data.ts,1)
                if size(data.Tsam{j2,j1},2) < 1.5
                    Tsam{j2,j1} = data.Tsam{j2,j1}(1,1)*(0:size(data.ts{j2,j1},2)-1);
                end
            end
        end
        data.Tsam = Tsam;
    end
else
    Tsam = {};
    for j1 = 1:size(data.ts,2)
        for j2 = 1:size(data.ts,1)
            Tsam{j2,j1} = 0:size(data.ts{j2,j1},2)-1;
        end
    end
    data.Tsam = Tsam;
end

%Check the dimension of the time series
for jd = 1:size(data.ts,1)
    n(jd) = size(data.ts{jd,1},1);
    for j = 2:size(data.ts,2)
        if size(data.ts{jd,j},1) ~= n(jd)
            error('Time series dimensions do not agree!')
        end
    end
end


if ~isfield(data,'pseudotime')
    data.pseudotime = false(size(data.ts,1),1);
end




%Perform linear interpolation if there is missing measurements. Note that
%the code also takes into account the missing measurements in the sampling.
if isfield(data,'missing')
    data.ts = missing_data_interpolation(data);
end


%Create times and indices for plotting the trajectory estimate.
indlast = 0;
for jd = 1:size(data.ts,1)
    for jser = 1:size(data.ts,2)
        data.plot_index{jd,jser} = (1:nstep*(size(data.ts{jd,jser},2)-1)+1) + indlast;
        indlast = nstep*(size(data.ts{jd,jser},2)-1) + 1 + indlast;   
        d = (data.Tsam{jd,jser}(2:end)-data.Tsam{jd,jser}(1:end-1))/nstep; 
        data.fine_times{jd,jser} = [0, cumsum(reshape(ones(nstep,1)*d, 1, []))];
    end
end


%% IMPORTANT! 
% Scaling data. Here each dimension (gene) is scaled so that the difference 
% between the maximum and minimum value of each gene is one. This does not 
% make sense if the data contains genes that are not expressed properly. 
% Therefore the data should be carefully pre-processed, if this scaling is 
% used!

allmaxs = [];
allmins = [];

for jd = 1:size(data.ts,1)

    maxs = zeros(size(data.ts{jd,1},1),1);
    mins = 1e5*ones(size(data.ts{jd,1},1),1);
    for jser = 1:size(data.ts,2)
        maxs = max(maxs,max(data.ts{jd,jser},[],2));
        mins = min(mins,min(data.ts{jd,jser},[],2));
    end

    for jser = 1:size(data.ts,2)
        data.ts{jd,jser} = data.ts{jd,jser}./(maxs-mins);
    end

    allmaxs = [allmaxs; maxs./(maxs-mins)];
    allmins = [allmins; maxs-1];
end
    
maxs = allmaxs;
mins = allmins;

%% Initialize the state of the sampler. Do not touch this part!

n_in = 0;
if isfield(data,'input')
    n_in = size(data.input{1},1);
end

%Initial connectivity guess
state.S = rand(sum(n),sum(n)+n_in)>.9;

%Initial state for hyperparameters beta
state.bets = abs(randn(sum(n),sum(n)+n_in)*.5) + 1e-3;

%Initial pseudoinputs. 
state.psi = [(maxs-mins).*rand(sum(n),nr_pi)+mins.*(ones(1,nr_pi));rand(n_in,nr_pi)];

%Initial trajectory x is obtained by linear interpolation of the time series data.

for jd = 1:size(data.ts,1)
    
    Ser = zeros(4,size(data.ts,2));
    Ser(1,1) = 1;
    Ser(2,1) = size(data.ts{jd,1},2);
    Ser(3:4,1) = [1; nstep*(Ser(2,1)-1)+1];
    for jser = 2:size(data.ts,2)
        Ser(1:2,jser) = [Ser(2,jser-1)+1;Ser(2,jser-1)+size(data.ts{jser},2)];
        Ser(3:4,jser) = [Ser(4,jser-1)+1; Ser(4,jser-1)+nstep*(Ser(2,jser)-Ser(1,jser))+1];
    end

    for jser = 1:size(Ser,2)
        for jj = 1:(Ser(2,jser)-Ser(1,jser))
            xs{jd}(:,Ser(3,jser)+(jj-1)*nstep:Ser(3,jser)+jj*nstep-1) = data.ts{jd,jser}(:,jj)*(1-(0:nstep-1)/nstep)+data.ts{jd,jser}(:,jj+1)*(0:nstep-1)/nstep;
        end
        xs{jd}(:,Ser(4,jser)) = data.ts{jd,jser}(:,end);
    end
end
state.xs = xs;


nry = zeros(sum(n),1);
nry_q = zeros(sum(n),1);
for jd = 1:size(data.ts,1)
    Ttot = 0;
    for l = 1:size(Ser,2)  
        nry(sum(n(1:jd-1))+1:sum(n(1:jd))) = nry(sum(n(1:jd-1))+1:sum(n(1:jd))) + sum((data.ts{jd,l}(:,2:end)-data.ts{jd,l}(:,1:end-1)).^2./(data.Tsam{jd,l}(2:end)-data.Tsam{jd,l}(1:end-1)),2);
        nry_q(sum(n(1:jd-1))+1:sum(n(1:jd))) = nry_q(sum(n(1:jd-1))+1:sum(n(1:jd))) + sum((data.ts{jd,l}(:,2:end)-data.ts{jd,l}(:,1:end-1)).^2,2);
        Ttot = Ttot + data.Tsam{jd,l}(end)-data.Tsam{jd,l}(1);    
    end 
    nry(sum(n(1:jd-1))+1:sum(n(1:jd))) = nry(sum(n(1:jd-1))+1:sum(n(1:jd)))/Ttot;
    nry_q(sum(n(1:jd-1))+1:sum(n(1:jd))) = nry_q(sum(n(1:jd-1))+1:sum(n(1:jd)))/Ttot;
end


%Initialise other hyperparameters
state.r = .0006*.1*ones(sum(n),1);
state.gamma = nry;        
state.q = nry_q/20;
state.ma = .1*ones(sum(n),1);
state.mb = .05*ones(sum(n),1);

%Initialise cost function values
state.P = -1e8*ones(sum(n),1);
state.J = 1e8*ones(sum(n),1);


