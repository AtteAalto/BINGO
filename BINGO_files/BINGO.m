function [Plink,chain,xstore,state,stats] = BINGO(data,state,parameters)
stats.start_time = clock;

if match_check(data.Tsam,data.ts)
    error('Sampling times are not consistent with the time series data!')
end

%Set the state of the sampler
q = state.q;
gamma = state.gamma;
r = state.r;
xs_old = state.xs;
Pold = state.P;
betsold = state.bets;
Jold = state.J;
Sold = state.S;
psiold = state.psi;
ma = state.ma;
mb = state.mb;

%Heuristic parameter to speed up topology sampling (default = 1)
if ~isfield(parameters,'Theur')
    parameters.Theur = 1;
end

%Reform data
Tsam = data.Tsam;
ts = data.ts;
u = [];
if isfield(data,'input')
    if match_check(data.input,data.ts)
        error('Input size is not consistent with the time series data!')
    end
    input = data.input;
    u = input{1}(:,1:end-1);
end

%TS data is put in one matrix y. The matrix Ser keeps track of the indices
%corresponding to different experiments.
y = ts{1};
Ser = zeros(4,size(ts,2));
Ser(1,1) = 1;
Ser(2,1) = size(y,2);
for j = 2:size(ts,2)
    y = [y, ts{j}];
    Ser(1:2,j) = [Ser(2,j-1)+1; Ser(2,j-1)+size(ts{j},2)];
    if isfield(data,'input')
        u = [u,input{j}(:,1:end-1)];
    end
end

%Check dimensions etc.
n = size(y,1);
n_in = size(u,1);
nstep = (size(xs_old,2)-size(Ser,2))/(size(y,2)-size(Ser,2));
M = size(psiold,2);
rany = [[min(y')',max(y')'];[min(u')',max(u')']];

%Check which variables are accounted for in each time series (excluding
%knocked-out genes in respective experiments)
Incl = ones(n,size(ts,2));
if isfield(data,'ko')
    for jser = 1:size(ts,2)
        Incl(data.ko{jser},jser) = zeros(length(data.ko{jser}),1);
    end
end 

%Check genes that are not knocked out in at least one experiment (normally
%geneList = 1:n)
geneList = find(sum(Incl,2)>0)';
excludedGenes = setdiff(1:n,geneList)

%Prior probability for the existence of a link (default = 1/n)
if isfield(parameters,'link_pr')
    log_link_pr = log(parameters.link_pr);
else
    log_link_pr = -log(n)*ones(n,n+n_in);
end
if numel(log_link_pr) == 1
    log_link_pr = log_link_pr*ones(n,n+n_in);
end

%Process the prior network information
S_aux = zeros(n,n+n_in);
if isfield(data,'sure')
    if size(data.sure,2) > n+.5
        S_aux = data.sure;
    else
        S_aux(1:n,1:n) = data.sure;
    end
end
Sold = max(Sold,S_aux);
Sold = min(Sold,1+S_aux);
S_aux = ones(n,n+n_in) - abs(S_aux);
S_aux(:,excludedGenes) = 0;
S_aux(excludedGenes,:) = 0;
Sold(:,excludedGenes) = 0;
Sold(excludedGenes,:) = 0;


%Rows 3-4 of Ser show the indices of different experiments in the finer grid.
Ser(3:4,1)  =[1; nstep*(Ser(2,1)-1)+1];
for jser = 2:size(Ser,2)
    Ser(3:4,jser) = [Ser(4,jser-1)+1; Ser(4,jser-1)+nstep*(Ser(2,jser)-Ser(1,jser))+1];
end

%Concatenate inputs as well
for j=1:size(u,2)
    u=[u(:,1:(j-1)*nstep),u(:,(j-1)*nstep+1)*ones(1,nstep),u(:,(j-1)*nstep+2:end)];  
end

%Variables miss/nomiss shows which measurements are/aren't missing
nomiss = ones(size(y));
if isfield(data,'missing')
    if norm(size(data.missing)-size(data.ts))>.1
        error('Missing measurements indices are not consistent with the time series data!')
    end   
    for i = 1:n
        for jser = 1:size(Ser,2)
            if size(data.missing{jser},1) > .5
                miss = find(abs(data.missing{jser}(:,1)-i)<.5);
                nomiss(i,data.missing{jser}(miss,2)+Ser(1,jser)-1) = zeros(1,length(miss));
            end
        end
    end
end

%Auxiliary indices showing which points in concatenated trajectory are used
%for probability distribution calculation and which points are compared
%with data
derind_full = (1:nstep*(Ser(2,1)-Ser(1,1)));
yind = 1 + (0:Ser(2,1)-1)*nstep;
d_full = (Tsam{1}(2:end)-Tsam{1}(1:end-1))/nstep;
for jser = 2:size(Ser,2)
    derind_full = [derind_full, (1:nstep*(Ser(2,jser)-Ser(1,jser)))+1+derind_full(end)];
    yind = [yind, yind(end)+1+(0:Ser(2,jser)-Ser(1,jser))*nstep];
    d_full = [d_full, (Tsam{jser}(2:end)-Tsam{jser}(1:end-1))/nstep];
end
d_full = reshape(ones(nstep,1)*d_full,1,[]);

%Calculate signal (quadratic) variation
nry = 0;
totvar = 0;
Total_time = 0;
for l = 1:size(Ser,2)
    nry = nry + sum(((y(:,Ser(1,l)+1:Ser(2,l))-y(:,Ser(1,l):Ser(2,l)-1)).^2./(Tsam{l}(2:end)-Tsam{l}(1:end-1))),2); 
    totvar = totvar + sum(abs(y(:,Ser(1,l)+1:Ser(2,l))-y(:,Ser(1,l):Ser(2,l)-1)),2);
    Total_time = Total_time + Incl(:,l).*(Tsam{l}(end)-Tsam{l}(1));
end 
nry = nry./Total_time;
totvar = totvar./Total_time;

%Embedding of the measurements to a piecewise linear function
mm = max(Ser(4,:)-Ser(3,:))+1;
Pr = zeros(mm,max(Ser(2,:)-Ser(1,:))+1);
Pr(1:nstep,1) = flipud((1:nstep)'/nstep);
Pr(mm-nstep+1:mm,end) = (1:nstep)'/nstep;
for j = 2:max(Ser(2,:)-Ser(1,:))
    Pr((j-2)*nstep+2:(j-1)*nstep+1,j) = (1:nstep)'/nstep;
    Pr((j-1)*nstep+1:j*nstep,j) = flipud((1:nstep)'/nstep);
end
Pintc = sin((1:nstep-1)'*(1:nstep-1)/nstep*pi)./(pi*ones(nstep-1,1)*(1:nstep-1))*2^.5;

%Store results and statistics
Plink = zeros(size(Sold));
chain = 0;
acctraj = 0;
xstore = 0*xs_old;
acctop = zeros(n,1);
acchyp = zeros(n,1);
accr = zeros(n,1);
yold = xs_old(:,yind);

%% Iterations
tic; time_mark = toc;
for k = 1:parameters.its
    %% Topology sampling
    for i = geneList
        S = Sold(i,:);
        %Check whether topology for row i is sampled on this time step (or 
        %only hyperparameters)
        top_change = (rand>.333);
        
        %Choose from two different types of moves
        inds = find(S_aux(i,:)>.5);
        topc = (rand>.5)*(sum(S(inds))>.5)*(sum(S(inds))<sum(S_aux(i,:))-.5);
       
        %Move type 1: Change one entry in S
        if top_change*(1-topc)
            indc = randi(sum(S_aux(i,:)),1);
            S(inds(indc)) = 1 - S(inds(indc));
        end  
        
        %Move type 2: Change one zero into one and one one into zero
        if top_change*topc
            ind1 = find(S(inds)>.5);
            ind0 = find(S(inds)<.5);
            indc01 = ind0(randi(length(ind0)));
            indc10 = ind1(randi(length(ind1)));
            S(inds(indc01)) = 1;
            S(inds(indc10)) = 0;
        end
        
        %Sample relevance parameters
        bets = (1-(parameters.ebeta)^2)^.5*betsold(i,:)+(parameters.ebeta)*randn(1,n+n_in);
        beta = .5+.45*bets;
        p_bets = exp(-abs(beta))./exp(-(beta-.5).^2/2/.45^2);
        beta = abs(beta);
        
        %Sample other hyperparameters
        gamma_tr = gamma(i)+parameters.egamma*nry(i)*randn;
        gamma_tr = 1e-4+abs(gamma_tr-1e-4);
        matr = ma(i)+parameters.ea*randn;
        matr = 1e-7+abs(matr-1e-7);
        mbtr = mb(i)+parameters.eb*randn;
        mbtr = 1e-7+abs(mbtr-1e-7);
    
        %Exclude the knockout data of gene i (if any)
        derind = derind_full;
        d = d_full;
        if sum(Incl(i,:)) < size(Ser,2) - .5
            derind = [];
            d = [];
            for jser = find(Incl(i,:)>.5)
                derind = [derind,derind_full((Ser(3,jser):Ser(4,jser)-1)-jser+1)];
                d = [d, d_full((Ser(3,jser):Ser(4,jser)-1)-jser+1)];
            end
        end   
        N = length(derind);
        
        %Form covariance matrices
        KM = zeros(M,M);
        KNM = zeros(N,M);
        for j = find(S(1:n)>.5)
            KM = KM + beta(j)*bsxfun(@minus,psiold(j,:),psiold(j,:)').^2;
            KNM = KNM + beta(j)*bsxfun(@minus,psiold(j,:),xs_old(j,derind)').^2;
        end
        for j = n + find(S(n+1:n+n_in)>.5)
            KM = KM + beta(j)*bsxfun(@minus,psiold(j,:),psiold(j,:)').^2;
            KNM = KNM + beta(j)*bsxfun(@minus,psiold(j,:),u(j-n,:)').^2;
        end
        KM = gamma_tr*exp(-KM);
        KNM = gamma_tr*exp(-KNM);
        
        %Compute the load
        %[k i]
        %size(d)
        %size(KNM)
        KC = chol(KM+1/q(i)*((KNM'*(d'.*KNM))) + 1e-5*eye(M));
        der = ((xs_old(i,derind+1)'-xs_old(i,derind)') - d'.*(mbtr-matr*xs_old(i,derind)'))/q(i);
        ld  =(KC'\(KNM'*der));
        
        
        %Part of the Wiener measure that is not in the CN sampler
        nrY = 0;
        for l = 1:size(Ser,2)   
            nrY = nrY + sum(((yold(i,Ser(1,l)+1:Ser(2,l))-yold(i,Ser(1,l):Ser(2,l)-1)).^2./(Tsam{l}(2:end)-Tsam{l}(1:end-1))),2);     
        end 
        
        %Cost function value
        J1 = .5*nrY/q(i) - .5*ld'*ld+log(prod(diag(KC))) - .5*log(det(KM+1e-5*eye(M))) - (mbtr-matr*xs_old(i,derind)')'*(xs_old(i,derind+1)'-xs_old(i,derind)')/q(i) + .5/q(i)*norm(d'.^.5.*(mbtr-matr*xs_old(i,derind)'))^2;
        PS = sum(S.*log_link_pr(i,:)) + log(prod(p_bets));
       
        %Acceptance of row i
        P_aux_ab = exp(.1*(ma(i)-matr+2*mb(i)-2*mbtr)/totvar(i));
        P_aux_gamma = gamma_tr/nry(i)*(30-gamma_tr/nry(i))/(gamma(i)/nry(i)*(30-gamma(i)/nry(i)))*exp(.2/nry(i)*(gamma(i)-gamma_tr));
        if P_aux_ab*P_aux_gamma*exp((PS-Pold(i)+Jold(i)-J1)/parameters.Theur) > rand
            Sold(i,:) = S;
            Pold(i) = PS;
            Jold(i) = J1;
            betsold(i,:) = beta;
            gamma(i) = gamma_tr;
            ma(i) = matr;
            mb(i) = mbtr;
            acctop(i) = acctop(i) + top_change;
            acchyp(i) = acchyp(i) + 1;
        end
        
        %Sampling of measurement noise variance r(i)
        rtr = r(i) + parameters.er*randn;
        rtr = 1e-8 + abs(rtr-1e-8);
        if (r(i)/rtr)^(1+sum(nomiss(i,:))/2)*exp(.00001./r(i)-.00001./rtr+sum((y(i,find(nomiss(i,:)>.5))-yold(i,find(nomiss(i,:)>.5))).^2)/2*(1/r(i)-1/rtr)) > rand
            r(i) = rtr;
            accr(i) = accr(i) + 1;
        end
    end
        
    
    %% Trajectory sampling
    
    %Sample the process noise variance
    qtr = q + parameters.eq*randn(size(q));
    qtr = .5e-5 + abs(qtr-.5e-5);
    
    %Sample the trajectory
    for l = 1:size(Ser,2)
        if isfield(data,'missing')  %Missing measurements
            if size(data.missing{l},1) > .5
                for i = 1:n
                    Csam = missing_data_sampler(data.missing{l},Tsam{l},r(i),qtr(i),i);
                    coef = (1-parameters.etraj^2)^.5*nomiss(i,Ser(1,l):Ser(2,l)) + (qtr(i)/q(i))^.5*(1-nomiss(i,Ser(1,l):Ser(2,l)));
                    yhat(i,Ser(1,l):Ser(2,l)) = y(i,Ser(1,l):Ser(2,l)) + coef.*(yold(i,Ser(1,l):Ser(2,l))-y(i,Ser(1,l):Ser(2,l))) + parameters.etraj*randn(1,Ser(2,l)-Ser(1,l)+1)*Csam';
                end
            else
                yhat(:,Ser(1,l):Ser(2,l)) = y(:,Ser(1,l):Ser(2,l)) + (1-parameters.etraj^2)^.5*(yold(:,Ser(1,l):Ser(2,l))-y(:,Ser(1,l):Ser(2,l))) + parameters.etraj*diag(r.^.5)*randn(n,Ser(2,l)-Ser(1,l)+1);   
            end
        else
            yhat(:,Ser(1,l):Ser(2,l)) = y(:,Ser(1,l):Ser(2,l)) + (1-parameters.etraj^2)^.5*(yold(:,Ser(1,l):Ser(2,l))-y(:,Ser(1,l):Ser(2,l))) + parameters.etraj*diag(r.^.5)*randn(n,Ser(2,l)-Ser(1,l)+1);   
        end        
        xs(:,Ser(3,l):Ser(4,l)) = sparse(diag((qtr./q).^.5))*(1-parameters.etraj^2)^.5*xs_old(:,Ser(3,l):Ser(4,l)) + (yhat(:,Ser(1,l):Ser(2,l))-sparse(diag((qtr./q).^.5))*(1-parameters.etraj^2)^.5*yold(:,Ser(1,l):Ser(2,l)))*Pr(1:(Ser(4,l)-Ser(3,l)+1),1:(Ser(2,l)-Ser(1,l)+1))'; 
        xs(:,Ser(3,l)+1:Ser(4,l)) = xs(:,Ser(3,l)+1:Ser(4,l)) + parameters.etraj*sparse(diag(qtr.^.5))*((nstep*d_full((Ser(3,l):Ser(4,l)-1)-l+1)).^.5.*reshape([Pintc*randn(nstep-1,n*(Ser(2,l)-Ser(1,l)));zeros(1,n*(Ser(2,l)-Ser(1,l)))],[],n)'); 
    end
    
    % Sample the pseudoinputs, and "mirror" them to the box containing the data. 
    % This corresponds to uniform prior in the box.
    psin = psiold + .025*randn(size(psiold));
    psin = min(psin,2*rany(:,2)-psin);
    psin = max(psin,2*rany(:,1)-psin);
    
    %Initialise the cost function value
    J1 = zeros(n,1);
    
    for i = geneList
        %Exclude the knockout data of gene i (if any)
        derind = derind_full;
        d = d_full;
        if sum(Incl(i,:)) < size(Ser,2) - .5
            derind = [];
            d = [];
            for jser = find(Incl(i,:)>.5)
                derind = [derind,derind_full((Ser(3,jser):Ser(4,jser)-1)-jser+1)];
                d = [d,d_full((Ser(3,jser):Ser(4,jser)-1)-jser+1)];
            end
        end   
        N = length(derind);
        
        %Compute the covariances
        KM = zeros(M,M);  %pseudoinputs
        KNM = zeros(N,M); %pseudoinput vs. trajectory
        for j = find(Sold(i,1:n)>.5)
            KM = KM + betsold(i,j)*bsxfun(@minus,psin(j,:),psin(j,:)').^2;
            KNM = KNM + betsold(i,j)*bsxfun(@minus,psin(j,:),xs(j,derind)').^2;
        end
        for j = n+find(Sold(i,n+1:n+n_in)>.5)
            KM = KM + betsold(i,j)*bsxfun(@minus,psin(j,:),psin(j,:)').^2;
            KNM = KNM + betsold(i,j)*bsxfun(@minus,psin(j,:),u(j-n,:)').^2;
        end
        KM = gamma(i)*exp(-KM);
        KNM = gamma(i)*exp(-KNM);
        
        %Compute the load 
        KC = chol(KM+1/qtr(i)*((KNM'*(d'.*KNM)))+1e-5*eye(M));
        der = ((xs(i,derind+1)'-xs(i,derind)') - (d'.*(mb(i)-ma(i)*xs(i,derind)')))/qtr(i);
        ld = (KC'\(KNM'*der));          

        %Part of the Wiener measure that is not in the CN sampler
        nrY = 0;
        for l = 1:size(Ser,2)
            nrY = nrY + sum(((yhat(i,Ser(1,l)+1:Ser(2,l))-yhat(i,Ser(1,l):Ser(2,l)-1)).^2./(Tsam{l}(2:end)-Tsam{l}(1:end-1))),2);     
        end 
        
        %Cost function value
        J1(i) = .5*nrY/qtr(i) - .5*ld'*ld+log(prod(diag(KC))) - .5*log(det(KM+1e-5*eye(M))) - (mb(i)-ma(i)*xs(i,derind)')'*(xs(i,derind+1)'-xs(i,derind)')/qtr(i) + .5/qtr(i)*norm(d'.^.5.*(mb(i)-ma(i)*xs(i,derind)'))^2;    
    end
       
    %Accept or reject the new sample
    P_aux_q = exp(sum(.00001./q-.00001./qtr+(log(q)-log(qtr)).*(1.001+.5*(size(y,2)-size(Ser,2)))));
    if P_aux_q*exp(sum(Jold-J1)) > rand
        Jold = J1;
        q = qtr;
        acctraj = acctraj + 1;
        xs_old = xs; 
        yold = yhat;
        psiold = psin;
    end
    
    %Notify the user if the simulation is estimated to take a long time (>15 min).
    if k == 100
        time_now = toc;
        left = (time_now-time_mark)*(parameters.its-k)/100;
        if left > 900
            h_left = floor(left/3600);
            left = left - 3600*h_left;
            min_left = floor(left/60);
            disp(['NOTE! Estimated time remaining for completion: ' num2str(h_left) ' h ' num2str(min_left) ' min'])
        end       
    end
    
    %Take every 10th sample to the actual distribution (thinning)
    if mod(k,10) < .5
        chain = chain + 1;
        Plink = Plink + Sold;
        xstore = xstore + xs_old;
        
        %Report progression
        if mod(k,10000) < .5
            time_now = toc;
            left = (time_now-time_mark)*(parameters.its-k)/10000;
            h_left = floor(left/3600);
            left = left - 3600*h_left;
            min_left = floor(left/60);
            left = floor(left-60*min_left);
            disp(['Iteration: ' num2str(k) ', Estimated time remaining: ' num2str(h_left) ' h ' num2str(min_left) ' min ' num2str(left) ' sec'])
            time_mark = time_now;
        end
    end
end

xstore = xstore/(parameters.its/10);
stats.acctraj = acctraj;
stats.acctop = acctop;
stats.acchyp = acchyp;
stats.accr = accr;
stats.fin_time = clock;

state.q = q;
state.gamma = gamma;
state.r = r;
state.xs = xs_old;
state.P = Pold;
state.bets = betsold;
state.J = Jold;
state.S = Sold;
state.psi = psiold;
state.ma = ma;
state.mb = mb;

