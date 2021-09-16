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
nrData = size(data.ts,1);
nstep = parameters.nstep;

%Heuristic parameter to speed up topology sampling (default = 1)
if ~isfield(parameters,'Theur')
    parameters.Theur = 1;
end

%Reform data
Tsam = data.Tsam;
ts = data.ts;
u = [];
rany = [];
if isfield(data,'input')
%     if match_check(data.input,data.ts)
%         error('Input size is not consistent with the time series data!')
%     end
    input = data.input;
    u = input{1}(:,1:end-1);
end

%TS data is put in one matrix y. The matrix Ser keeps track of the indices
%corresponding to different experiments.
dataInd = [];
for jd = 1:nrData
    y{jd} = ts{jd,1};
    Ser{jd} = zeros(4,size(ts,2));
    Ser{jd}(1,1) = 1;
    Ser{jd}(2,1) = size(y{jd},2);
    Ser{jd}(3:4,1) = [1; nstep*(Ser{jd}(2,1)-1)+1];
    for j = 2:size(ts,2)
        y{jd} = [y{jd}, ts{jd,j}];
        Ser{jd}(1:2,j) = [Ser{jd}(2,j-1)+1; Ser{jd}(2,j-1)+size(ts{jd,j},2)];
        Ser{jd}(3:4,j) = [Ser{jd}(4,j-1)+1; Ser{jd}(4,j-1)+nstep*(Ser{jd}(2,j)-Ser{jd}(1,j))+1];
        if jd == 1 && isfield(data,'input')
            u = [u,input{j}(:,1:end-1)];
        end
    end
    n(jd) = size(y{jd},1);
    rany = [rany; [min(y{jd}')',max(y{jd}')']];
    dataInd = [dataInd; jd*ones(n(jd),1)];
end


%Check dimensions etc.
n_in = size(u,1);
M = size(psiold,2);
rany = [rany; [min(u')',max(u')']];

%Check which variables are accounted for in each time series (excluding
%knocked-out genes in respective experiments)
Incl = ones(sum(n),size(ts,2));
if isfield(data,'ko')
    for jd = 1:nrData
        for jser = 1:size(ts,2)
            Incl(sum(n(1:jd-1))+data.ko{jd,jser},jser) = zeros(length(data.ko{jd,jser}),1);
        end
    end
end 

%Check genes that are not knocked out in at least one experiment (normally
%geneList = 1:n)
geneList = find(sum(Incl,2)>0)';
excludedGenes = setdiff(1:sum(n),geneList);

%Prior probability for the existence of a link (default = 1/n)
if isfield(parameters,'link_pr')
    log_link_pr = log(parameters.link_pr);
else
    log_link_pr = -log(sum(n))*ones(sum(n),sum(n)+n_in);
end
if numel(log_link_pr) == 1
    log_link_pr = log_link_pr*ones(sum(n),sum(n)+n_in);
end

%Process the prior network information
S_aux = zeros(sum(n),sum(n)+n_in);
if isfield(data,'sure')
    if size(data.sure,2) > sum(n) + .5
        S_aux = data.sure;
    else
        S_aux(1:sum(n),1:sum(n)) = data.sure;
    end
end
Sold = max(Sold,S_aux);
Sold = min(Sold,1+S_aux);
S_aux = ones(size(S_aux)) - abs(S_aux);
S_aux(:,excludedGenes) = 0;
S_aux(excludedGenes,:) = 0;
Sold(:,excludedGenes) = 0;
Sold(excludedGenes,:) = 0;


%Concatenate inputs as well
for j = 1:size(u,2)
    u = [u(:,1:(j-1)*nstep),u(:,(j-1)*nstep+1)*ones(1,nstep),u(:,(j-1)*nstep+2:end)];  
end

%Variables miss/nomiss shows which measurements are/aren't missing
nry = zeros(sum(n),1);
totvar = zeros(sum(n),1);
for jd = 1:nrData
    nomiss{jd} = ones(size(y{jd}));
    if isfield(data,'missing')
        if norm(size(data.missing)-size(data.ts))>.1
            error('Missing measurements indices are not consistent with the time series data!')
        end   
        for i = 1:n(jd)
            for jser = 1:size(Ser{jd},2)
                if size(data.missing{jd,jser},1) > .5
                    miss = find(abs(data.missing{jd,jser}(:,1)-i)<.5);
                    nomiss{jd}(i,data.missing{jd,jser}(miss,2)+Ser{jd}(1,jser)-1) = zeros(1,length(miss));
                end                
            end            
        end
    end
    
    %Auxiliary indices showing which points in concatenated trajectory are used
    %for probability distribution calculation and which points are compared
    %with data
    derind_full{jd} = (1:nstep*(Ser{jd}(2,1)-Ser{jd}(1,1)));
    yind{jd} = 1 + (0:Ser{jd}(2,1)-1)*nstep;
    d = reshape(ones(nstep,1)*(Tsam{jd,1}(2:end)-Tsam{jd,1}(1:end-1))/nstep,1,[]);
    d_full{jd} = d;    
    tx{jd} = [Tsam{jd,1}(1), Tsam{jd,1}(1) + cumsum(d)];
    for jser = 2:size(Ser{jd},2)
        derind_full{jd} = [derind_full{jd}, (1:nstep*(Ser{jd}(2,jser)-Ser{jd}(1,jser)))+1+derind_full{jd}(end)];
        yind{jd} = [yind{jd}, yind{jd}(end)+1+(0:Ser{jd}(2,jser)-Ser{jd}(1,jser))*nstep];
        d = reshape(ones(nstep,1)*(Tsam{jd,jser}(2:end)-Tsam{jd,jser}(1:end-1))/nstep,1,[]);
        d_full{jd} = [d_full{jd}, d];
        tx{jd} = [tx{jd}, [Tsam{jd,jser}(1), Tsam{jd,jser}(1) + cumsum(d)]];
    end
    
    
    %Calculate signal (quadratic) variation
    Total_time = 0;
    for l = 1:size(Ser{jd},2)
        nry(sum(n(1:jd-1))+(1:size(y{jd},1))) = nry(sum(n(1:jd-1))+(1:size(y{jd},1))) + sum(((y{jd}(:,Ser{jd}(1,l)+1:Ser{jd}(2,l))-y{jd}(:,Ser{jd}(1,l):Ser{jd}(2,l)-1)).^2./(Tsam{jd,l}(2:end)-Tsam{jd,l}(1:end-1))),2); 
        totvar(sum(n(1:jd-1))+(1:size(y{jd},1))) = totvar(sum(n(1:jd-1))+(1:size(y{jd},1))) + sum(abs(y{jd}(:,Ser{jd}(1,l)+1:Ser{jd}(2,l))-y{jd}(:,Ser{jd}(1,l):Ser{jd}(2,l)-1)),2);
        Total_time = Total_time + Incl((sum(n(1:jd-1))+1):sum(n(1:jd)),l).*(Tsam{jd,l}(end)-Tsam{jd,l}(1));
    end 
    nry(sum(n(1:jd-1))+(1:size(y{jd},1))) = nry(sum(n(1:jd-1))+(1:size(y{jd},1)))./Total_time;
    totvar(sum(n(1:jd-1))+(1:size(y{jd},1))) = totvar(sum(n(1:jd-1))+(1:size(y{jd},1)))./Total_time;
   
    
    %Embedding of the measurements to a piecewise linear function
    mm = max(Ser{jd}(4,:)-Ser{jd}(3,:))+1;
    Pr{jd} = zeros(mm,max(Ser{jd}(2,:)-Ser{jd}(1,:))+1);
    Pr{jd}(1:nstep,1) = flipud((1:nstep)'/nstep);
    Pr{jd}(mm-nstep+1:mm,end) = (1:nstep)'/nstep;
    for j = 2:max(Ser{jd}(2,:)-Ser{jd}(1,:))
        Pr{jd}((j-2)*nstep+2:(j-1)*nstep+1,j) = (1:nstep)'/nstep;
        Pr{jd}((j-1)*nstep+1:j*nstep,j) = flipud((1:nstep)'/nstep);
    end
    Pintc = sin((1:nstep-1)'*(1:nstep-1)/nstep*pi)./(pi*ones(nstep-1,1)*(1:nstep-1))*2^.5;
     
    yold{jd} = xs_old{jd}(:,yind{jd});
    xstore{jd} = zeros(size(xs_old{jd}));
     
end

%Interpolation matrices between data types (if nrData > 1)
for jd = 1:nrData
    for jin = [1:jd-1 jd+1:nrData]
        Pip{jd,jin} = interpol(tx{jin},tx{jd},Ser{jin},Ser{jd});
    end
end

%Store results and statistics
Plink = zeros(size(Sold));
chain = 0;
acctraj = 0;
acctop = zeros(sum(n),1);
acchyp = zeros(sum(n),1);
accr = zeros(sum(n),1);


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
        bets = (1-(parameters.ebeta)^2)^.5*betsold(i,:)+(parameters.ebeta)*randn(1,sum(n)+n_in);
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
        derind = derind_full{dataInd(i)};
        d = d_full{dataInd(i)};
        sers = 1:size(Ser{dataInd(i)},2);
        if sum(Incl(i,:)) < size(Ser{dataInd(i)},2) - .5
            derind = [];
            d = [];
            sers = find(Incl(i,:)>.5);
            for jser = sers
                derind = [derind,derind_full{dataInd(i)}((Ser{dataInd(i)}(3,jser):Ser{dataInd(i)}(4,jser)-1)-jser+1)];
                d = [d, d_full{dataInd(i)}((Ser{dataInd(i)}(3,jser):Ser{dataInd(i)}(4,jser)-1)-jser+1)];
            end
        end   
        N = length(derind);
        
        %Form covariance matrices
        KM = zeros(M,M);
        KNM = zeros(N,M);
        for j = find(S(1:sum(n))>.5)
            KM = KM + beta(j)*bsxfun(@minus,psiold(j,:),psiold(j,:)').^2;
            if dataInd(i) == dataInd(j)
                KNM = KNM + beta(j)*bsxfun(@minus,psiold(j,:),xs_old{dataInd(j)}(j-sum(n(1:dataInd(i)-1)),derind)').^2;
            else
                x1 = xs_old{dataInd(j)}(j-sum(n(1:dataInd(j)-1)),:)*Pip{dataInd(i),dataInd(j)};
                KNM = KNM + beta(j)*bsxfun(@minus,psiold(j,:),x1').^2; 
            end        
        end
        for j = sum(n) + find(S(sum(n)+1:sum(n)+n_in)>.5)
            KM = KM + beta(j)*bsxfun(@minus,psiold(j,:),psiold(j,:)').^2;
            KNM = KNM + beta(j)*bsxfun(@minus,psiold(j,:),u(j-sum(n),:)').^2;
        end
        KM = gamma_tr*exp(-KM);
        KNM = gamma_tr*exp(-KNM);
        
        %Compute the load
        KC = chol(KM+1/q(i)*((KNM'*(d'.*KNM))) + 1e-5*eye(M));
        der = ((xs_old{dataInd(i)}(i-sum(n(1:dataInd(i)-1)),derind+1)'-xs_old{dataInd(i)}(i-sum(n(1:dataInd(i)-1)),derind)') - d'.*(mbtr-matr*xs_old{dataInd(i)}(i-sum(n(1:dataInd(i)-1)),derind)'))/q(i);
        ld  =(KC'\(KNM'*der));
        
        
        %Part of the Wiener measure that is not in the CN sampler
        nrY = 0;
        for l = 1:size(Ser{dataInd(i)},2)   
            nrY = nrY + sum(((yold{dataInd(i)}(i-sum(n(1:dataInd(i)-1)),Ser{dataInd(i)}(1,l)+1:Ser{dataInd(i)}(2,l))-yold{dataInd(i)}(i-sum(n(1:dataInd(i)-1)),Ser{dataInd(i)}(1,l):Ser{dataInd(i)}(2,l)-1)).^2./(Tsam{dataInd(i),l}(2:end)-Tsam{dataInd(i),l}(1:end-1))),2);     
        end 
        
        %Cost function value
        J1 = .5*nrY/q(i) - .5*ld'*ld+log(prod(diag(KC))) - .5*log(det(KM+1e-5*eye(M))) - (mbtr-matr*xs_old{dataInd(i)}(i-sum(n(1:dataInd(i)-1)),derind)')'*(xs_old{dataInd(i)}(i-sum(n(1:dataInd(i)-1)),derind+1)'-xs_old{dataInd(i)}(i-sum(n(1:dataInd(i)-1)),derind)')/q(i) + .5/q(i)*norm(d'.^.5.*(mbtr-matr*xs_old{dataInd(i)}(i-sum(n(1:dataInd(i)-1)),derind)'))^2;
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
        if (r(i)/rtr)^(1+sum(nomiss{dataInd(i)}(i-sum(n(1:dataInd(i)-1)),:))/2)*exp(.00001./r(i)-.00001./rtr+sum((y{dataInd(i)}(i-sum(n(1:dataInd(i)-1)),find(nomiss{dataInd(i)}(i-sum(n(1:dataInd(i)-1)),:)>.5))-yold{dataInd(i)}(i-sum(n(1:dataInd(i)-1)),find(nomiss{dataInd(i)}(i-sum(n(1:dataInd(i)-1)),:)>.5))).^2)/2*(1/r(i)-1/rtr)) > rand
            r(i) = rtr;
            accr(i) = accr(i) + 1;
        end
    end
        
    
    %% Trajectory sampling
    
    %Sample the process noise variance
    qtr = q + parameters.eq*randn(size(q));
    qtr = .5e-5 + abs(qtr-.5e-5);
    
    %Sample the trajectory
    for jd = 1:nrData
        for l = 1:size(Ser{jd},2)
            if isfield(data,'missing')  %Missing measurements
                if size(data.missing{jd,l},1) > .5
                    for i = sum(n(1:jd-1))+1:sum(n(1:jd))
                        Csam = missing_data_sampler(data.missing{jd,l},Tsam{jd,l},r(i),qtr(i),i);
                        coef = (1-parameters.etraj^2)^.5*nomiss{dataInd(i)}(i-sum(n(1:dataInd(i)-1)),Ser{jd}(1,l):Ser{jd}(2,l)) + (qtr(i)/q(i))^.5*(1-nomiss{dataInd(i)}(i-sum(n(1:dataInd(i)-1)),Ser{jd}(1,l):Ser{jd}(2,l)));
                        yhat{jd}(i-sum(n(1:dataInd(i)-1)),Ser{jd}(1,l):Ser{jd}(2,l)) = y{dataInd(i)}(i-sum(n(1:dataInd(i)-1)),Ser{jd}(1,l):Ser{jd}(2,l)) + coef.*(yold{dataInd(i)}(i-sum(n(1:dataInd(i)-1)),Ser{jd}(1,l):Ser{jd}(2,l))-y{dataInd(i)}(i-sum(n(1:dataInd(i)-1)),Ser{jd}(1,l):Ser{jd}(2,l))) + parameters.etraj*randn(1,Ser{jd}(2,l)-Ser{jd}(1,l)+1)*Csam';
                    end
                else
                    yhat{jd}(:,Ser{jd}(1,l):Ser{jd}(2,l)) = y{jd}(:,Ser{jd}(1,l):Ser{jd}(2,l)) + (1-parameters.etraj^2)^.5*(yold{jd}(:,Ser{jd}(1,l):Ser{jd}(2,l))-y{jd}(:,Ser{jd}(1,l):Ser{jd}(2,l))) + parameters.etraj*diag(r(sum(n(1:jd-1))+1:sum(n(1:jd))).^.5)*randn(n(jd),Ser{jd}(2,l)-Ser{jd}(1,l)+1);   
                end
            else
                yhat{jd}(:,Ser{jd}(1,l):Ser{jd}(2,l)) = y{jd}(:,Ser{jd}(1,l):Ser{jd}(2,l)) + (1-parameters.etraj^2)^.5*(yold{jd}(:,Ser{jd}(1,l):Ser{jd}(2,l))-y{jd}(:,Ser{jd}(1,l):Ser{jd}(2,l))) + parameters.etraj*diag(r(sum(n(1:jd-1))+1:sum(n(1:jd))).^.5)*randn(n(jd),Ser{jd}(2,l)-Ser{jd}(1,l)+1);   
            end        
            xs{jd}(:,Ser{jd}(3,l):Ser{jd}(4,l)) = sparse(diag((qtr(sum(n(1:jd-1))+1:sum(n(1:jd)))./q(sum(n(1:jd-1))+1:sum(n(1:jd)))).^.5))*(1-parameters.etraj^2)^.5*xs_old{jd}(:,Ser{jd}(3,l):Ser{jd}(4,l)) + (yhat{jd}(:,Ser{jd}(1,l):Ser{jd}(2,l))-sparse(diag((qtr(sum(n(1:jd-1))+1:sum(n(1:jd)))./q(sum(n(1:jd-1))+1:sum(n(1:jd)))).^.5))*(1-parameters.etraj^2)^.5*yold{jd}(:,Ser{jd}(1,l):Ser{jd}(2,l)))*Pr{jd}(1:(Ser{jd}(4,l)-Ser{jd}(3,l)+1),1:(Ser{jd}(2,l)-Ser{jd}(1,l)+1))'; 
            xs{jd}(:,Ser{jd}(3,l)+1:Ser{jd}(4,l)) = xs{jd}(:,Ser{jd}(3,l)+1:Ser{jd}(4,l)) + parameters.etraj*sparse(diag(qtr(sum(n(1:jd-1))+1:sum(n(1:jd))).^.5))*((nstep*d_full{jd}((Ser{jd}(3,l):Ser{jd}(4,l)-1)-l+1)).^.5.*reshape([Pintc*randn(nstep-1,n(jd)*(Ser{jd}(2,l)-Ser{jd}(1,l)));zeros(1,n(jd)*(Ser{jd}(2,l)-Ser{jd}(1,l)))],[],n(jd))'); 
        end
    end
    
    % Sample the pseudoinputs, and "mirror" them to the box containing the data. 
    % This corresponds to uniform prior in the box.
    psin = psiold + .025*randn(size(psiold));
    psin = min(psin,2*rany(:,2)-psin);
    psin = max(psin,2*rany(:,1)-psin);
    
    %Initialise the cost function value
    J1 = zeros(sum(n),1);
    
    for i = geneList
        jd = dataInd(i);
        
        %Exclude the knockout data of gene i (if any)
        derind = derind_full{jd};
        d = d_full{jd};
        if sum(Incl(i,:)) < size(Ser{jd},2) - .5
            derind = [];
            d = [];
            sers = find(Incl(i,:)>.5);
            for jser = sers
                derind = [derind,derind_full{jd}((Ser{jd}(3,jser):Ser{jd}(4,jser)-1)-jser+1)];
                d = [d,d_full{jd}((Ser{jd}(3,jser):Ser{jd}(4,jser)-1)-jser+1)];
            end
        end   
        N = length(derind);
        
        %Compute the covariances
        KM = zeros(M,M);  %pseudoinputs
        KNM = zeros(N,M); %pseudoinput vs. trajectory
        for j = find(Sold(i,1:sum(n))>.5)
            KM = KM + betsold(i,j)*bsxfun(@minus,psin(j,:),psin(j,:)').^2;
            if dataInd(j) == jd
                KNM = KNM + betsold(i,j)*bsxfun(@minus,psin(j,:),xs{jd}(j-sum(n(1:dataInd(j)-1)),derind)').^2;
            else
                x1 = xs{dataInd(j)}(j-sum(n(1:dataInd(j)-1)),:)*Pip{dataInd(i),dataInd(j)};
                KNM = KNM + betsold(i,j)*bsxfun(@minus,psin(j,:),x1').^2; 
            end 
        end
        for j = sum(n)+find(Sold(i,sum(n)+1:sum(n)+n_in)>.5)
            KM = KM + betsold(i,j)*bsxfun(@minus,psin(j,:),psin(j,:)').^2;
            KNM = KNM + betsold(i,j)*bsxfun(@minus,psin(j,:),u(j-sum(n),:)').^2;
        end
        KM = gamma(i)*exp(-KM);
        KNM = gamma(i)*exp(-KNM);
        
        %Compute the load 
        KC = chol(KM+1/qtr(i)*((KNM'*(d'.*KNM)))+1e-5*eye(M));
        der = ((xs{jd}(i-sum(n(1:jd-1)),derind+1)'-xs{jd}(i-sum(n(1:jd-1)),derind)') - (d'.*(mb(i)-ma(i)*xs{jd}(i-sum(n(1:jd-1)),derind)')))/qtr(i);
        ld = (KC'\(KNM'*der));          

        %Part of the Wiener measure that is not in the CN sampler
        nrY = 0;
        for l = 1:size(Ser{jd},2)
            nrY = nrY + sum(((yhat{jd}(i-sum(n(1:jd-1)),Ser{jd}(1,l)+1:Ser{jd}(2,l))-yhat{jd}(i-sum(n(1:jd-1)),Ser{jd}(1,l):Ser{jd}(2,l)-1)).^2./(Tsam{jd,l}(2:end)-Tsam{jd,l}(1:end-1))),2);     
        end 
        
        %Cost function value
        J1(i) = .5*nrY/qtr(i) - .5*ld'*ld+log(prod(diag(KC))) - .5*log(det(KM+1e-5*eye(M))) - (mb(i)-ma(i)*xs{jd}(i-sum(n(1:jd-1)),derind)')'*(xs{jd}(i-sum(n(1:jd-1)),derind+1)'-xs{jd}(i-sum(n(1:jd-1)),derind)')/qtr(i) + .5/qtr(i)*norm(d'.^.5.*(mb(i)-ma(i)*xs{jd}(i-sum(n(1:jd-1)),derind)'))^2;    
    end
       
    %Accept or reject the new sample
    P_aux_q = 1;
    for jd = 1:nrData
        P_aux_q = P_aux_q*exp(sum(.00001./q(sum(n(1:jd-1))+1:sum(n(1:jd)))-.00001./qtr(sum(n(1:jd-1))+1:sum(n(1:jd)))+(log(q(sum(n(1:jd-1))+1:sum(n(1:jd))))-log(qtr(sum(n(1:jd-1))+1:sum(n(1:jd))))).*(1.001+.5*(size(y{jd},2)-size(Ser{jd},2)))));
    end
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
        for jd = 1:nrData
            xstore{jd} = xstore{jd} + xs_old{jd};
        end
        
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
for jd = 1:nrData
    xstore{jd} = xstore{jd}/(parameters.its/10);
end
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

