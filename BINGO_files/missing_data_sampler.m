function C=missing_data_sampler(missing,Tsam,r,q,jgene)

C=r^.5*eye(length(Tsam));

miss=find(abs(missing(:,1)-jgene)<.5);
ind_miss=sort(missing(miss,2));
if length(ind_miss)>.5

    miss_start=sum(ind_miss'-(1:length(ind_miss))<.5);
    miss_end=sum((length(Tsam)-length(ind_miss)+1:length(Tsam))-ind_miss'<.5);
    
    for jj=1:miss_start
        C(1:jj,jj)=q.^.5*(Tsam(jj+1)-Tsam(jj))^.5*ones(jj,1);
    end
    for jj=1:miss_end
        C(end-jj+1,end-jj+1:end)=q.^.5*(Tsam(end-jj+1)-Tsam(end-jj))^.5*ones(jj,1);
    end

    jtime=miss_start+1;
    
    while length(ind_miss(jtime:end-miss_end))>.5
        nr_miss=sum(abs(ind_miss(jtime:end)'-(ind_miss(jtime):(ind_miss(jtime)+length(ind_miss(jtime:end))-1)))<.5);    %    jtime:jtime+length(ind_miss)-miss_start-1))<.5);
        
        cov=q*(Tsam(ind_miss(jtime)+nr_miss)-max(Tsam(ind_miss(jtime):ind_miss(jtime)+nr_miss-1),Tsam(ind_miss(jtime):ind_miss(jtime)+nr_miss-1)')).*(min(Tsam(ind_miss(jtime):ind_miss(jtime)+nr_miss-1),Tsam(ind_miss(jtime):ind_miss(jtime)+nr_miss-1)')-Tsam(ind_miss(jtime)-1))/(Tsam(ind_miss(jtime)+nr_miss)-Tsam(ind_miss(jtime)-1));       
        C(ind_miss(jtime):ind_miss(jtime)+nr_miss-1,ind_miss(jtime):ind_miss(jtime)+nr_miss-1)=chol(cov)';
            
        jtime=jtime+nr_miss;
    end
end

