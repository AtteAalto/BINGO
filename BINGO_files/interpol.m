function P =  interpol(tx,tpol,Serx,Serpol)

%xj is a trajectory given at time points tx, and it is interpolated to
%times tpol.

%xpol = zeros(1,length(tpol)-size(Serx,2));

P = sparse(zeros(length(tx),length(tpol)-size(Serpol,2)));

jaux = 0;
for jser = 1:size(Serx,2)

    txser = tx(Serx(3,jser):Serx(4,jser));
    tpolser = tpol(Serpol(3,jser):Serpol(4,jser)-1);

    if jser == 1
        add = 0;
    else
        add = Serx(4,jser-1);
    end
    
    
    for jx = 1:length(tpolser)
        jaux = jaux+1;
        
        tind = sum(tpolser(jx) > txser);
        
        if tind == 0
            P(add + 1,jaux) = 1;
        elseif tind == length(txser)
            P(Serx(4,jser),jaux) = 1;
        else
            cc = (tpolser(jx)-txser(tind))/(txser(tind+1)-txser(tind));
            
            P(add + tind + 1, jaux) = cc;
            P(add + tind, jaux) = 1 - cc;
        end

    end
       
end
        
        
    










