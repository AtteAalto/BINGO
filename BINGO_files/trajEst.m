function data = trajEst(data)

%Set value for outlier tolerance (default = 5)
outlierTolerance = 5;
if isfield(data,'outlierTolerance')
    outlierTolerance = data.outlierTolerance;
end

dispAux = true;

ndim = size(data.ts{1},1);
for jser = 1:length(data.ts)
    OL = zeros(size(data.ts{jser}));
    
    for iEst = data.notMeasured{jser}
        yest = zeros(1,size(data.ts{jser},2));
        naux = 0;
        
        %Estimate the missing gene based on each time series separetely
        for jl = [1:jser-1 jser+1:length(data.ts)]
            
            %Check if gene iEst is measured in time series jl
            if isempty(data.notMeasured{jl}) || min(abs(iEst-data.notMeasured{jl})) > .5
            
                %Check common genes in the two time series
                ii1 = setdiff([1:iEst-1 iEst+1:ndim],data.notMeasured{jser});
                ii2 = setdiff([1:iEst-1 iEst+1:ndim],data.notMeasured{jl});
                iCommon = intersect(ii1,ii2);

                %If iEst is the only common gene, it cannot be estimated
                %using time seriel jl. In that case, skip that time series.
                if isempty(iCommon)
                    continue
                end
                    
                %Form regression matrices
                Scale = mean(mean(data.ts{jser}(iCommon,:)))/mean(mean(data.ts{jl}(iCommon,:)));
                XX = [Scale*data.ts{jl}(iCommon,:); ones(1,size(data.ts{jl},2))];
                XXc = XX;
                XXc(1:end-1,:) = XXc(1:end-1,:) - sum(XXc(1:end-1,:),2)/size(XXc,2);
                YY = Scale*data.ts{jl}(iEst,:);
                
                %Estimate the trajectory by linear regression. The Tikhonov
                %regularisation parameter is selected using the method
                %presented in G. Golub, M. Heath, and G. Wahba (1979). 
                %"Generalized cross-validation as a method for choosing a 
                %good ridge parameter", Technometrics 21(2): 215–223.
                
                %Tikhonov parameter determination
                [U,S,~] = svd(XXc');
                q = size(XXc,1);
                m = size(XXc,2);
                U = U(:,1:q);
                S = diag(S);
                RSS0 = sum((YY' - U*U'*YY').^2);
                UTb = U'*YY';
                Jcost = @(x)((RSS0 + x^4*sum(UTb.^2./(S.^2+x^2).^2))/(m - q + x^2*sum(1./(S.^2+x^2)))^2);
                alpha = fminsearch(Jcost,.1);
                
                %Outlier handling
                Xin = [data.ts{jser}(iCommon,:); ones(1,size(data.ts{jser},2))];
                Xinc = Xin;
                Xinc(1:end-1,:) = Xinc(1:end-1,:) - sum(Xinc(1:end-1,:),2)/size(Xinc,2);
                B = zeros(length(iCommon),length(iCommon)+1);
                Cdiff = zeros(length(iCommon),size(data.ts{jser},2));
                for jin = 1:length(iCommon)
                    aaAux = (XXc([1:jin-1 jin+1:end],:)*XXc([1:jin-1 jin+1:end],:)' + alpha^2*eye(size(XXc,1)-1))^-1*XXc([1:jin-1 jin+1:end],:)*XX(jin,:)';
                    B(jin,[1:jin-1 jin+1:end]) = aaAux';
                end
                Xdiff = - Xin(1:end-1,:) + B*Xinc;
                for jt = 1:size(Xdiff,2)
                    Cdiff(:,jt) = (Xdiff(:,jt)'*(B(:,1:end-1)*diag(Xinc(1:end-1,jt))))';
                end
                Cdiff = Cdiff./sum(Xdiff.^2,1).^.5*size(Xdiff,1)^.5/mean(mean(Xin(1:end-1,:)));
                Xcoeff = Cdiff / outlierTolerance;
                Xcoeff = max(1,Xcoeff);
                Xinc(1:end-1,:) = Xinc(1:end-1,:)./Xcoeff;

                %Keep track of outliers
                OL(iCommon,:) = max(OL(iCommon,:),(Xcoeff - 1));
                
                %Trajectory estimation with the linear (static) model
                aa = (XXc*XXc' + alpha^2*eye(size(XXc,1)))^-1*XXc*YY';
                yestTemp = aa'*Xinc;
                WW = sum((YY - aa'*XX).^2);
                
                %Calculate the weighted average based on "training data" 
                %residuals of individual time series
                yest = yest + yestTemp/WW;
                naux = naux + 1/WW;
                
            end
            
        end
        
        if naux == 0
            error(['Data is not suffcicient for estimating gene ' num2str(iEst) '!'])
        end
        data.ts{jser}(iEst,:) = yest/naux;
        
        %Optional smoothing
        if isfield(data,'smoothing')
            ySmooth = zeros(1,size(data.ts{jser},2));
            ySmooth(1) = (1-data.smoothing)*data.ts{jser}(iEst,1) + data.smoothing*data.ts{jser}(iEst,2);
            ySmooth(end) = (1-data.smoothing)*data.ts{jser}(iEst,end) + data.smoothing*data.ts{jser}(iEst,end-1);
            for jt = 2:length(ySmooth)-1
                ySmooth(jt) = (1-2*data.smoothing)*data.ts{jser}(iEst,jt) + data.smoothing*data.ts{jser}(iEst,jt-1) + data.smoothing*data.ts{jser}(iEst,jt+1);
            end
            data.ts{jser}(iEst,:) = ySmooth;
        end
        
        %Enforcing positivity
        if isfield(data,'positive') && data.positive
            if min(data.ts{jser}(iEst,:)) < 0
                negMean = sum(data.ts{jser}(iEst,data.ts{jser}(iEst,:) < 0))/size(data.ts{jser},2);
                data.ts{jser}(iEst,:) = data.ts{jser}(iEst,:) - negMean;
                data.ts{jser}(iEst,:) = max(data.ts{jser}(iEst,:),0);
                
            end
        end
        
    end
    
    %Display outlier information
    if max(max(OL))>0
        if dispAux
            dispAux = false;
            disp('Following outliers handled in estimating missing genes:')
            disp('* Series * Gene * Time point * Score')
        end
        
        for jg = 1:size(OL,1)
            iAux = find(OL(jg,:)>0);
            for is = 1:length(iAux)
                jserStr = num2str(jser);
                fill1 = repmat(' ',1,6-length(jserStr));
                
                jgStr = num2str(jg);
                fill2 = repmat(' ',1,5-length(jgStr));
                
                iaStr = num2str(iAux(is));
                fill3 = repmat(' ',1,11-length(iaStr));
                
                disp(['* ' fill1 jserStr ' *' fill2 jgStr ' *' fill3 iaStr ' * ' num2str(1+OL(jg,iAux(is)))])
            end
        end
    end
    
    
end


