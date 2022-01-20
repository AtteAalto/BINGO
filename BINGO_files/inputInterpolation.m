function u = inputInterpolation(data, nrTS, nstep)

u = [];

%Zero-order-hold
if strcmp(data.inputInterpolation,'ZOH')
    for j = 1:nrTS
        u = [u,data.input{j}(:,1:end-1)];
    end
    for j = 1:size(u,2)
        u = [u(:,1:(j-1)*nstep),u(:,(j-1)*nstep+1)*ones(1,nstep),u(:,(j-1)*nstep+2:end)];  
    end
end

%Linear interpolation
if strcmp(data.inputInterpolation,'linear')
    for j = 1:nrTS
        for jt = 1:size(data.input{j},2)-1
            u = [u, data.input{j}(:,jt) + (data.input{j}(:,jt+1)-data.input{j}(:,jt))*(0:nstep-1)/nstep];
        end
    end
end
         
%Minimal Sobolev-2 norm    
if strcmp(data.inputInterpolation,'Sobolev')
    maxdeg = 50;
    for j = 1:nrTS
        tFine = data.fine_times{j};
        tCoarse = tFine(1:nstep:end);
        
        VDM = [ones(1,length(tCoarse)); tCoarse-mean(tCoarse)];
        VDMbig = [ones(1,length(tFine)); tFine-mean(tFine)];
        ranT = max(tCoarse) - min(tCoarse);
        minT = min(tCoarse);
        
        for jdeg = 1:maxdeg
            VDM = [VDM; sin(2*jdeg*pi*(tCoarse-minT)/ranT); cos(2*jdeg*pi*(tCoarse-minT)/ranT)];
            VDMbig = [VDMbig; sin(2*jdeg*pi*(tFine-minT)/ranT); cos(2*jdeg*pi*(tFine-minT)/ranT)];
        end
        
        u0 = zeros(size(data.input{j},1),length(tFine));
        for jn = 1:size(data.input{j},1)

            aa = (VDM*VDM' + 1e-10*(2*pi)^4*ranT/2*diag([0 0 kron(1:maxdeg,[1 1])].^4))^-1*VDM*data.input{j}(jn,:)'; 
            u0(jn,:) = aa'*VDMbig;
        end
    
        u = [u, u0(:,1:end-1)];
    end
    
end

