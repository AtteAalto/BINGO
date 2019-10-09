function y_out=missing_data_interpolation(data)

for jser=1:size(data.ts,2)
    y=data.ts{jser};
    if size(data.missing{jser},1)>.5    
        for jgene=1:size(y,1)
            miss=find(abs(data.missing{jser}(:,1)-jgene)<.5);
            ind_miss=sort(data.missing{jser}(miss,2));
            if length(ind_miss)>.5

                miss_start=sum(ind_miss'-(1:length(ind_miss))<.5);
                miss_end=sum((size(y,2)-length(ind_miss)+1:size(y,2))-ind_miss'<.5);
                y(jgene,1:miss_start)=y(jgene,miss_start+1)*ones(1,miss_start);
                y(jgene,size(y,2)-miss_end+1:end)=y(jgene,end-miss_end)*ones(1,miss_end);

                jtime=miss_start+1;
                while jtime<size(y,2)-miss_end+.5
                    nr_miss=sum(abs(ind_miss(miss_start+1:end)'-(jtime:jtime+length(ind_miss)-miss_start-1))<.5);

                    if nr_miss>.5
                        y(jgene,jtime:jtime+nr_miss-1)=y(jgene,jtime-1)*(1-(1:nr_miss)/(nr_miss+1))+y(jgene,jtime+nr_miss)*(1:nr_miss)/(nr_miss+1);

                        miss_start=miss_start+nr_miss;
                        jtime=jtime+nr_miss;
                    end
                    jtime=jtime+1;
                end
            end
        end
    end   
    y_out{jser}=y;
end

    
   






