function no_output=disp_stats(message,stats,chain,its)

n=size(stats.acctop,1);
disp(' ')
disp('*******************')
disp(message)    
disp(['***********************************************************************************************'])
a=stats.start_time;
m1=':'; m2=':'; if floor(a(5))<9.5; m1=':0'; end;  if floor(a(6))<9.5; m2=':0'; end 
disp(['* Start time: ' num2str(floor(a(4))) m1 num2str(floor(a(5))) m2 num2str(floor(a(6)))])
a=stats.fin_time;
m1=':'; m2=':'; if floor(a(5))<9.5; m1=':0'; end;  if floor(a(6))<9.5; m2=':0'; end 
disp(['* End time: ' num2str(floor(a(4))) m1 num2str(floor(a(5))) m2 num2str(floor(a(6)))])
clear('m1','m2','a')
disp(['* Number of samples: ' num2str(chain)])
disp(['* Acceptance probability for trajectory and q: ' num2str(round(stats.acctraj/its*100,1))])
disp(['* Average/minimum acceptance probability for topology: ' num2str(round(sum(stats.acctop)/its/n*100,1)) ' / ' num2str(round(min(stats.acctop)/its*100,1))])
disp(['* Average/minimum acceptance probability for (gamma, beta, a, b): ' num2str(round(sum(stats.acchyp)/its/n*100,1)) ' / ' num2str(round(min(stats.acchyp)/its*100,1)) '  r: ' num2str(round(sum(stats.accr)/its/n*100,1)) ' / ' num2str(round(min(stats.accr)/its*100,1))])
disp(['***********************************************************************************************'])
disp([' '])


