function[power2]= Copy_of_water_filling(Ptot,noise)
% output: the power levels of each parallel channel

Pused=0;
Pleft=Ptot;
%N=input('# of channels:');
noise=1./noise;
temp=noise; %store unsorted version
noise=sort(noise);
N= length(noise);
power=zeros(1,N);

for i=1:N
   if Pleft == 0
      i=N+1; 
   elseif(i==N)
     for k=1:i
     power(k)=power(k)+ Pleft/i;
     end 
   elseif(Pleft >= i*(noise(i+1)-noise(i)))
    for k=1:i
     power(k)=power(k)+ noise(i+1)-noise(i);
    end  
   else
    for k=1:i
     power(k)=power(k)+ Pleft/i;
    end
   end
 Pused=sum(power);
 Pleft=Ptot-Pused;
end 

% convert the power vector to initial version;
power2=zeros(1,N);
for t=1:N
   for h=1:N
      if temp(t)==noise(h)
         power2(t)=power(h);
     end   
   end    
end      
a(2,:)=power2;a(1,:)=temp;
b=a';
%bar(b,'stack')


