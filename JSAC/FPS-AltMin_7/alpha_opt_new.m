function [alpha,objt,S] = alpha_opt_new( x )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% y=x(1:3,1:3);
% [z,ord1] = sort([y(:);0],'ascend');
[x,ord] = sort([x(:);0],'ascend'); %ord---元素在x中的位置，列
n = length(x);
ave = zeros(n-1,1);

ave(1) = x(1);
ave(end) = x(end);
k = find(x==0,1);
for i = 2:k-1
    ave(i) = (ave(i-1)*(i-1)+x(i))/i;
end
for i = n-2:-1:k
    ave(i) = (ave(i+1)*(n-i-1)+x(i+1))/(n-i);
end

idx = find(ave>=2*x(1:n-1) & ave<=2*x(2:n));  %(16)中找i

obj = zeros(2,length(idx));
for i = 1:length(idx)
    alpha = ave(idx(i));
    if(alpha<=0)
        s(:,i) = double(x-alpha/2<=0);
        obj(:,i) = [alpha;norm(x-alpha*s(:,i),2)^2];
    else
        s(:,i) = double(x-alpha/2>0);
        obj(:,i) = [alpha;norm(x-alpha*s(:,i),2)^2];
    end
end
[objt,loc] = min(obj(2,:));  %找到函数
% objt = objt-norm(x,2)^2;
alpha = obj(1,loc);
temp = s(:,loc);
S(ord) = temp;
S(end)=[];
end