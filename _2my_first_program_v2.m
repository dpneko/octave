T = 4
x = zeros(T,1);
x(1) = 10;
for i = 2:T
  x(i) = x(i-1) + 1;
end
% x
sum(x)



