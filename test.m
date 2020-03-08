x = randn(1,1000, 'int8');
y = zeros(1,400);
for i_x = 1 : 1000
  if abs(x(1,i_x) <= 200
    y(1, x(1,i_x)+200) = y(1, x(1,i_x)+200) .+ 1;
  end;
end;
figure(1);
plot(1:400, y);