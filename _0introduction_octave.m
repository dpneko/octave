%% Introduction to Octave
% -----------------------

% Note: You can type in any of those commands
% in the command window of Octave,
% this is good practice for learning the syntax.
% 
% Alternatively, you can highlight sections of commands in
% this document, and use "Run Selection" in the context menu,
% or press the function key F9 to do so.

% ASSIGNMENTS
x = 1;
y = 2;
x = x + 3;

% Column vector
%创建由0组成的列矩阵
x = zeros(10,1);

% Inspect elements
x(3)
% x(20)  % meant to generate an instructive error

% Row vector
x = zeros(1,5);

% Matrix
y = zeros(5,5);

% Inspect elements
y(4,4)

% Assign a new value to an element
x   % inspect
x(2) = 1
x   % inspect again
y(2,4) 
y(2,4) = 7;
y   % inspect again

i_row = 2;
i_col = 4; 
y(i_row,i_col)
y(i_row,i_col) = 8;

i_row = 3;
y(i_row,i_col) = 8;
y        % inspect again

i_col = 2;
y(i_row,i_col) = 17;
y        % inspect again

% Systematically combining assignments
x = zeros(5,1);

x(1) = 1;
x(2) = x(1) + 1;
x(3) = x(2) + 1;

i = 4;
x(i) = x(i-1) + 1;
i = 5;
x(i) = x(i-1) + 1;

% A simple for loop
for i = 4:5;
    x(i) = x(i-1) + 1;
end;

% A simple for loop with a different starting value
x = zeros(5,1);
x(1) = 21;
for i = 2:5;
    x(i) = x(i-1) + 1;
end;


% Dimensions of a matrix
length(x)
size(x,1)
size(x,2)
size(y)     % get size in all dimensions

% Note: Make sure to check, use, and experiment with,
% the separate small programs
% - my_first_program.m
% - my_first_program_v2.m
% which illustrate the concepts of
% for-loops, and if-conditions
% ----------------------------


% There are specialized commands for many things,
% as you will see in the programs we distribute for models in the course
% For instance, calculating an evenly spaced grid can conveniently be done
% with a specialized command
my_favorite_grid = linspace(1,5,10);

% MATRICES
% Generate a 5x5 matrix consisting of random numbers
X = rand(5,5);

% Matrix operations are particularly nice in Octave
% Matrix multiplication
X2 = X*X;

% Multiplication of a matrix with a vector
y = ones(1,5)
% y2 = X*y          % this generates an instructive error

% Use a transposed vector/matrix
y = y' 
y2 = X*y

% Inverse of a matrix
X_inv = inv(X);
check_identity = X*X_inv;

% format long
% if the previous line is not commented out,
% you can see errors some digits after the decimal
check_identity  % This should correspond to the identity matrix.

format short
% Element-by-element multiplication
W = X.*X;

