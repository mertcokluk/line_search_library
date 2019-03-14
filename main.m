% Illustration of several renowned line search methods 
% Straightforward objective function J=0.5 x'Ax - b'x and descent direction d=-dJ/dx=grad(J(x)) = Ax-b

A=[2,1;0,1];            %Define the objective function and initial guess vector
b=[0.5;0.3];
x=[1;1];
cost = 0.5*x'*A*x - b'*x;

for iter=1:20
	g = A*x-b;                %Calculate gradient
	d = -1*g;                 %Descent direction: -grad
	x = x+0.5*d;              %update the solution
	cost = 0.5*x'*A*x - b'*x
end
