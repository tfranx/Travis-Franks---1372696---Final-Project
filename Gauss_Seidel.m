%Travis Franks 1372696 Helmholtz Equation Final Project Gauss Seidel Method
clear all
clc
%Setting the value of X_Internal_Nodes, the number of elements along X-domain:
X_Internal_Nodes = input('Enter value of X_Internal_Nodes, the number of internal nodes for the X-domain: ');
%Setting the value of Y_Internal_Nodes, the number of elements along
%Y-domain:
Y_Internal_Nodes = input('Enter value of Y_Internal_Nodes, the number of internal nodes for the Y-domain: ');
%Setting the value of C (capital lambda in problem statement):
C = input('Enter value of C, the given constant for capital lambda: ');
%Setting the value of Es, the acceptable limit of error for system convergence:
Es = 10^-10;
%Defining L, the length of the X and Y domains:
L = 2 * pi();
%Determining DX and DY, the change of X and Y respectively from element to
%element:
DX = L / (X_Internal_Nodes + 1);
DY = L / (Y_Internal_Nodes + 1);
%Defining coefficients, A and B, that remain constant throughout matrix operations:
A = C - (2 * (DX^2)) - (2 * (DY^2));
B = (DX^2) * (DY^2);
%Determining number of elements  on both X and Y domains, so as to only
%perform the calculation one time for optimization:
N = X_Internal_Nodes + 2; %N = number of X domain values
M = Y_Internal_Nodes + 2; %M = number of Y domain values
%Performing setup for Gauss-Seidel Approximation that will solve for U values
%(unknown solution values over X and Y domains):
U = zeros(N,M);
W = zeros(N,M);
Z = 0; %Z functions as a counter for number of iterations performed during Gauss-Seidel
Error = zeros(1, N * (M-2)); %(M-2) rather than just M because the Dirichlet Boundary Conditions cause two rows to have constant values and therefore will have an error of 0 per iteration (optimal to exclude unnecessary repeated calculations) 
Ea = 100; %Provides initial value of Ea, or the relative iterative error, for the Gauss_Seidel Approximation
%Solving for U values defined by Dirichlet Boundary Conditions that will
%remain constant as Gauss Seidel iterations are performed:
for i = 1:N
    U(i,1) = cos(pi() * DX * (i-1)) * cosh((2 * pi()) - (DX * (i-1)));
    U(i,N) = ((i-1) * DX)^2 * sin(((i-1) * DX) / 4);
end
%Defining node indexing points to be used for the general form of the
%discretization:
NN = N-1;
MM = M-1;
%Performing Gauss Seidel Approximation to solve for U values:
if (A ~= 0)%Set condition for just in case the variable coefficient is equal to zero from a poor choice in nodes along X and Y domains, as it will function as a denominator
    while (Ea > Es)
        %Evaluating for general expression untouched by boundary
        %conditions:
        for j = 2:MM
            for i = 2:NN
                U(i,j) = (B
