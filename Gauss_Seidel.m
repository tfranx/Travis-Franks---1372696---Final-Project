%http://searchdl.org/public/book_series/elsevierst/7/7.pdf
%Travis Franks 1372696 Helmholtz Equation Final Project Gauss Seidel Method
clear all
clc
%Setting the value of N, the number of elements per row (along x domain):
N = input('Enter value of N, the number of elements per row: ');
%Setting the value of M, the number of elements per column (along y
%domain):
M = input('Enter value of M, the number of elements per column: ');
%Setting the value of C (capital lambda in problem statement):
C = input('Enter value of C, the given constant for capital lambda: ');
%Setting the value of Es, the acceptable limit of error for system convergence:
Es = 10^-10;
%Defining L, the length of the X and Y domains:
L = 2 * pi();
%Determining DX and DY, the change of X and Y respectively from element to
%element:
DX = L / (N + 1);
DY = L / (M + 1);
%Defining coefficients, A and B, that remain constant throughout matrix operations:
A = C - (2 * (DX^2)) - (2 * (DY^2));
B = (DX^2) * (DY^2);
%Determining the sizes, Element_Count, of the matrices to be evaluated:
Element_Count = (N+2) * M; %+2 comes from the two Neumann boundary conditions imposing ghost nodes
%Solving for F(X,Y) values, "F", arrayed as a matrix:
F = zeros(Element_Count, 1);
for j = 1:M
    for i = 1:(N+2)
        F(i + (N+2) * (j-1), 1) = cos((pi() / 2) * ((((i - 1) * DX) / pi()) + 1)) * sin((j * DY) / 2);
    end
end
% Determining values of U determined from the Dirichlet Boundary Conditions
% imposed onto the Y-domain:
Dirichlet_Val_Lower = zeros(N+2, 1);
Dirichlet_Val_Upper = zeros(N+2, 1);
for i = 1:(N+2)
    Dirichlet_Val_Lower(i, 1) = cos(pi() * DX * (i-1)) * cosh(2 * pi() - (DX * (i-1))); %Applies to known U values determined from the Dirichlet Boundary Condition imposed at Y = -pi
    Dirichlet_Val_Upper(i, 1) = (((i-1) * DX)^2) * sin(((i-1) * DX) / 4); %Applies to known U values determined from the Dirichlet Boundary Condition imposed at Y = pi
end
%Generating original K (left hand side coefficients) matrix:
K = zeros(Element_Count, Element_Count);
for j = 1:M %i + (N+2) * (j-1) for indexing the column values
    for i = 1:(N+2)
        K(i + (N+2) * (j-1), i + (N+2) * (j-1)) = A; %Fills the diagonal with the A values
        if (i == 1)
            K(i + (N+2) * (j-1), i + (N+2) * (j-1) + 1) = 2 * (DY^2); %Fills all points in matrix where the Neumann Boundary Condition for the left end of the x domain produces 2*DY^2
        elseif (i == N+2)
            K(i + (N+2) * (j-1), i + (N+2) * (j-1) - 1) = 2 * (DY^2); %Fills all points in matrix where the Neumann Boundary Condition for the right end of the x domain produces 2*DY^2
        else
            K(i + (N+2) * (j-1), i + (N+2) * (j-1) + 1) = DY^2; %Fills the upper portion of the center tridiagonal with DY^2, where applicable
            K(i + (N+2) * (j-1), i + (N+2) * (j-1) - 1) = DY^2; %Fills the lower portion of the center tridiagonal with DY^2, where applicable
        end
        if (j > 1)
            K(i + (N+2) * (j-2), i + (N+2) * (j-1)) = DX^2; %Fills the uppermost diagonal with values of DX^2
        end
        if (j < M)
            K(i + (N+2) * (j), i + (N+2) * (j-1)) = DX^2; %Fills the lowest diagonal with values of DX^2
        end
    end
end
%Generating original P (right hand side known values) matrix:
P = zeros(Element_Count, 1);
for i = 1:Element_Count
    if (i <= (N+2))
        P(i,1) = B * F(i, 1) - ((DX^2) * Dirichlet_Val_Lower(i, 1)); %All values impacted by lower Y domain (at Y = -pi) Dirichlet Boundary Condition
    elseif (i > (Element_Count - (N+2)))
        P(i,1) = B * F(i, 1) - ((DX^2) * Dirichlet_Val_Upper(i - ((N+2) * (M-1)), 1));  %All values impacted by upper Y domain (at Y = pi) Dirichlet Boundary Condition
    else
        P(i,1) = B * F(i, 1); %All values not impacted by Dirichlet Boundary Conditions
    end
end
save variables.mat
%%
load variables.mat
%Performing iterations of Gauss-Seidel Method to solve for U values
%(unknown solution values over X and Y domains):
U = zeros(Element_Count, 1);
W = zeros(Element_Count, 1);
Z = 0; %Z functions as a counter for number of iterations performed during Gauss-Seidel
Error = zeros(1, Element_Count);
Ea = 100; %Provides initial value of Ea, or the relative iterative error, for the Gauss_Seidel Approximation
if (A ~= 0)%Set condition for just in case the center diagonal is zero from a poor choice in nodes along X and Y domains
while (Ea > Es)
    for i = 1:Element_Count
        W(i,1) = U(i,1); %Provides value of U from previous iteration to be used in evaluating the relative error
        U(i,1) = P(i,1) / K(i,i); %Known value divided by the desired U-value's associated coefficient, from which iterations will subtract from for Gauss-Seidel Approximation
        for j = 1:Element_Count
            if (i==j)
                U(i,1) = U(i,1); %Keeps the loop from subtracting the necessary P-value from itself when performing iterations
            else
                U(i,1) = U(i,1) - (K(i,j) * U(j,1) / K(i,i)); %Subtracts each successive U value (divided by the coefficient of the desired U value) from the equation (in essence the heart of the Gauss-Seidel Approximation)
            end
        end
        Error(1,i) = abs((U(i,1) - W(i,1)) / U(i,1)); %Computes the relative error of each U value individually
    end
    Ea = max(Error); %Uses maximum relative error of iteration as operational relative error value to be used to satisfy the loop condition
    Z = Z + 1; %Counts number of iterations
end
else
    disp('Cannot compute as coefficient matrix center diagonal = 0. Consider entering different number of nodes for X and/or Y domain.')
end
Answer = P\K; %Provides exact solution for comparison (need to verify that it is correct as well)
%%
%Gauss Seidel Approximation is outputting divergent (approches infinity)
%values for the U matrix. Consider indexing loike the matrix setup and
%separating into independent blocks to aid in grid convergence.
