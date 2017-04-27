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
