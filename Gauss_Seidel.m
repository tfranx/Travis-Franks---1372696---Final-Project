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
%Generating original K (left hand side coefficients) and P (right hand side
%known values) matrices:
K = zeros(Element_Count, Element_Count);
P = zeros(Element_Count, 1);
for j = 1:M %i + (N+2) * (j-1) for indexing the column values
    for i = 1:(N+2)
        K(i + (N+2) * (j-1), i + (N+2) * (j-1)) = A; %Fills the diagonal with the A values
        if (i == 1)
            K(i + (N+2) * (j-1), i + (N+2) * (j-1) + 1) = 2 * (DY^2);
        elseif (i == N+2)
            K(i + (N+2) * (j-1), i + (N+2) * (j-1) - 1) = 2 * (DY^2);
        else
            K(i + (N+2) * (j-1), i + (N+2) * (j-1) + 1) = DY^2;
            K(i + (N+2) * (j-1), i + (N+2) * (j-1) - 1) = DY^2;
        end
    end
end
