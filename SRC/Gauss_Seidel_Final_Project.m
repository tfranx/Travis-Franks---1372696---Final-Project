%Travis Franks 1372696 Helmholtz Equation Final Project Gauss Seidel Method
clearvars
clc

%Defining ax, bx, ay, and by:
ax = -pi;
bx = pi;
ay = -pi;
by = pi;

%Setting the value of X_Internal_Nodes, the number of elements along X-domain:
X_Internal_Nodes = input('Enter value of X_Internal_Nodes, the number of internal nodes for the X-domain: ');

%Setting the value of Y_Internal_Nodes, the number of elements along
%Y-domain:
Y_Internal_Nodes = input('Enter value of Y_Internal_Nodes, the number of internal nodes for the Y-domain: ');

%Setting the value of C (capital lambda in problem statement):
C = input('Enter value of C, the given constant for capital lambda: ');

%Setting the value of Es, the acceptable limit of error for system convergence:
Es = 10^-9;

%Defining L, the length of the X and Y domains:
L = 2 * pi();

%Determining DX and DY, the change of X and Y respectively from element to
%element:
DX = L / (X_Internal_Nodes + 1);
DY = L / (Y_Internal_Nodes + 1);

%Creating matrices for X and Y domains:
[X,Y] = meshgrid(ax:DX:bx, ay:DY:by);
X = X';
Y = Y';

%Defining coefficients, A and B, that remain constant throughout matrix operations:
B = (DX^2) * (DY^2);
A = (B * C) - (2 * (DX^2)) - (2 * (DY^2));

%Determining number of elements  on both X and Y domains, so as to only
%perform the calculation one time for optimization:
N = X_Internal_Nodes + 2; %N = number of X domain values
M = Y_Internal_Nodes + 2; %M = number of Y domain values

%Performing setup for Gauss-Seidel Approximation that will solve for U values
%(unknown solution values over X and Y domains):
U = zeros(N,M); %Provides initial guesses and preallocates for optimization

%for i = 1:N
%    for j = 1:M
%        U(i,j) = 1 + X(i,j)^2 + 2 * Y(i,j)^2; %To be used only for method of manufactured solutions
%    end
%end

W = zeros(N,M);

Z = 0; %Z functions as a counter for number of iterations performed during Gauss-Seidel

Error = zeros(N,M-2); %(M-2) rather than just M because the Dirichlet Boundary Conditions cause two rows to have 
%constant values and therefore will have an error of 0 per iteration (optimal to exclude unnecessary repeated calculations) 

Ea = 100; %Provides initial value of Ea, or the relative iterative error, for the Gauss_Seidel Approximation

%Evaluating for F values:
F = zeros(N,M);
for i = 1:N
    for j = 1:M
        %F(i,j) = 6 + C * (1 + X(i,j)^2 + 2 * Y(i,j)^2); %To be used only for Method of Manufactured Solutions
        %F(i,j) = 0; %To be commented out unless simulating Laplace or debugging
        F(i,j) = cos((pi/2) * (2 * ((X(i,j) - ax) / (bx - ax)) + 1)) * sin(pi * ((Y(i,j) - ay) / (by - ay)));
    end
end

%Solving for U values defined by Dirichlet Boundary Conditions that will
%remain constant as Gauss Seidel iterations are performed:
for i = 1:N
    U(i,1) = cos(pi() * DX * (i-1)) * cosh((2 * pi()) - (DX * (i-1))); %Dirichlet BC at Y = -pi
    U(i,M) = ((i-1) * DX)^2 * sin(((i-1) * DX) / 4); %Dirichlet BC at Y = pi
end

%Defining node indexing points to be used for the general form of the
%discretization so as to only perform calculation once for optimization:
NN = N-1;
MM = M-1;

%Defining initial Tcheck value for checkpointing every set period of time of
%the Gauss-Seidel loop:
Tcheck = 60; %checkpoint every 60 seconds

%Defining initial timer variable value to be greater than start timer initial value to start timer loop in Gauss-Seidel loop:
T = 0;
Tcounter = 0; %Stores time values for combined loop times, if loop takes less than a minute to run
Ttotal = 0; %Stores the total time elapsed while performing Gauss-Seidel iterations
save('Variables.mat') %Resets variable values before load checkpoint, so code can be run from either the start or load point
%%
%*****To load variables for checkpointing, run starting from this block*****
load('Variables.mat')

%Performing Gauss Seidel Approximation to solve for U values:
if (A ~= 0)%Set condition for just in case the variable coefficient is equal to zero from a poor choice in nodes along X and Y domains, as it will function as a denominator
    while (Ea > Es)
        tic;
        if (Tcounter >= Tcheck) %Saves at least every 60 seconds or every Gauss-Seidel loop iteration, if the iterations take longer
            Tcounter = 0; %Resets time counter to 0
            save('Variables') %Saves variables for checkpointing for start/restart capability
        end
        
        %Evaluating for general expression (internal nodes):
        for j = 2:MM
            for i = 2:NN
                W(i,j) = U(i,j); %W saves value of U for error calculation
                
                %OPTIMIZES via loop unrolling by evaluating from multiple
                %start points on domain within loop:
                
                %Evaluates starting from bottom left corner of domain:
                U(i,j) = (B * F(i,j) - (DY^2) * U(i-1,j) - (DY^2) * U(i+1,j) - (DX^2) * U(i,j-1) - (DX^2) * U(i,j+1)) / A;
                
                %Evaluates starting from top right corner of domain:
                U(N-i+1,M-j+1) = (B * F(N-i+1, M-j+1) - (DY^2) * U(N-i,M-j+1) - (DY^2) * U(N-i+2,M-j+1) - (DX^2) * U(N-i+1,M-j) - (DX^2) * U(N-i+1,M-j+2)) / A;
                Error(i,j) = abs((U(i,j) - W(i,j)) / U(i,j)); %Computes relative error for this calculation inside this iteration
            end
            
            %Evaluating for U (solution) values defined by Neumann Boundary
            %Conditions that are evaluated by Gauss-Seidel Method:
            W(1,j) = U(1,j); %W saves value of U for error calculation
            U(1,j) = (B * F(1,j) - 2 * (DY^2) * U(2,j) - (DX^2) * U(1,j-1) - (DX^2) * U(1,j+1)) / A;  %cos((pi() / 2) * (0 + 1)) = 0 so Fi,j = 0, Neumann condition for X = -pi
            Error(1,j) = abs((U(1,j) - W(1,j)) / U(1,j)); %Computes relative error for this calculation inside this iteration
            
            W(N,j) = U(N,j); %W saves value of U for error calculation
            U(N,j) = (B * F(N,j) - 2 * (DY^2) * U(NN,j) - (DX^2) * U(N,j-1) - (DX^2) * U(N,j+1)) / A; %cos(3 * pi() / 2) = 0 so Fi,j = 0, Neumann condition for X = pi
            Error(N,j) = abs((U(N,j) - W(N,j)) / U(N,j)); %Computes relative error for this calculation inside this iteration
            
            %Averaging contributions from boundary conditions at corner
            %points to minimize large data spikes at the corners:
            if (j==2)
                U(1,1) = (U(1,2) + U(2,1)) / 2;
                U(N,1) = (U(NN,1) + U(N,2)) / 2;
                U(1,M) = (U(2,M) + U(1,MM)) / 2;
                U(N,M) = (U(NN,M) + U(N,MM)) / 2;
            end
        end
        
        Ea = max(max(Error)); %Computes the maximum relative error of the Gauss-Seidel iteration
        
        Z = Z + 1; %Counts the number of loop iterations
        
        T = toc; %Determines how long the Gauss-Seidel iteration took
        
        Tcounter = Tcounter + T; %Counts how much time has elapsed since the last variable checkpoint save
        
        Ttotal = Ttotal + T; %Stores the total time elapsed while performing Gauss-Seidel iterations
    end
else
    disp('Select a different number of nodes for X or Y domain or change the value of C, the given constant for capital lambda.')
end

%Performing statisical analysis of U solution matrix for Grid Independence
%Study:
Grid = mean(mean(U.^2));

save('Variables.mat')
%%
%Plotting visualizations for ease of interpretation of results:
load('Variables.mat') %Provides the option of loading the variables directly onto the plot section, if a failure occurs while plotting

%Displaying the total time elapsed during the Gauss-Seidel approximation
%iterations, Ttotal (in seconds):
disp('Total time elapsed, in seconds = ')
disp(Ttotal)

%Plotting surface plot of U, the solution matrix for the Helmholtz
%equation:
figure
subplot(1,2,1)
mesh(X, Y, U) %Produces surface plot of solution matrix for the Helmholtz equation:
xlabel('X axis'), ylabel('Y axis'), title('Surface Plot of U for Lambda = 0 and F = 0'), colorbar

%Plotting contour lines of U, the solution matrix for 2-D interpretation of
%results obtained:
subplot(1,2,2)
[Matrix, Object] = contourf(X, Y, U, 20); %Plots the contour of the solution matrix
xlabel('X axis'), ylabel('Y axis'), title('Contour Plot of U for Lambda = 0 and F = 0'), colorbar
clabel(Matrix, Object) %Labels the peak values for all of the contour lines
