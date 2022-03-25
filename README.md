# homework2
Financial Computation &amp; Simulation


1.	Write a function in R to convert a decimal (base 10) to binary number.  Please extend the code we started in class to handle fractional numbers (e.g., 200.65).  [NOTE: your code should truncate the precision of the fractional part of the number].  See section 3.1.1 for help.

2.	Using R, confirm that the LU factorization out-performs R’s Gaussian elimination (x = A\b) procedure in terms of speed for the following system of equations:

A = [4 2 3; 5 -8 1; 4 7 -9]
b = [1 4 0]’

[HINT: Compare the speed of each methodology for a number of iterations].   

3.	Write an R function to solve the following systems of equations using the Jacobi iterative method discussed in class.  The code should check to ensure the system will converge and if not, it should report a warning message.  
A1 = [3 1 1 0; 1 5 -1 2; 1 0 3 1; 0 1 1 4];
A2 = [2.5 1 1 0; 1 4.1 -1 2; 1 0 2.1 1; 0 1 1 2.1];
A3 = [2 1 1 0; 1 3.5 -1 2; 1 0 2.1 1; 0 1 1 2.1];
b1 = b2 = b3 = [1 4 -2 1]';


4.	Write an R function to solve the above systems of equations using the Gauss-Seidel iterative method discussed in class.  The code should check to ensure the system will converge and if not, it should report an error message.
A1 = [3 1 1 0; 1 5 -1 2; 1 0 3 1; 0 1 1 4];
A2 = [2.5 1 1 0; 1 4.1 -1 2; 1 0 2.1 1; 0 1 1 2.1];
A3 = [2 1 1 0; 1 3.5 -1 2; 1 0 2.1 1; 0 1 1 2.1];
b1 = b2 = b3 = [1 4 -2 1]';


