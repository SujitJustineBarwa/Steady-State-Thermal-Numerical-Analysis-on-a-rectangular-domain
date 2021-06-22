# Steady State Thermal Numerical Analysis on a rectangular domain
 An explicit time stepping method is used find the steady state thermal contour on a rectangular domain.

Special Thanks to Dr. Ashwani Assam (IIT,PATNA)
 
# How to use:
1. Set L,W which are the length and width of the rectangular domain.\
2. Set the step size the dx.\
3. Set the T1(the temperature for left right and bottom edges) and T2(the temperature for top edge).\
4. Set temperature in line 34(Recomended : avg(T1,T2)).\
5. Set a value of n(the step size for analytical temperature calculation (Recomended : n = min(L,W)/dx).\
6. Hit the Run button

# Interpreting the Results
1.The plot 1 is the analytical temp plot.\
2.The plot 2 is the solution we get from iterative calculation.\
3.The plot 3 is the center point temperature along with iteration number(in current setting it becomes constant at 120 deg approx.\
4.Residue in the command window is the heat input minus heat output. 
