;Define a vector of weights.
weights=1.0/e4p5ex_45
 
;Provide an initial guess of the function’s parameters.
a=[10.^result_45[0],result_45[1]]
 
;Compute the parameters.
yfit = CURVEFIT(n_45,l4p5ex_45,weights,A,SIGMA,FUNCTION_NAME='powfunct')
 
;Print the parameters returned in A.
PRINT, 'Function parameters: ', A