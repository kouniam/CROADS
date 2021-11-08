beta = 2;

i = 1;

for rho=0.1:0.05:1.1
    
    mu = (1/beta)*log(rho);
    
    M(i) = mu;
    
    i = i+1;

end