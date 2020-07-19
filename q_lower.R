library(pracma)

fdrlasso =function (tpp, delta, epsi){
#--------------------------------------------------------------------------
# This function essentially adapts from the MATLAB code from https://github.com/wjsu/fdrlasso.
# It calculates the lower bound of FDP of Lasso given TPP (true
# positive proportion), DELTA = n/p (shape of the design matrix, or
# subsampling rate), and EPSI = k/p (sparsity ratio).
# All TPP, DELTA, and EPSI are between 0 and 1; if the
# pair (DELTA, EPSI) is above the Donoho-Tanner phase transition, TPP
# should be no larger than u^\star = powermax(DELTA, EPSI)
#--------------------------------------------------------------------------

if (tpp == 0){
  q = 0;
  return (q)
}

## make stepsize smaller for higher accuracy
stepsize = 0.1;
tmax = max(10, sqrt(delta/epsi/tpp) + 1);
tmin = tmax - stepsize;

while (tmin > 0){
  if (lsandwich(tmin, tpp, delta, epsi) < rsandwich(tmin, tpp)){
    break
  }
  
  tmax = tmin;
  tmin = tmax - stepsize;
}

if (tmin <= 0){
  stepsize = stepsize/100;
  tmax = max(10, sqrt(delta/epsi/tpp) + 1);
  tmin = tmax - stepsize;
  while (tmin > 0){
    if (lsandwich(tmin, tpp, delta, epsi) < rsandwich(tmin, tpp)){
      break
    }
     
    tmax = tmin;
    tmin = tmax - stepsize;
  }
}
  
diff = tmax - tmin;
while (diff > 1e-6){
  tmid = 0.5*tmax + 0.5*tmin;
  if (lsandwich(tmid, tpp, delta, epsi) > rsandwich(tmid, tpp)){
    tmax = tmid;
  }
  else   
  { 
    tmin = tmid;
  }
  diff = tmax - tmin;
}

t = (tmax + tmin)/2;

q = 2*(1-epsi)*pnorm(-t)/(2*(1-epsi)*pnorm(-t) + epsi*tpp);

return (q);
}




####################################
lsandwich=function(t, tpp, delta, epsi){
Lnume = (1-epsi)*(2*(1+t^2)*pnorm(-t) - 2*t*dnorm(t)) + epsi*(1+t^2) - delta;
Ldeno = epsi*((1+t^2)*(1-2*pnorm(-t)) + 2*t*dnorm(t));
L = Lnume/Ldeno;
return (L);
}
####################################
rsandwich=function (t, tpp){
R = (1 - tpp)/(1 - 2*pnorm(-t));
return (R);
}


####################################
## The function of solving the maximum power for delta < 1 and epsilon > epsilon_phase
powermax=function(delta, epsilon){
if (delta >= 1){
  power = 1;
  return;
}
epsilon_star = epsilonDT(delta);
if (epsilon <= epsilon_star){
  power = 1;
  return;
}
power = (epsilon - epsilon_star)*(delta - epsilon_star)/epsilon/(1 - epsilon_star) + epsilon_star/epsilon;
return (power)
}

############# 
## The function of solving the critical epsilon
epsilonDT=function (delta){
minus_f = function(x) -(1+2/delta*x*dnorm(x) - 2/delta*(1+x^2)*pnorm(-x))/(1+x^2-2*(1+x^2)*pnorm(-x)+2*x*dnorm(x))*delta;
alpha_phase = fminbnd(minus_f, 0, 8)$xmin;
epsi = -minus_f(alpha_phase);
return (epsi)
}

