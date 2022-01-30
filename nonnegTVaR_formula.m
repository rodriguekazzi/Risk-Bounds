function nonnegTVaR_formula()
% Corollary 5
%Choose the values for mu, s, and alpha

mu=0.00675;
s=0.01824;
alpha=0.95;
[res,zone]=supTVaR(alpha,s,mu);

disp(['For alpha = ' num2str(alpha) ', Zone = A' num2str(zone) ' and MaxTVaR = ' num2str(res)]) 

end

function [res,zone]=supTVaR(alpha,s,mu)

if (alpha>=0) && (alpha<1) && (s>=mu*sqrt((alpha+1/3)/(1-alpha)))
    % set A1
    res=mu/(1-alpha);
    zone=1;

elseif (alpha>=1/3) && (alpha<1) && (s>=mu*sqrt((alpha-1/9)/(1-alpha))) && (s<mu*sqrt((alpha+1/3)/(1-alpha)))
    % first part of set A2
    res=3/16*mu*(-3*(s/mu)^4*(1-alpha)+2*(s/mu)^2*(3*alpha+1)+3*alpha+5);
    zone=2;
   
elseif (alpha>0) && (alpha<1/3) && (s>=mu/sqrt(3)) && (s<=mu*sqrt((alpha+1/3)/(1-alpha)))
     % second part of set A2 
   res=3/16*mu*(-3*(s/mu)^4*(1-alpha)+2*(s/mu)^2*(3*alpha+1)+3*alpha+5);
   zone=2;  
 
 elseif (alpha>=1/2) && (alpha<1) && (s>=0) && (s<mu*sqrt((alpha-1/9)/(1-alpha)))
     % first part of set A3 
    res=mu+s*sqrt(8/(9*(1-alpha))-1);
    zone=3;
   
 elseif (alpha>=1/3) && (alpha<1/2)  && (s>=mu/sqrt(3)) && (s<mu*sqrt((alpha-1/9)/(1-alpha)))
    % 2nd part of set A3
    res=mu+s*sqrt(8/(9*(1-alpha))-1);
    zone=3;  
    
 elseif (alpha>=1/3) && (alpha<1/2) && (s>=mu*sqrt(alpha*(8-9*alpha))/(4-3*alpha)) && (s<mu/sqrt(3))
    %  set A4
    b2=2/3*(mu^2+3*s^2-mu*sqrt(mu^2-3*s^2))/(mu^2+s^2);
   res=max(mu+s*sqrt(8/(9*(1-alpha))-1),mu+s/(1-alpha)*sqrt(3)*(b2^2*(-alpha)+2*b2*alpha-alpha^2)/sqrt(b2^3*(4-3*b2))); 
    zone=4; 
    
 elseif (alpha>0) && (alpha<1/3) && (s>=mu*sqrt(alpha*(8-9*alpha))/(4-3*alpha))&& (s<mu/sqrt(3))
      % set A5
  b2=2/3*(mu^2+3*s^2-mu*sqrt(mu^2-3*s^2))/(mu^2+s^2);
   res=mu+s/(1-alpha)*sqrt(3)*(b2^2*(-alpha)+2*b2*alpha-alpha^2)/sqrt(b2^3*(4-3*b2)); 
   zone=5;   
    
 elseif (alpha>0) && (alpha<1/2) && (s>=0)&& (s<mu*sqrt(alpha*(8-9*alpha))/(4-3*alpha))
      %  set A6
    res=mu+s*sqrt(alpha*(8-9*alpha))/(3*(1-alpha));
    zone=6; 

else
    res=0;
    zone=7;
    disp(['error with mu=' num2str(mu) ', s=' num2str(s) ', alpha=' num2str(alpha) ])
end


end

