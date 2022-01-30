function nonnegVaR_formula()
% Theorem 3
%Choose the values for mu, s, and alpha
mu=0.00675;
s=0.01824;
alpha=0.6;
[res,zone]=supVaR(alpha,s,mu);

disp(['For alpha = ' num2str(alpha,2) ', Zone = A' num2str(zone) ' and MaxVaR = ' num2str(res)]) 


end


function [res,zone]=supVaR(alpha,s,mu)
if (alpha>=2/3) && (alpha<1) && (s>=mu*sqrt((alpha-1/3)/(1-alpha)))
   
    % set A1
    res=mu/(2*(1-alpha));
    zone=1;
elseif (alpha>=2/3) && (alpha<5/7) && (s>=mu*sqrt((3*alpha-1-2*sqrt((3*alpha-2)/(2-alpha)))/(3*(1-alpha)))) && (s<mu*sqrt((alpha-1/3)/(1-alpha)))
    % first part of set A2
    res=3/8*mu*(-3*(s/mu)^4*(1-alpha)+2*(s/mu)^2*(3*alpha-1)+3*alpha+1);
    zone=2;
elseif (alpha>=5/7) && (alpha<1) && (s>=mu*sqrt((alpha-5/9)/(1-alpha))) && (s<mu*sqrt((alpha-1/3)/(1-alpha)))
    % second part of set A2
    res=3/8*mu*(-3*(s/mu)^4*(1-alpha)+2*(s/mu)^2*(3*alpha-1)+3*alpha+1);
    zone=2;
elseif (alpha>=5/6) && (alpha<1) && (s>=0) && (s<mu*sqrt((alpha-5/9)/(1-alpha)))
    % first part of set A3
    res=mu+s*sqrt(4/(9*(1-alpha))-1);
    zone=3;
elseif (alpha>=5/7) && (alpha<5/6) && (s>=mu*alpha/(2-alpha)*sqrt(9*(1-alpha)/(9*alpha-5))) && (s<mu*sqrt((alpha-5/9)/(1-alpha)))
    % second part of set A3
    res=mu+s*sqrt(4/(9*(1-alpha))-1);
    zone=3;
elseif (alpha>0) && (alpha<2/3) && (s>=mu/(2-alpha)*sqrt(alpha*(4/3-alpha)))
    % first part of set A4
    res=2*mu/(2-alpha);
    zone=4;
elseif (alpha>=2/3) && (alpha<5/7) && (s>=mu/(2-alpha)*sqrt(alpha*(4/3-alpha))) && (s<mu*sqrt((3*alpha-1-2*sqrt((3*alpha-2)/(2-alpha)))/(3*(1-alpha))))
    % second part of set A4
    res=2*mu/(2-alpha);
    zone=4;
elseif (alpha>=5/7) && (alpha<5/6) && (s>=mu/(2-alpha)*sqrt(alpha*(4/3-alpha))) && (s<mu*alpha/(2-alpha)*sqrt((9*(1-alpha)/(9*alpha-5))))
    % second part of set A4
    res=2*mu/(2-alpha);
    zone=4;
elseif (alpha>0) && (alpha<5/6) && (s>=0) && (s<mu*sqrt(alpha*(4/3-alpha))/(2-alpha))
    % set A5
    res=mu+s*sqrt(3*alpha/(4-3*alpha));
    zone=5;
else
    
    disp(['error with mu=' num2str(mu) ', s=' num2str(s) ', alpha=' num2str(alpha) ])
end


end
