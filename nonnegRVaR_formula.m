function nonnegRVaR_formula()
%Theorem 4

%Choose the values for mu, s, alpha, and beta
mu=0.00675;
s=0.01824;
alpha=0.2;
beta=0.6;
[res,zone]=supRVaR(alpha,beta,s,mu);

if zone <=10
disp(['For alpha = ' num2str(alpha) ' and beta = ' num2str(beta) ', Zone = A' num2str(zone)   ' and MaxRVaR = ' num2str(res)]) 
elseif  zone == 11
disp(['For alpha = ' num2str(alpha) ' and beta = ' num2str(beta) ', Zone = B1'  ' and MaxRVaR = ' num2str(res)])     
elseif zone == 12
disp(['For alpha = ' num2str(alpha) ' and beta = ' num2str(beta) ', Zone = B2'   ' and MaxRVaR = ' num2str(res)])     
elseif zone == 13
disp(['For alpha = ' num2str(alpha) ' and beta = ' num2str(beta) ', Zone = B3'   ' and MaxRVaR = ' num2str(res)])     
end

end





function [res,zone]=supRVaR(alpha,beta,s,mu)
b1= alpha*(3*alpha+2-sqrt((3*alpha-2)^2+12*(1-beta)))/(2*(2*alpha+beta-1));
sstar=mu*sqrt((alpha/3)*(2*sqrt(alpha^2+4*(1-beta))-alpha)/(alpha^2+4*(1-beta)));
b3=(4*(1-beta)+sqrt(alpha^2+4*(1-beta))*(2*alpha-sqrt(2)*sqrt(alpha^2+2*(1-beta)-alpha*sqrt(alpha^2+4*(1-beta)))))/(6*(1-beta)+alpha^2+alpha*sqrt(alpha^2+4*(1-beta)));
Qb3sstar= mu+sstar/(beta-alpha)*sqrt(3)*(b3^2*(beta-alpha-1)+2*b3*alpha-alpha^2)/sqrt(b3^3*(4-3*b3));
L= (Qb3sstar-mu)*sqrt((2-alpha-beta)/(alpha+beta-10/9));
K= mu*sqrt((3*(alpha+beta)-2-4*sqrt(1-(2-alpha-beta)*Qb3sstar/mu))/(3*(2-alpha-beta)));
if (alpha>=1/2) &&(alpha<5/6)
    fun = @(beta) auxzero(beta,alpha);    % function of x alone
   betastar= fzero(fun,[alpha,1]); % we solve for beta within alpha,1
end


if (2*alpha+beta==1) && (alpha>0) && (alpha<=1/3) && (s>=sstar) 
    % set B1
    res=Qb3sstar;
    zone=13;
    
elseif (2*alpha+beta==1) && (alpha>0) && (alpha<=1/3) && (s>=mu*sqrt(alpha*(3*alpha+8))/(3*alpha+4)) && (s<sstar) 
    % set B2
    b2=2/3*(mu^2+3*s^2-mu*sqrt(mu^2-3*s^2))/(mu^2+s^2);
    res=mu+s/(beta-alpha)*sqrt(3)*(b2^2*(beta-alpha-1)+2*b2*alpha-alpha^2)/sqrt(b2^3*(4-3*b2)); 
    zone=12;

elseif (2*alpha+beta==1) && (alpha>0) && (alpha<=1/3) && (s>=0) && (s<mu*sqrt(alpha*(3*alpha+8))/(3*alpha+4)) 
    % set B3
    res=mu+s/3*sqrt(alpha*(3*alpha+8));
    zone=11;

elseif (alpha>0) && (alpha<2/3) && (beta>=max(alpha,1-alpha)) && (beta<=min(4/3-alpha,1)) && (s>=mu*sqrt((alpha+beta-2/3)/(2-(alpha+beta))))
      % set A10   
    res=max(mu/(2-(alpha+beta)),Qb3sstar);
    zone=10; 
    
elseif (alpha>0) && (alpha<2/3) && (beta>=max(alpha,1-alpha)) && (beta<=min(4/3-alpha,1)) && (s>=mu/sqrt(3))&& (s<mu*sqrt((alpha+beta-2/3)/(2-(alpha+beta))))
      % set A9  
   res=max(3/16*mu*(-3*(s/mu)^4*(2-(alpha+beta))+2*(s/mu)^2*(3*(alpha+beta)-2)+3*(alpha+beta)+2),Qb3sstar);
   zone=9;
   
elseif (alpha>0) && (alpha<2/3) && (beta>=alpha) && (beta<=min(4/3-alpha,1)) && (s>=mu*sqrt(b1*(4/3-b1))/(2-b1))&& (s<sstar)
      % set A8   
   b2=2/3*(mu^2+3*s^2-mu*sqrt(mu^2-3*s^2))/(mu^2+s^2);
   res=mu+s/(beta-alpha)*sqrt(3)*(b2^2*(beta-alpha-1)+2*b2*alpha-alpha^2)/sqrt(b2^3*(4-3*b2)); 
   zone=8;
   
elseif (alpha>=1/2) && (alpha<5/6) && (beta>=max(betastar,4/3-alpha)) && (beta<=1) && (s>=0) && (s<mu*sqrt(b1*(4/3-b1))/(2-b1))
      % set A7    
    res=max(mu+s*sqrt(8/(9*(2-(alpha+beta)))-1),mu+s/(beta-alpha)*sqrt(3)*(b1^2*(beta-alpha-1)+2*b1*alpha-alpha^2)/sqrt(b1^3*(4-3*b1)));
    zone=7;
    
 elseif (alpha>=1/3) && (alpha<5/6) && (beta>=max(alpha,4/3-alpha)) && (beta<=1) && (s>=mu*sqrt(b1*(4/3-b1))/(2-b1)) && (s<sstar)
    % set A6
    b2=2/3*(mu^2+3*s^2-mu*sqrt(mu^2-3*s^2))/(mu^2+s^2);
    res=max(mu+s*sqrt(8/(9*(2-(alpha+beta)))-1),mu+s/(beta-alpha)*sqrt(3)*(b2^2*(beta-alpha-1)+2*b2*alpha-alpha^2)/sqrt(b2^3*(4-3*b2))); 
    zone=6;    
      
elseif (alpha>0) && (alpha<1/2) && (beta>=alpha) && (beta<=1) && (s>=0)&& (s<mu*sqrt(b1*(4/3-b1))/(2-b1))
    % third part of set A5  
    res=mu+s/(beta-alpha)*sqrt(3)*(b1^2*(beta-alpha-1)+2*b1*alpha-alpha^2)/sqrt(b1^3*(4-3*b1));
    zone=5;
    
elseif (alpha>=1/2) && (alpha<2/3) && (beta>=alpha) && (beta<max(betastar,4/3-alpha)) && (s>=0)&& (s<mu*sqrt(b1*(4/3-b1))/(2-b1))
    % second part of set A5 
    res=mu+s/(beta-alpha)*sqrt(3)*(b1^2*(beta-alpha-1)+2*b1*alpha-alpha^2)/sqrt(b1^3*(4-3*b1));
    zone=5;

elseif (alpha>=2/3) && (alpha<5/6) && (beta>=alpha) && (beta<betastar) && (s>=0) && (s<mu*sqrt(b1*(4/3-b1))/(2-b1))
    % first part of set A5
    res=mu+s/(beta-alpha)*sqrt(3)*(b1^2*(beta-alpha-1)+2*b1*alpha-alpha^2)/sqrt(b1^3*(4-3*b1));
    zone=5;
    
elseif (alpha>0) && (alpha<1/2) && (beta>=alpha) && (beta<1-alpha) && (s>=mu/sqrt(3))
    % fourth part of set A4  
  res=Qb3sstar; 
  zone=4;
  
 elseif (alpha>0) && (alpha<2/3) && (beta>=alpha) && (beta<=min(1,4/3-alpha)) && (s>=sstar)&& (s<mu/sqrt(3))
    % third part of set A4   
   res=Qb3sstar; 
  zone=4;
  
 elseif (alpha>=1/3) && (alpha<5/7) && (beta>=max(4/3-alpha,alpha)) && (beta<=1) && (s>=mu*sqrt((alpha+beta-10/9)/(2-(alpha+beta)))) && (s<K)
    % second part of set A4 
  res=Qb3sstar;
  zone=4;
 
  elseif (alpha>=1/3) && (alpha<5/6) && (beta>=max(4/3-alpha,alpha)) && (beta<=1) && (s>=sstar) && (s<=min(L,mu*sqrt((alpha+beta-10/9)/(2-(alpha+beta)))))
    % first part of set A4
  res=Qb3sstar;
  zone=4;
  
 elseif (alpha>=1/3) && (alpha<5/6) && (beta>=max(4/3-alpha,alpha)) && (beta<=1) && (s>=max(sstar,L)) && (s<mu*sqrt((alpha+beta-10/9)/(2-(alpha+beta))))
    % 2nd part of set A3
    res=mu+s*sqrt(8/(9*(2-(alpha+beta)))-1);
    zone=3;
    
elseif (alpha>=5/6) && (alpha<1) && (beta>=alpha) && (beta<=1) && (s>=0) &&  (s<mu*sqrt((alpha+beta-10/9)/(2-(alpha+beta))))
    % 1st part of set A3
    res=mu+s*sqrt(8/(9*(2-(alpha+beta)))-1);
    zone=3; 
    
elseif (alpha>=1/3) && (alpha<5/7) && (beta>=max(4/3-alpha,alpha)) && (beta<=1)&& (s>=max(K,mu*sqrt((alpha+beta-10/9)/(2-(alpha+beta))))) && (s<mu*sqrt((alpha+beta-2/3)/(2-(alpha+beta))))
    % second part of set A2
    res=3/16*mu*(-3*(s/mu)^4*(2-(alpha+beta))+2*(s/mu)^2*(3*(alpha+beta)-2)+3*(alpha+beta)+2);
    zone=2;
   
elseif (alpha>=5/7) && (alpha<1) && (beta>=alpha) && (beta<=1) && (s>=mu*sqrt((alpha+beta-10/9)/(2-(alpha+beta)))) && (s<mu*sqrt((alpha+beta-2/3)/(2-(alpha+beta))))
    % first part of set A2
    res=3/16*mu*(-3*(s/mu)^4*(2-(alpha+beta))+2*(s/mu)^2*(3*(alpha+beta)-2)+3*(alpha+beta)+2);
    zone=2;
   
elseif (alpha>=1/3) && (alpha<1) && (beta>=max(4/3-alpha,alpha)) && (beta<=1) && (s>=mu*sqrt((alpha+beta-2/3)/(2-(alpha+beta))))
    % set A1
    res=mu/(2-(alpha+beta));
    zone=1;
    
elseif beta<=alpha
    res=0;
    zone=0;
    disp(['error with mu=' num2str(mu) ', s=' num2str(s) ', alpha=' num2str(alpha) ',beta=' num2str(beta) ])
else
    res=0;
    zone=14;
    disp(['error with mu=' num2str(mu) ', s=' num2str(s) ', alpha=' num2str(alpha) ',beta=' num2str(beta) ])
end


end

function res=auxzero(beta,alpha)
    
    res=14*beta.^4/9-(4/3*alpha+95/27)*beta.^3-(2*alpha^2-5*alpha-2)*(beta.^2)+alpha*(alpha-4)*beta-alpha^3+2*alpha^2;

end