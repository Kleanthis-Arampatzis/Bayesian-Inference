% Lasso with Gibbs sampling

function [Final_Beta_mean]=q1_BayesLasso(y,x)
    %Center the data
    y_til=y-mean(y);
    x_til=x-mean(x);
    [n,p]=size(x);
    %Burn_in
    burn_in=100;
    niter=1000;
    nsim=burn_in+niter +1;
    keep = (burn_in + 2);

    %Time saving operations
    xtx = x_til'*x_til;
    xty = x_til.'*y;
    %prior sigma2
    sigma2 = 1/3;
    %prior lambda2
    r = 1;
    delta = 1.78;
    lambda2=gamrnd(r,delta); %proposed by slides
    
    %prior tau2
    invtau2=exprnd(lambda2/2,p,1);
    %prior beta
    prior_Dt=invtau2.*eye(p);
    prior_beta_mu=zeros(p,1);
    prior_beta_var = sigma2.*prior_Dt;
    %beta is px1
    beta = mvnrnd(prior_beta_mu,prior_beta_var)';
    
    postBeta =zeros(p,nsim);
    postSigma2=zeros(1,nsim);
    postinvTau2=zeros(p,nsim);
    postLambda2=zeros(1,nsim);
    postBeta(:,1) = beta;
    postSigma2(:,1) = sigma2;
    postinvTau2(:,1) = invtau2;
    postLambda2(:,1) = lambda2;
    index = linspace(1,p,p);

    for iter = 2:nsim
         %update beta
         invDt = diag(invtau2);
         invV = (xtx + invDt) \ eye(p);
         mu = invV * xty;
         varcov = sigma2 * invV;
         beta = mvnrnd(mu,varcov)';
         postBeta(:,iter) = beta;

         %Sample sigma2 
         shape = (n+p-1)/2;
         %psi is the greek letter ø
         psi = (y_til - x_til*beta)'*(y_til - x_til*beta);
         scale = (psi + beta'*invDt*beta) / 2;
         sigma2 = 1 / gamrnd(shape,1/scale); 
         postSigma2(:,iter) = sigma2;

         %Sample tau2
         mu = sqrt((lambda2 * sigma2) ./ (beta.^2));
         invtau2(index) = random('InverseGaussian', mu(index) , lambda2);
         postinvTau2(:,iter) = invtau2;

         %Update lambda
         shape = r + p;
         scale = 1/(delta + (sum(1./invtau2)/2)); 
         lambda2 = gamrnd(shape, scale);
         postLambda2(:,iter) = lambda2;

    end

    finalBeta=postBeta(:,keep:end);
    Final_Beta_mean=mean(finalBeta,2);
end
