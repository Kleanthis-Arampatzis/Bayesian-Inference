% Group Lasso with Gibbs sampling
function [Final_Beta_mean]=q2_BayesGroupLasso(y,x)
    
    [n,p]=size(x);
    %Center the data
    y_til=y-mean(y);
    x_til=x-mean(x);
    
    %Burn_in
    burn_in=100;
    niter=1000;
    nsim=burn_in+niter +1;
    keep = (burn_in + 2);

    %Grouping
    %Setting the size of each group 
    %the groups function returns the number of
    %groups and a vector that demonstrates in 
    %which group each covariate belongs to
    size_g=3;
    [group_size, mk, K]=groups(x_til,size_g);

    %Priors
    %prior sigma2
    sigma2 = 1/3;
    %prior lambda2
    r = 1;
    delta = 1.78;
    lambda2=gamrnd(r,delta); 
    %prior tau2

    invtau2=1./gamrnd((mk+1)/2,lambda2/2,K,1);
    %prior beta
    beta=zeros(p,1);
    for i=1:K
        prior_beta_mu=zeros(mk(i),1);
        prior_beta_var = sigma2*invtau2(i)*eye(mk(i));
        %beta is px1
        beta(group_size==i) = mvnrnd(prior_beta_mu,prior_beta_var)';
    end

    %Initialization of posteriors
    postBeta =zeros(p,nsim);
    postSigma2=zeros(1,nsim);
    postinvTau2=zeros(K,nsim);
    postLambda2=zeros(1,nsim);
    postBeta(:,1) = beta;
    postSigma2(1,1) = sigma2;
    postinvTau2(:,1) = invtau2;
    postLambda2(1,1) = lambda2;
    

    for iter = 2:nsim
         %update beta
         for i=1:K
             Xk=x_til(:,group_size==i);
             Ak_inv=inv(Xk'*Xk + (invtau2(i)).*eye(mk(i)));
             sxkbgk=0.5.*(x_til(:,group_size~=i)*beta(group_size~=i));%check
             mu = Ak_inv*Xk'*(y_til-sxkbgk);
             varcov = sigma2*Ak_inv;
             beta(group_size==i) = mvnrnd(mu,varcov)';
             postBeta(group_size==i,iter) = beta(group_size==i);
         end

         %Sample sigma2 
         shape = (n+p-1)/2;
         %psi is the greek letter ø
         psi = (y_til - x_til*beta)'*(y_til - x_til*beta);
         bgk2=bgk(beta,K,group_size);
         scale = (psi + sum(invtau2.*bgk2)) / 2; 
         sigma2 = 1 / gamrnd(shape,1/scale); 
         postSigma2(:,iter) = sigma2;

         %Sample tau2
         mu = sqrt((lambda2 * sigma2) ./ bgk2);
         invtau2 = random('InverseGaussian', mu , lambda2);
         postinvTau2(:,iter) = invtau2;

         %Update lambda
         shape = (p+K)/2+r;
         scale = 1/((sum(1./invtau2)/2)+delta );
         lambda2 = gamrnd(shape, scale);
         postLambda2(:,iter) = lambda2;

    end

    finalBeta=postBeta(:,keep:end);
    Final_Beta_mean=mean(finalBeta,2);
end