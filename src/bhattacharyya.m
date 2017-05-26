function bdist = bhattacharyya(mu1,mu2,sigma1,sigma2)
sigma = (sigma1 +sigma2) / 2;
bdist = 1/8*(mu1-mu2)/sigma*(mu1-mu2)' + 1/2*log(det(sigma)/sqrt(det(sigma1)*det(sigma2)));