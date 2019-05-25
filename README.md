# Weibull
fitting the weibull distribution

I used the Excel description by Charles Zaiontz in his Excel blog

SimpleFit is this

https://www.real-statistics.com/distribution-fitting/method-of-moments/method-of-moments-weibull/

which the Excel implementation of using the Mean & variance to estimate the parameters as described in this paper

SIMPLIFIED  METHOD-OF-MOMENTS  ESTIMATION  FOR  THE  WEIBULL  DISTRIBUTION  

OSCAR   GARCIA

https://www.scionresearch.com/__data/assets/pdf_file/0010/36874/NZJFS1131981GARCIA304-306.pdf

Then Fit is this refinement, which Charles uses Solver to calculat the MLE

https://www.real-statistics.com/distribution-fitting/distribution-fitting-via-maximum-likelihood/fitting-weibull-parameters-mle/


Whereas I use a NLP Solver in JuMP
