DefensePresentation.Rproj
========================================================
author: Kirk Gosik
date: Sys.Date()
autosize: true

First Slide
========================================================

For more details on authoring R presentations please visit <https://support.rstudio.com/hc/en-us/articles/200486468>.

- Bullet 1
- Bullet 2
- Bullet 3

Background
========================================================
type: section


Background 1
========================================================

Background 2
========================================================

Background 3
========================================================

Motivation
========================================================
type: section

Motivation 1
========================================================


Motivation 1
========================================================


Model
========================================================
type: section


Model Properies
========================================================

Model Properies
========================================================

Model Theoretical Results
========================================================
$$ BIC_1 = n*log(RSS/n) + k*log(n) $$  
$$ BIC_2 = n*log(RSS/n) + k*(log(n) + 2*log(d^{\star})) $$

derived BIC2 by controlling the false discovery rate (FDR) and showed that it is selection consistent if $$d = O(n\xi) for some \xi > 0$$

Wang (2009) showed its selection consistency for FS under ultra-high dimensional setup d = O(exp(n<U+03BE>)).

Simulation Results
========================================================
type: section

Simulation Results
========================================================

Simulation Results
========================================================


Application 1
========================================================
type: section

Application 1 results
========================================================


Application 2 results
========================================================


Application 2
========================================================


Application 2 as Motivation
========================================================

Use Mei Tree results and data to set up motivation for functional response variable and using the iform procedure.  How can we incorporate all information from the response variable into one model?


Background of Functional Reponse
========================================================
type: section

 - Regression as Linear Combination of Basis Functions
 - Legendre Polynomials to model genetic effect of each SNP over time
 - Generalized Least Squares for correlated errors (AR1 model assumed)
 
In these sections talk about Treating linear regression as basis functions, legendre polynomials and generalized least squares. 

Linear regression when there is a certain degree of correlation between the residuals in a regression model. In these cases, ordinary least squares and weighted least squares can be statistically inefficient, or even give misleading inferences



Generalized Least Squares
========================================================

[GLS](https://en.wikipedia.org/wiki/Generalized_least_squares)

The GLS estimator is unbiased, consistent, efficient, and asymptotically normal:


$\sqrt n $

AR(1) model assumed for our purposes.  


Legendre Polynomials
========================================================

 - Orthogonal and therefore variables will not be correlated
 - variety of fits up to the order of the researcher's choosing
 - easily implemented as basis functions 
 -- saves on computational cost and time over using splines or other basis functions learned from the data
 

$$\begin{equation}
\begin{split}
P_n(x) & = \frac{1}{2^n}\sum_{k=0}^{n}{{n}\choose{k}}^2(x-1)^{n-k}(x+1)^k \\
& = \sum_{k=0}^{n}{{n}\choose{k}}{{-n-1}\choose{k}}{\left(\frac{1-x}{2}\right)}^k \\
& = 2^{-n}\sum_{k=0}^{n} x^k {{n}\choose{k}}{{\frac{n+k+1}{2}}\choose{k}} \\
\end{split}
(\#eq:leg-eq)
\end{equation}$$


Legendre Polynomials
========================================================

show toy example of how the polynomials could model the genetic effect well


Generalized Least Squares
========================================================

We have assumed that $$var(\epsilon) = \sigma^2I$$ when the resposne is statit but if we have correlated errors like when we have repeated measurements over time $$var(\epsilon) = \sigma^2\Sigma$$ where sigma^2 is unknown but Sigma is known.  We can use Generalized least squares, Instead of minimizing the ordinary least squares estimate we have to find the arg minimum of 
$$ (y - X\beta)^T\Sigma^{-1}(y - X\beta) $$
which is solved by 
$$ \hat\beta = (X^T\Sigma^{-1}X)^{-1}X^T\Sigma^{-1}y $$
since we can write $$\Sigma = SS^T$$, where S is a triangualar matrix using the Choleski Decomposition, we have
$$ (y - X\beta)^T(S^{-T}S^{-1})(y - X\beta) = (S^{-1}y - S^{-1}X\beta)^T(S^{-1}y - S^{-1}X\beta)$$

So GLS is like regressing $S^{-1}X$ on $S^{-1}y$. Furthermore
$$y = X\beta + \epsilon$$
$$S^{-1}y = S^{-1}X\beta + S^{-1}\epsilon$$
$$y' = X'\beta + \epsilon'$$


So we have a new regression equation $y' = X'\beta + \epsilon'$ where if we examine the variance of the new errors, $\epsilon'$
we find

$$var(\epsilon') = var(S^{-1}\epsilon) = S^{-1}(var(\epsilon))S^{-T} = S^{-1}\sigma^2SS^TS^{-T} = \sigma^2I$$

So the new variables $y'$ and $X'$ are related by a regression equation which has uncorrelated errors with
equal variance. Of course, the practical problem is that $\Sigma$ may not be known.
We find that,

$var(\hat\beta) = (X^T\Sigma^{-1}X)^{-1}\sigma^2$

Simulation Results
========================================================

initial results by 100 replications


Simulation Results
========================================================

perturbed the data to see robustness of the model more


Application
========================================================
Comparison of results from running the model this way


Application
========================================================
More notes on the comparison

Uses all information together
more statistically powerful this way

maybe also perturb the data to see how robust the estimates are


Conclusions
========================================================
Power model
many applications
flexible model, can handle other growth equeations or biologically relevant functional equations



Future Aims
========================================================

 - Aim 1
incorporate other error structures and mathematical functions that are biologically meaningful

 - Aim 2
 use other type of multivariate responses that are correlated.  Like gene expression and protein expression on same genes.  What SNPs predict these .
 
 Or use expression levels of different cell lines or tissue types.  
 
 You could also maybe incorporate a functional componenet to the expression levels over time. 
 
 __Analysis of Time-Series Gene Expression Data: Methods, Challenges, and Opportunities__
 Monitoring the change in expression patterns over time provides the distinct possibility of unraveling the mechanistic drivers characterizing cellular responses. Gene arrays measuring the level of mRNA expression of thousands of genes simultaneously provide a method of high-throughput data collection necessary for obtaining the scope of data required for understanding the complexities of living organisms. Unraveling the coherent complex structures of transcriptional dynamics is the goal of a large family of computational methods aiming at upgrading the information content of time-course gene expression data. In this review, we summarize the qualitative characteristics of these approaches, discuss the main challenges that this type of complex data present, and, finally, explore the opportunities in the context of developing mechanistic models of cellular response.
 
 
 - Aim 3


Thank you
========================================================

