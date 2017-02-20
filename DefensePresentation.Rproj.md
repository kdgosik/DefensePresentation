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


Instead of minimizing
========================================================

Instead of minimizing the ordinary least squares estimate we have to find the arg minimum of ...


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


Future Aims
========================================================

 - Aim 1
 - Aim 2
 - Aim 3


Thank you
========================================================

