STATISTICAL MODELS FOR HIGH DIMENSIONAL SCREENING OF GENETIC AND EPIGENETIC EFFECTS
========================================================
author: Kirk Gosik
date: 02/27/2017
css: custom.css
autosize: true
transition: fade
type: exclaim

<div class="midcenter" style="margin-left:0px; margin-top:0px;">
<img src="images/PennStateCollegeOfMedicine.png" style="background-color:transparent; border:0px; box-shadow:none;"> </img>
</div>


Background
========================================================
type: exclaim

- High Throughput Data 
  - Biotechnologies that continually create large scale and complex data
- Biological System
  - Use entire system when predicting a phenotype
  - Network of interaction effects
     - Gene-gene interactions
     - DNA Methylation Analysis
     - Gene-environment interactions
- Statistical Models and analysis are playing an increasing role in  mapping and identifying important quantitative trait loci and other genetic traits


Genotype-Phenotype Map
========================================================
type: exclaim

<div class="midcenter" style="margin-left:-300px; margin-top:-300px;">
<img src='images/SystemsMapping.png' width="120%" height=100%> </img>
</div>

<div class="footer" style="margin-top:-150px;font-size:80%;">
Genotype-Phenotype mapping of various quantitative trait loci</div>



Genetic Interactions
========================================================
type: incomplete

**HGIs**
Genetic interactions (sometimes referred to as epistatic
interactions) contribute to many complex traits (see Glossary).
Despite widespread recognition of this point [1–6],
relatively little is known about the specific forms of genetic
interactions that are important to heritable phenotypic
variation. To date, researchers have mainly reported genetic
interactions involving only two loci (e.g., [7–11]).
However, this emphasis on gene–gene interactions over
HGIs involving three or more loci (Figure 1) is rooted in
technical issues, rather than biology.




Motivation
========================================================
type: incomplete

<div class="midcenter" style="margin-left:-300px; margin-top:-300px;">
<img src='images/EpistasisScissors.png' style="background-color:transparent; border:0px; box-shadow:none;"> </img>
</div>

<div class="footer" style="margin-top:0px;font-size:80%;">
Epistasis Works Like Scissors</div>



Motivation
========================================================
type: incomplete

**There are two important considerations in designing a screening operator. One pinnacle consideration is the low computational requirement. After all, screening is predominantly used to quickly reduce the dimensionality. The other is that the resulting estimator must possess the sure screening property under reasonable assumptions. Otherwise, the very purpose of variable screening is defeated. SIS operates by evaluating the correlations between the response and one predictor at a time, and retaining the features with top correlations.**
\cite{Xiangyu Wang and Chenlei Leng, HOLP paper}




Motivation
========================================================
type: incomplete

<h3> Variable selection is usually implemented in order to handle the high dimensionality of the data.  Many techniques exist including, </h3>

 - LASSO
 - SCAD
 - Elasticnet
 - Dantzig selector 
 
<h3> Interaction selection Methods </h3>
 - Glinternet (TREVOR HASTIE)
 - hierNet (Robert Tibshirani)


Notation
========================================================
type: exclaim

$\mathcal{P_1}$ - main effects predictors  
$\mathcal{P_2}$ - interaction effects predictors  
$\mathcal{M}$ -  Model Set  
$\mathcal{C}$ -  Candidate Set  
$\mathcal{S}$ -  Solution Set  

$\mathcal{T}$ - True Model  
$\mathcal{F}$ - The Full Model

Generic Model and Notation
========================================================
type: incomplete

Assume a linear model

$$
\begin{equation}
Y_i = \mathbf{X_i}^T\beta + \epsilon_i
\end{equation}
$$
Where $(\mathbf{X_i}, Y_i)$ are independent observations
$$
\epsilon \sim N(0, \sigma^2)
$$  
$$
E(Y_i) = 0~~and~~Var(Y_i) = 1
$$



Motivation
========================================================
type: incomplete

 FORWARD REGRESSION METHOD
 
 - Algorithm
  - Step 1: (Initialization) Set $\mathcal{S}^{(0)} = \emptyset$
  - Step 2: (Forward Regression)
    - Evaluation. In the kth step (k ≥ 1), we are given $S^{(k−1)}$. 
    Then, for every $j \in \mathcal{F}/S^{(k−1)}$, we construct a
    candidate model $\mathcal{M}^{(k−1)} = \mathcal{S}^{(k-1)} \cup j$. 
    We then compute $RSS^{(k−1)}$
    - Screen. We then find $a_k = argmin(RSS_{j}^{(k-1)})$ and update 
    $\mathcal{S}^{(k)}=\mathcal{S}^{(k-1)} \cup {a_k}$ accordingly.
  - Step 3: (Solution Path). Iterating Step for n times, which leads to a total of n nested candidate models.  We then collect those models by a solution path $\mathbb{S} = \{\mathcal{S}^{(k)}: 1 \le k \le n\}$
    

Motivation
========================================================
type: incomplete

 FORWARD REGRESSION METHOD
 
 - Assumptions (properties) Standard technical conditions are needed
 
 - (C1) Normality assumption.  Assume that both X and $\epsilon$ follow normal distributions.
 - (C2) Covariance matrix: $\lambda_{min}(\mathbf{A}) and \lambda_{max}(\mathbf{A})$ represent, respectively the smallest and largest eigenvalues of an arbitrary positive definite matrix $\Sigma$.  We assume that there exist two positive constants $0 \lt \tau_{min} \lt \tau_{max} \lt \infty$, such that $2\tau_{min} \lt \lambda_{min}(\Sigma) \lt \lambda_{max}(\Sigma) \lt \frac{1}{2}\tau_{max}$/
 - (C3) Regression coefficients.  We assume that $||\beta|| \le \mathcal{C_{\beta}}$ for some constant $\mathcal{C_{\beta}} \gt 0$ and $\beta_{min} \ge \nu_{\beta}n^{\xi_{min}}$ for some $\xi_{min} \gt 0$
 - (C4) Divergence speed of d and $d_0$.  There exists constants 
 $\xi$, $\xi_0$, and $\nu$ such that $log(d) \le \nu n^{\xi_{0}}$, and $\xi + 6\xi_0 + 12\xi_{min} \ lt 1$.



Motivation
========================================================
type: incomplete

Screening Consistency

Note that, it is unrealistic to require $\mathcal{T} \in \mathcal{S}$ because this is
not guaranteed even in the fixed dimension situation. However,
it is indeed possible to have $\mathcal{T} \subset \mathcal{S}^{(k)}$ for some 1 ≤ k ≤ n (Fan
and Lv 2008). Otherwise, there exists at least one relevant predictor
completely missed by the solution path $\mathcal{S}$.


Motivation
========================================================
type: incomplete

solution path $\mathcal{S}$ to be screening consistent, if  
$$P(\mathcal{T} \subset \mathcal{S}^{(k)} \in \mathcal{S} for some 1 \le k \le n) \rightarrow 1$$

Theorem 1. Under model (2.1) and conditions (C1) - (C4), we have as $n \rightarrow \infty$
$$P(\mathcal{T} \subset \mathcal{S}^{([K \nu n^{2\xi_0 + 4\xi_{min}}])}) \rightarrow 1$$

within
$O(n^{2\xi_0 + 4\xi_{min}})$ steps which is much smaller than the samples size, n under (C4).
 






Background
========================================================
type: incomplete

relate to genetic data and the necessity of including interactions
 - iForm procedure
  - marginality
  - heredity (both strong and weak)
 - Higher order iForm
 - Functional Mapping iForm
 
 
Background
========================================================
type: incomplete

Why assume heredity principle?

 - Statistical efficiency
 - Computational efficient
 - Cost Effective


Background
========================================================
type: incomplete



Model
========================================================
type: incomplete

$$\mathbf{Y} = \mathbf{X^T}\beta^{(1)} + \mathbf{Z^T}\beta^{(2)}$$

Assumptions
 - X_i and X_jX_k are marginally and jointly normal
 - Constants two constants $0 \lt \tau_{min} \lt \frac{1}{4} \lt 1 \lt \tau_{max} \lt \infty$ s.t $2\tau_{min}  \lambda_{min}(\Sigma_{(1)})$
 


Model Properties
========================================================
type: incomplete

Computational Complexity is linear in p.

Computational complexity is O(mn) with n being the sample size and m being the number of predictors in the candidate set at that iteration of the procedure

order-2 m <= 


Model Properties
========================================================
type: incomplete

Remark 2. **Beyond normality.** Lemmas 6, 7, 10 play important roles in the proofs of Theorems
1 and 2. A key assumption is $E(e^{T_0|W_i|^{\alpha}}) \le A_0$ where $W_i$ is (higher) product of predictors. It is easy to see that the condition still holds, using the argument of Lemma 9, if the marginal distributions of
X is sub Gaussian. In particular, Theorem 2 is still true if (C1') holds and the total covariance matrix  $\Sigma$  has bounded eigenvalues asymptotically


Model Theoretical Results
========================================================
type: incomplete

Lemmas 1 - 10 
Theorem 2

show it holds for all three scenarios
 - iForm with order-2 interactions
  - (C2a) and (C4a)
 - iForm with higher order interactions
  - (C2b) and (C4b)
 - iForm with generalized least squares approach
  - show transformation from general case with correlated errors to making have uncorrelated errors

Model Theoretical Results
========================================================
type: incomplete

$$ BIC_1 = n*log(RSS/n) + k*log(n) $$  
$$ BIC_2 = n*log(RSS/n) + k*(log(n) + 2*log(d^{\star})) $$

derived BIC2 by controlling the false discovery rate (FDR) and showed that it is selection consistent if $d_0 = O(n\xi)$ for some $\xi > 0$

Wang (2009) showed its selection consistency for FS under ultra-high dimensional setup $d_0 = O(e^(n\xi))$.



Model (iForm)
========================================================
type: incomplete

 iForm, still needs updated [Go to FS Algorithm](#/7)
 
 - Algorithm
  - Step 1: (Initialization) 
  Set $\mathcal{S}^{(0)} = \emptyset$, $\mathcal{M}_0 = \emptyset$ and $\mathcal{C}_0 = \mathcal{P_1}$
  
  - Step 2: (Selection)
    - Evaluation. In the kth step (k ≥ 1), we are given $S^{(k−1)}$. 
    Then, for every $j \in \mathcal{F}/S^{(k−1)}$, we construct a
    candidate model $\mathcal{M}^{(k−1)} = \mathcal{S}^{(k-1)} \cup j$. 
    We then compute $RSS^{(k−1)}$
    - Screen. We then find $a_k = argmin(RSS_{j}^{(k-1)})$ and update 
    $\mathcal{S}^{(k)}=\mathcal{S}^{(k-1)} \cup {a_k}$ accordingly.
    
    (Selection) In the kth step with given $\mathcal{S}^{(k-1)}$, $\mathcal{C}^{k−1}$ and $\mathcal{M}^{k−1}$, forward regression is used to select one more predictor from $\mathcal{C}^{k−1}/ \mathcal{S}^{k−1}$ into the model. We add the selected one into $\mathcal{S}^{k−1}$ to get $\mathcal{S}^k$. We also update $\mathcal{C}^k$ and $\mathcal{M}^k$ if the newly selected predictor is a main effect. Otherwise, $\mathcal{C}^k = \mathcal{C}^{k−1}$ and $\mathcal{M}^k = \mathcal{M}^{k−1}$.
    
  - Step 3: (Solution Path). Iterating Step 2, for d times, which leads to a total of n nested candidate models.  We then collect those models by a solution path $\mathbb{S} = \{\mathcal{S}^{(k)}: 1 \le k \le n\}$
 
 
 
 
 
 Model (iForm HGI)
========================================================
type: incomplete

 iForm, still needs updated
 
 - Algorithm
  - Step 1: (Initialization) 
  Set $\mathcal{S}^{(0)} = \emptyset$, $\mathcal{M}_0 = \emptyset$ and $\mathcal{C}_0 = \mathcal{P_1}$
  
  - Step 2: (Selection)
    - Evaluation. In the kth step (k ≥ 1), we are given $S^{(k−1)}$. 
    Then, for every $j \in \mathcal{F}/S^{(k−1)}$, we construct a
    candidate model $\mathcal{M}^{(k−1)} = \mathcal{S}^{(k-1)} \cup j$. 
    We then compute $RSS^{(k−1)}$
    - Screen. We then find $a_k = argmin(RSS_{j}^{(k-1)})$ and update 
    $\mathcal{S}^{(k)}=\mathcal{S}^{(k-1)} \cup {a_k}$ accordingly.
    
    (Selection) In the kth step with given $\mathcal{S}^{(k-1)}$, $\mathcal{C}^{k−1}$ and $\mathcal{M}^{k−1}$, forward regression is used to select one more predictor from $\mathcal{C}^{k−1}/ \mathcal{S}^{k−1}$ into the model. We add the selected one into $\mathcal{S}^{k−1}$ to get $\mathcal{S}^k$. We also update $\mathcal{C}^k$ and $\mathcal{M}^k$ if the newly selected predictor is a main effect. Otherwise, $\mathcal{C}^k = \mathcal{C}^{k−1}$ and $\mathcal{M}^k = \mathcal{M}^{k−1}$.
    
  - Step 3: (Solution Path). Iterating Step 2, for d times, which leads to a total of n nested candidate models.  We then collect those models by a solution path $\mathbb{S} = \{\mathcal{S}^{(k)}: 1 \le k \le n\}$
 



Simulation Results
========================================================
type: incomplete

Show Results from HigherOrder Simulation results comparing all models
  - show for different heredity structures
    - leave out anti-heredity
    - leave out only interaction model


Simulation Results
========================================================
type: incomplete

second simulation results of different heredity structure


Application 1 (C Elegans)
========================================================
type: incomplete

<div class="midcenter" style="margin-left:-300px; margin-top:-300px;">
<img src='images/OutputTable_HighDeQTL.png'> </img>
</div>

<div class="footer" style="margin-top:-150px;font-size:80%;">
This text is my footer text<br>
Normally, the top of the text block would appear 90% of the way down the slide (i.e., at the bottom), because that is how the 'top' property is defined in the 'footer' style. However this particular footer wraps onto multiple lines, so we shift it up by specifying a negative value for the 'margin-top' property.</div>


Application 1 (C Elegans)
========================================================
type: incomplete

<div class="midecenter" style="margin-left:-300px; margin-top:-300px;">
<img src='images/iForm_Network_Chr1.jpeg'> </img>
</div>

<div class="footer" style="margin-top:-150px;font-size:80%;">
</div>


Application 2 (Mei Trees)
========================================================
type: incomplete

Mei tree comparison between running the order-2 compared to the order-3 version of the iForm selection

Application 2 (Mei Trees)
========================================================
type: incomplete

Show results in graphical form of the different SNPs and how they impact Tree height growth over time learned from running the selection procedure.  

Show Epistasis equations that were solved to figure out the epistasis.  

$\mu_{111} = \mu(t) + \alpha_1(t) + \alpha_2(t) + \alpha_3(t) + i_{12}(t) + i_{13}(t) + i_{23}(t) + i_{123}(t)$
$\mu_{112} = \mu(t) + \alpha_1(t) + \alpha_2(t) - \alpha_3(t) + i_{12}(t) - i_{13}(t) - i_{23}(t) - i_{123}(t)$
$\mu_{121} = \mu(t) + \alpha_1(t) - \alpha_2(t) + \alpha_3(t) - i_{12}(t) + i_{13}(t) - i_{23}(t) - i_{123}(t)$
$\mu_{122} = \mu(t) + \alpha_1(t) - \alpha_2(t) - \alpha_3(t) - i_{12}(t) - i_{13}(t) + i_{23}(t) + i_{123}(t)$
$\mu_{211} = \mu(t) - \alpha_1(t) + \alpha_2(t) + \alpha_3(t) - i_{12}(t) - i_{13}(t) + i_{23}(t) - i_{123}(t)$
$\mu_{212} = \mu(t) - \alpha_1(t) + \alpha_2(t) - \alpha_3(t) - i_{12}(t) + i_{13}(t) - i_{23}(t) + i_{123}(t)$
$\mu_{221} = \mu(t) - \alpha_1(t) - \alpha_2(t) + \alpha_3(t) + i_{12}(t) - i_{13}(t) - i_{23}(t) + i_{123}(t)$
$\mu_{222} = \mu(t) - \alpha_1(t) - \alpha_2(t) - \alpha_3(t) + i_{12}(t) + i_{13}(t) + i_{23}(t) - i_{123}(t)$

Show the epistasis graph. 

$$i_{12}(t) = [(\mu_{111}(t) + \mu_{112}(t) + \mu_{221}(t) + \mu_{222}(t)) - (\mu_{121}(t) + \mu_{122}(t) + \mu_{211}(t) + \mu_{212}(t))]$$
$$i_{13}(t) =  [(\mu_{111}(t) + \mu_{121}(t) + \mu_{212}(t) + \mu_{222}(t)) - (\mu_{112}(t) + \mu_{122}(t) + \mu_{211}(t) + \mu_{221}(t))]$$
$$i_{23}(t) =  [(\mu_{111}(t) + \mu_{122}(t) + \mu_{211}(t) + \mu_{222}(t)) - (\mu_{112}(t) + \mu_{121}(t) + \mu_{212}(t) + \mu_{221}(t))]$$
$$i_{123}(t) =  [(\mu_{111}(t) + \mu_{122}(t) + \mu_{212}(t) + \mu_{122}(t)) - (\mu_{112}(t) + \mu_{121}(t) + \mu_{211}(t) + \mu_{222}(t))]$$



Application 2 (Mei Trees)
========================================================
type: incomplete

switch out to showing R code, 
 - maybe add the actual data points to this graph
 - how many trees make up each group
 
<img src='images/GrowthCurveComparison.png', width = "1000px">   
<img src='images/EpistasisComparison.png', width = "1000px">



Mei Tree as Motivation
========================================================
type: incomplete

Previous results were from pre-fitted growth curves fit to each tree height, following the asymptotic logistic growth curve.  

Use Mei Tree results and data to set up motivation for functional response variable and using the iform procedure.  How can we incorporate all information from the response variable into one model?  



Background of Functional Response
========================================================
type: incomplete

 - Regression as Linear Combination of Basis Functions
 - Legendre Polynomials to model genetic effect of each SNP over time
 - Generalized Least Squares for correlated errors (AR1 model assumed)
 
In these sections talk about Treating linear regression as basis functions, legendre polynomials and generalized least squares. 

Linear regression when there is a certain degree of correlation between the residuals in a regression model. In these cases, ordinary least squares and weighted least squares can be statistically inefficient, or even give misleading inferences



Generalized Least Squares
========================================================
type: incomplete

[GLS](https://en.wikipedia.org/wiki/Generalized_least_squares)

The GLS estimator is unbiased, consistent, efficient, and asymptotically normal:

AR(1) model assumed for our purposes.  




Legendre Polynomials
========================================================
type: incomplete

 - Orthogonal and therefore variables will not be correlated
 - variety of fits up to the order of the researcher's choosing
 - easily implemented as basis functions 
 -- saves on computational cost and time over using splines or other basis functions learned from the data
 

$$
\begin{equation}
\begin{split}
P_n(x) & = \frac{1}{2^n}\sum_{k=0}^{n}{{n}\choose{k}}^2(x-1)^{n-k}(x+1)^k \\
& = \sum_{k=0}^{n}{{n}\choose{k}}{{-n-1}\choose{k}}{\left(\frac{1-x}{2}\right)}^k \\
& = 2^{-n}\sum_{k=0}^{n} x^k {{n}\choose{k}}{{\frac{n+k+1}{2}}\choose{k}} \\
\end{split}
\end{equation}
$$



<img src="DefensePresentation.Rproj-figure/Legendre-fig-1.png" title="First 10 Legendre Polynomials" alt="First 10 Legendre Polynomials" width="80%" style="display: block; margin: auto;" />

Legendre Polynomials
========================================================
type: incomplete

show toy example of how the polynomials could model the genetic effect well


Legendre Polynomials
========================================================
type: incomplete

<img src='images/GrowthCurveExample.png', width = "1000px"> 


Legendre Polynomials
========================================================
type: incomplete

<img src='images/LegendreFit.png', width = "1000px"> 


Legendre Polynomials
========================================================
type: incomplete

<img src='images/ExampleDataGrowthCurve.png', width = "1000px"> 


Generalized Least Squares
========================================================
type: incomplete

We have assumed that $$var(\epsilon) = \sigma^2I$$ when the response is static but if we have correlated errors like when we have repeated measurements over time $$var(\epsilon) = \sigma^2\Sigma$$ where sigma^2 is unknown but Sigma is known.  We can use Generalized least squares, Instead of minimizing the ordinary least squares estimate we have to find the arg minimum of 
$$ (y - X\beta)^T\Sigma^{-1}(y - X\beta) $$
which is solved by 
$$ \hat\beta = (X^T\Sigma^{-1}X)^{-1}X^T\Sigma^{-1}y $$
since we can write $$\Sigma = SS^T$$, where S is a triangular matrix using the Choleski Decomposition, we have
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


iForm with Legendre Model
========================================================
type: incomplete

$$ \alpha_j(t) = (L_0(t), L_1(t), ... , L_s(t))*(u_{j0}, u_{j1},...,u_{js})^T $$
$$ \beta_j(t) = (L_0(t), L_1(t), ... , L_{s'}(t))*(v_{j0}, v_{j1},...,v_{js'})^T $$

$y(t) = \mu(t) + \sum_{j=1}^{J}\alpha_j(t)\xi_j + \sum_{k=1}^{K}\beta_k(t)\zeta_k + \sum_{I_1<I_2=1}^{I}\gamma_I^{aa}(t) \xi_{I_1}\xi_{I_2} + \sum_{I_1<I_2=1}^{I}\gamma_I^{ad}(t) \xi_{I_1}\zeta_{I_2} + \sum_{I_1<I_2=1}^{I} \gamma_I^{da}(t) \zeta_{I_1}\xi_{I_2} + \sum_{I_1<I_2=1}^{I} \gamma_I^{dd}(t)\zeta_{I_1}\zeta_{I_2} + \epsilon(t)$



 
 Model (iForm with Legendre Model)
========================================================
type: incomplete

 iForm with Legendre Model, still needs updated
 
 - Algorithm
  - Step 1: (Initialization) 
  Set $\mathcal{S}^{(0)} = \emptyset$, $\mathcal{M}_0 = \emptyset$ and $\mathcal{C}_0 = \mathcal{P_1}$
  
  - Step 2: (Selection)
    - Evaluation. In the kth step (k ≥ 1), we are given $S^{(k−1)}$. 
    Then, for every $j \in \mathcal{F}/S^{(k−1)}$, we construct a
    candidate model $\mathcal{M}^{(k−1)} = \mathcal{S}^{(k-1)} \cup j$. 
    We then compute $RSS^{(k−1)}$
    - Screen. We then find $a_k = argmin(RSS_{j}^{(k-1)})$ and update 
    $\mathcal{S}^{(k)}=\mathcal{S}^{(k-1)} \cup {a_k}$ accordingly.
    
    (Selection) In the kth step with given $\mathcal{S}^{(k-1)}$, $\mathcal{C}^{k−1}$ and $\mathcal{M}^{k−1}$, forward regression is used to select one more predictor from $\mathcal{C}^{k−1}/ \mathcal{S}^{k−1}$ into the model. We add the selected one into $\mathcal{S}^{k−1}$ to get $\mathcal{S}^k$. We also update $\mathcal{C}^k$ and $\mathcal{M}^k$ if the newly selected predictor is a main effect. Otherwise, $\mathcal{C}^k = \mathcal{C}^{k−1}$ and $\mathcal{M}^k = \mathcal{M}^{k−1}$.
    
  - Step 3: (Solution Path). Iterating Step 2, for d times, which leads to a total of n nested candidate models.  We then collect those models by a solution path $\mathbb{S} = \{\mathcal{S}^{(k)}: 1 \le k \le n\}$
 




Simulation Results
========================================================
type: incomplete

initial results by 100 replications


Simulation Results
========================================================
type: incomplete

perturbed the data to see robustness of the model more


Application
========================================================
type: incomplete

Comparison of results from running the model this way

Mei Trees results with running functional model


Application
========================================================
type: incomplete

More notes on the comparison

Uses all information together
more statistically powerful this way

maybe also perturb the data to see how robust the estimates are


Conclusions
========================================================
type: incomplete

Power model
many applications
flexible model, can handle other growth equations or biologically relevant functional equations



Future Aims
========================================================
type: incomplete

 - Aim 1
incorporate other error structures and mathematical functions that are biologically meaningful

 - Aim 2
 use other type of multivariate responses that are correlated.  Like gene expression and protein expression on same genes.  What SNPs predict these .
 
 Or use expression levels of different cell lines or tissue types.  
 
 You could also maybe incorporate a functional component to the expression levels over time. 
 
 __Analysis of Time-Series Gene Expression Data: Methods, Challenges, and Opportunities__
 Monitoring the change in expression patterns over time provides the distinct possibility of unraveling the mechanistic drivers characterizing cellular responses. Gene arrays measuring the level of mRNA expression of thousands of genes simultaneously provide a method of high-throughput data collection necessary for obtaining the scope of data required for understanding the complexities of living organisms. Unraveling the coherent complex structures of transcriptional dynamics is the goal of a large family of computational methods aiming at upgrading the information content of time-course gene expression data. In this review, we summarize the qualitative characteristics of these approaches, discuss the main challenges that this type of complex data present, and, finally, explore the opportunities in the context of developing mechanistic models of cellular response.
 
 
 - Aim 3


Referneces
========================================================
type: incomplete

Acknowledgments
========================================================
type: incomplete

Thank you
========================================================
type: incomplete
