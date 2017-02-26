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
<img src='images/SystemsMapping.png' width="1500" height="500"> </img>
<footer>
Genotype-Phenotype mapping of various quantitative trait loci
</footer>
</div>



Experimental Design
========================================================
type: incomplete

<small>
Suppose that our data contains n progeny, each of which is genotyped by p markers, such as single nucleotide polymorphisms (SNPs), distributed over different chromosomes
</small>
**update picture**
<div class="midcenter" style="margin-left:-500px; margin-top:0px;">
<img src='images/SNP_Picture.png' style="background-color:transparent; border:0px; box-shadow:none;"> </img>
</div>

***

<br>
$$
\begin{equation}
Phenotype_i = \mathbf{SNP_i}^T\beta_j + \epsilon_i \\
for~j~in~1,\dots,p
\end{equation}
$$
<br>
Where $~\beta_j~$ are the genetic/epistatic effects  
<br>
$$
\epsilon_i \sim N(0, \sigma^2)
$$  


Important Considerations 
========================================================
type: exclaim
incremental: true

 Consider all SNPs at once as possible QTLs (p >> n)  

Screening is predominantly used to quickly reduce the dimensionality. We would need a relativley efficient computational procedure and the resulting model should posses sure screening properties under reasonable assumptions
<small>
 - Including order-2 epistatic effects increases to ${{p}\choose{2}}+p=\frac{p^2+p}{2}$
  - p = 500, gives 125,250 possible candidate SNPs and epistatic effects
 - Including order-3 epistatic effects increases to ${{p}\choose{3}}+{{p}\choose{2}}+p=\frac{p^3+5p}{6}$
  - p = 500 gives 20,833,750 possibie candidate SNPs and epistatic effects
</small>


Variable Selection
========================================================
type: exclaim

<h4> Variable selection is usually implemented in order to handle the high dimensionality of the data.  Many techniques exist including, </h3>

 - LASSO
 - SCAD
 - Elasticnet
 - Dantzig selector 
 
<h4> Interaction selection Methods </h3>
 - Glinternet[[2]](#/48)
 - hierNet [[3]](#/48)


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


FORWARD REGRESSION METHOD
========================================================
type: exclaim

<small>
Algorithm 
 - Step 1: (Initialization) Set $\mathcal{S}^{(0)} = \emptyset$
 - Step 2: (Forward Regression)
   - Evaluation. In the kth step (k ≥ 1), we are given $S^{(k−1)}$. 
    Then, for every $j \in \mathcal{F}/S^{(k−1)}$, we construct a
    candidate model $\mathcal{M}^{(k−1)} = \mathcal{S}^{(k-1)} \cup j$. 
    We then compute $RSS^{(k−1)}$
    - Screen. We then find $a_k = argmin(RSS_{j}^{(k-1)})$ and update 
    $\mathcal{S}^{(k)}=\mathcal{S}^{(k-1)} \cup {a_k}$ accordingly.
 - Step 3: (Solution Path). Iterating Step for n times, which leads to a total of n nested candidate models.  We then collect those models by a solution path $\mathbb{S} = \{\mathcal{S}^{(k)}: 1 \le k \le n\}$
</small>    

FORWARD REGRESSION METHOD
========================================================
type: exclaim

<small>
Assumptions for standard technical conditions needed to show screening consistency of forward regression [[4]](#/48)
 - (C1) Normality assumption.  Assume that both X and $\epsilon$ follow normal distributions.
 - (C2) Covariance matrix: $\lambda_{min}(\mathbf{A})~and~\lambda_{max}(\mathbf{A})$ represent, respectively the smallest and largest eigenvalues of an arbitrary positive definite matrix $\Sigma$.  We assume that there exist two positive constants $0 \lt \tau_{min} \lt \tau_{max} \lt \infty$, such that $2\tau_{min} \lt \lambda_{min}(\Sigma) \lt \lambda_{max}(\Sigma) \lt \frac{1}{2}\tau_{max}$
 - (C3) Regression coefficients.  We assume that $||\beta|| \le \mathcal{C_{\beta}}$ for some constant $\mathcal{C_{\beta}} \gt 0$ and $\beta_{min} \ge \nu_{\beta}n^{\xi_{min}}$ for some $\xi_{min} \gt 0$
 - (C4) Divergence speed of d and $d_0$.  There exists constants 
 $\xi$, $\xi_0$, and $\nu$ such that $log(d) \le \nu n^{\xi_{0}}$, and $\xi + 6\xi_0 + 12\xi_{min} \lt 1$.
</small>


Screening Consistency
========================================================
type: incomplete
<small>
Note that, it is unrealistic to require $\mathcal{T} \in \mathcal{S}$ because this is not guaranteed even in the fixed dimension situation. However, it is indeed possible to have $\mathcal{T} \subset \mathcal{S}^{(k)}$ for some 1 ≤ k ≤ n [[6]](#/48). Otherwise, there exists at least one relevant predictor completely missed by the solution path $\mathcal{S}$.

(Lemma 3 and Lemma 4)

solution path $\mathcal{S}$ to be screening consistent, if  
$$P(\mathcal{T} \subset \mathcal{S}^{(k)} \in \mathcal{S} for some 1 \le k \le n) \rightarrow 1$$

Theorem 1. Under model (2.1) and conditions (C1) - (C4), we have as $n \rightarrow \infty$
$$P(\mathcal{T} \subset \mathcal{S}^{([K \nu n^{2\xi_0 + 4\xi_{min}}])}) \rightarrow 1$$

within
$O(n^{2\xi_0 + 4\xi_{min}})$ steps which is much smaller than the samples size, n under (C4).
</small> 


Genetic Interactions
========================================================
type: exclaim

<small>
 - Referred to as epistatic interactions, they contribute to many complex traits.  
 - Relatively little is known about the specific forms of genetic interactions that are important to heritable phenotypic variation.
 - Mainly reported genetic interactions involving only two loci.
   - Emphasis on gene–gene interactions over highr-order genetic interactions (HGIs) involving three or more loci is rooted in technical issues, rather than biology. [[1]](#/48)
</small>
<div class="midcenter" style="margin-left:150px; margin-top:150px;">
<img src='images/EpistasisScissors.png' height="200" widith="200" style="background-color:transparent; border:0px; box-shadow:none;"> </img>
<footer>
Epistatic Effects Work Like Scissors
</footer>
</div>


Marginality
========================================================
type: incomplete

Forward Selection already includes the marginality principle.  
Relate to genetic data and the necessity of including interactions
 - iForm procedure
  - marginality
  - heredity (both strong and weak)
 - Higher order iForm
 - Functional Mapping iForm
 
 
Heredity Principle
========================================================
type: incomplete

Why assume heredity principle?

 - Statistical efficiency
 - Computational efficient
 - Cost Effective


Model
========================================================
type: incomplete

$$\mathbf{Y} = \mathbf{X^T}\beta^{(1)} + \mathbf{Z^T}\beta^{(2)} + \epsilon$$

Assumptions
 - $X_i$ and $X_jX_k$ are marginally and jointly normal
 - Constants two constants $0 \lt \tau_{min} \lt \frac{1}{4} \lt 1 \lt \tau_{max} \lt \infty$ s.t $2\tau_{min}  \lambda_{min}(\Sigma_{(1)})$
 

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
 


Model Properties
========================================================
type: exclaim

Computational complexity is O(mn) with n being the sample size and m being the number of predictors in the candidate set at that iteration of the procedure, 
 - order-2 $m \le p+n(n+1)~$ for any step.  Making it $O(n^2p+n^4)$
 - order-3 $m \le p+n(n^2-3n+8)~$ for any step. Making it $O(n^2p+n^5-3n^4+8n^3)$

Therefore the complexity grows linearly with p


Model Properties
========================================================
type: incomplete

Remark 2. **Beyond normality.** Lemmas 6, 7, 10 play important roles in the proofs of Theorems
1 and 2. A key assumption is $E(e^{T_0|W_i|^{\alpha}}) \le A_0$ where $W_i$ is (higher) product of predictors. It is easy to see that the condition still holds, using the argument of Lemma 9, if the marginal distributions of
X is sub Gaussian. In particular, Theorem 2 is still true if (C1') holds and the total covariance matrix  $\Sigma$  has bounded eigenvalues asymptotically


Model Theoretical Results
========================================================
type: incomplete

Lemmas 1 - 10 (Mostly 3 and 4)
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



Simulation (truth obeys strong heredity)
========================================================
type: incomplete

<small>

|Model                |T1_tpr |T1_tnr |T1_fpr |T1_fnr |T2_tpr |T2_tnr |T2_fpr |T2_fnr |T3_tpr |T3_tnr |T3_fpr |T3_fnr | Train_MSE| Train_Rsq| Test_MSE| Test_Rsq|RunTime |
|:--------------------|:------|:------|:------|:------|:------|:------|:------|:------|:------|:------|:------|:------|---------:|---------:|--------:|--------:|:-------|
|forward_select       |1.000  |0.999  |0.001  |0.000  |0.000  |1.000  |0.000  |1.000  |0.000  |1.000  |0.000  |1.000  |     3.330|     0.727|    3.490|    0.711|0.757   |
|iform_order_2_weak   |1.000  |1.000  |0.000  |0.000  |0.000  |1.000  |0.000  |1.000  |0.000  |1.000  |0.000  |1.000  |     1.128|     0.907|    1.252|    0.895|5.896   |
|iform_order_2_strong |1.000  |1.000  |0.000  |0.000  |0.000  |1.000  |0.000  |1.000  |0.000  |1.000  |0.000  |1.000  |     1.102|     0.909|    1.198|    0.900|1.557   |
|forward_select_2     |1.000  |1.000  |0.000  |0.000  |0.000  |1.000  |0.000  |1.000  |0.000  |1.000  |0.000  |1.000  |     1.086|     0.910|    1.198|    0.900|25.481  |
|forward_select_3     |1.000  |1.000  |0.000  |0.000  |0.000  |1.000  |0.000  |1.000  |0.000  |1.000  |0.000  |1.000  |     0.992|     0.918|    1.121|    0.906|471.881 |
|iform_order_3_weak   |1.000  |1.000  |0.000  |0.000  |0.000  |1.000  |0.000  |1.000  |0.000  |1.000  |0.000  |1.000  |     1.020|     0.916|    1.135|    0.905|11.346  |
|iform_order_3_strong |1.000  |1.000  |0.000  |0.000  |0.000  |1.000  |0.000  |1.000  |0.000  |1.000  |0.000  |1.000  |     0.968|     0.920|    1.060|    0.911|1.872   |
|glinternet           |1.000  |0.559  |0.441  |0.000  |1.000  |0.982  |0.018  |0.000  |0.000  |1.000  |0.000  |1.000  |     1.246|     0.898|    1.446|    0.880|208.167 |
|hierNet              |1.000  |0.697  |0.303  |0.000  |1.000  |0.976  |0.024  |0.000  |0.000  |1.000  |0.000  |1.000  |     0.906|     0.925|    1.421|    0.882|27.521  |
|Oracle               |NA     |NA     |NA     |NA     |NA     |NA     |NA     |NA     |NA     |NA     |NA     |NA     |     0.953|     0.921|    1.050|    0.912|NA      |
</small>


Simulation (truth obeys weak heredity)
========================================================
type: incomplete

<small>

|Model                |T1_tpr |T1_tnr |T1_fpr |T1_fnr |T2_tpr |T2_tnr |T2_fpr |T2_fnr |T3_tpr |T3_tnr |T3_fpr |T3_fnr | Train_MSE| Train_Rsq| Test_MSE| Test_Rsq|RunTime |
|:--------------------|:------|:------|:------|:------|:------|:------|:------|:------|:------|:------|:------|:------|---------:|---------:|--------:|--------:|:-------|
|forward_select       |1.000  |0.999  |0.001  |0.000  |0.000  |1.000  |0.000  |1.000  |0.000  |1.000  |0.000  |1.000  |  3.325580|  0.731288| 3.480321| 0.715774|4.355   |
|iform_order_2_weak   |1.000  |1.000  |0.000  |0.000  |0.000  |1.000  |0.000  |1.000  |0.000  |1.000  |0.000  |1.000  |  1.119445|  0.909529| 1.200191| 0.900796|8.342   |
|iform_order_2_strong |1.000  |0.992  |0.008  |0.000  |0.000  |1.000  |0.000  |1.000  |0.000  |1.000  |0.000  |1.000  |  1.580267|  0.872314| 1.706937| 0.859151|2.952   |
|forward_select_2     |1.000  |1.000  |0.000  |0.000  |0.000  |1.000  |0.000  |1.000  |0.000  |1.000  |0.000  |1.000  |  1.082909|  0.912358| 1.166813| 0.903782|38.872  |
|forward_select_3     |1.000  |1.000  |0.000  |0.000  |0.000  |1.000  |0.000  |1.000  |0.000  |1.000  |0.000  |1.000  |  0.978966|  0.920820| 1.089431| 0.910194|569.983 |
|iform_order_3_weak   |1.000  |1.000  |0.000  |0.000  |0.000  |1.000  |0.000  |1.000  |0.000  |1.000  |0.000  |1.000  |  1.003161|  0.919024| 1.079398| 0.910740|13.054  |
|iform_order_3_strong |1.000  |0.992  |0.008  |0.000  |0.000  |1.000  |0.000  |1.000  |0.000  |1.000  |0.000  |1.000  |  1.578001|  0.872493| 1.704871| 0.859294|2.787   |
|glinternet           |1.000  |0.469  |0.531  |0.000  |1.000  |0.980  |0.020  |0.000  |0.000  |1.000  |0.000  |1.000  |  0.905914|  0.926629| 1.424976| 0.883336|29.975  |
|hierNet              |1.000  |0.657  |0.343  |0.000  |1.000  |0.973  |0.027  |0.000  |0.000  |1.000  |0.000  |1.000  |  0.856364|  0.930865| 1.412004| 0.884118|33.302  |
|Oracle               |NA     |NA     |NA     |NA     |NA     |NA     |NA     |NA     |NA     |NA     |NA     |NA     |  0.940407|  0.923918| 1.034358| 0.914689|NA      |
</small>


Simulation (truth is anti-heredity)
========================================================
type: incomplete

<small>

|Model                |T1_tpr |T1_tnr |T1_fpr |T1_fnr |T2_tpr |T2_tnr |T2_fpr |T2_fnr |T3_tpr |T3_tnr |T3_fpr |T3_fnr | Train_MSE| Train_Rsq| Test_MSE| Test_Rsq|RunTime |
|:--------------------|:------|:------|:------|:------|:------|:------|:------|:------|:------|:------|:------|:------|---------:|---------:|--------:|--------:|:-------|
|forward_select       |1.000  |1.000  |0.000  |0.000  |0.000  |1.000  |0.000  |1.000  |0.000  |1.000  |0.000  |1.000  |  3.283812|  0.728960| 3.509981| 0.713558|1.005   |
|iform_order_2_weak   |1.000  |0.996  |0.004  |0.000  |0.000  |1.000  |0.000  |1.000  |0.000  |1.000  |0.000  |1.000  |  3.139589|  0.740934| 3.435077| 0.719197|7.866   |
|iform_order_2_strong |1.000  |1.000  |0.000  |0.000  |0.000  |1.000  |0.000  |1.000  |0.000  |1.000  |0.000  |1.000  |  3.283812|  0.728960| 3.509981| 0.713558|2.386   |
|forward_select_2     |1.000  |1.000  |0.000  |0.000  |0.000  |1.000  |0.000  |1.000  |0.000  |1.000  |0.000  |1.000  |  1.081473|  0.910605| 1.171091| 0.903961|29.095  |
|forward_select_3     |1.000  |1.000  |0.000  |0.000  |0.000  |1.000  |0.000  |1.000  |0.000  |1.000  |0.000  |1.000  |  0.988849|  0.918303| 1.095295| 0.910102|548.617 |
|iform_order_3_weak   |1.000  |0.997  |0.003  |0.000  |0.000  |1.000  |0.000  |1.000  |0.000  |1.000  |0.000  |1.000  |  3.155176|  0.739209| 3.447518| 0.718685|13.216  |
|iform_order_3_strong |1.000  |1.000  |0.000  |0.000  |0.000  |1.000  |0.000  |1.000  |0.000  |1.000  |0.000  |1.000  |  3.283812|  0.728960| 3.509981| 0.713558|2.703   |
|glinternet           |1.000  |0.290  |0.710  |0.000  |1.000  |0.971  |0.029  |0.000  |0.000  |1.000  |0.000  |1.000  |  0.844189|  0.930543| 1.577582| 0.871187|26.564  |
|hierNet              |1.000  |0.142  |0.858  |0.000  |1.000  |0.915  |0.085  |0.000  |0.000  |1.000  |0.000  |1.000  |  0.306779|  0.974754| 2.215589| 0.818759|3.417   |
|Oracle               |NA     |NA     |NA     |NA     |NA     |NA     |NA     |NA     |NA     |NA     |NA     |NA     |  0.951558|  0.921322| 1.031376| 0.915384|NA      |
</small>


Simulation (truth is constructed of pure interactions)
========================================================
type: incomplete

<small>

|Model                |T1_tpr |T1_tnr |T1_fpr |T1_fnr |T2_tpr |T2_tnr |T2_fpr |T2_fnr |T3_tpr |T3_tnr |T3_fpr |T3_fnr | Train_MSE| Train_Rsq| Test_MSE|   Test_Rsq|RunTime |
|:--------------------|:------|:------|:------|:------|:------|:------|:------|:------|:------|:------|:------|:------|---------:|---------:|--------:|----------:|:-------|
|forward_select       |NaN    |0.980  |0.020  |NaN    |0.000  |1.000  |0.000  |1.000  |0.000  |1.000  |0.000  |1.000  |  3.316366| 0.0251390| 3.445319| -0.0388060|1.177   |
|iform_order_2_weak   |NaN    |0.972  |0.028  |NaN    |0.000  |1.000  |0.000  |1.000  |0.000  |1.000  |0.000  |1.000  |  3.006933| 0.1151560| 3.180697|  0.0404880|5.840   |
|iform_order_2_strong |NaN    |0.979  |0.021  |NaN    |0.000  |1.000  |0.000  |1.000  |0.000  |1.000  |0.000  |1.000  |  3.294337| 0.0309730| 3.428952| -0.0339540|2.081   |
|forward_select_2     |NaN    |1.000  |0.000  |NaN    |0.000  |1.000  |0.000  |1.000  |0.000  |1.000  |0.000  |1.000  |  1.117353| 0.6694520| 1.169791|  0.6435790|26.396  |
|forward_select_3     |NaN    |1.000  |0.000  |NaN    |0.000  |1.000  |0.000  |1.000  |0.000  |1.000  |0.000  |1.000  |  1.004500| 0.7025090| 1.080699|  0.6705220|530.362 |
|iform_order_3_weak   |NaN    |0.975  |0.025  |NaN    |0.000  |1.000  |0.000  |1.000  |0.000  |1.000  |0.000  |1.000  |  3.043130| 0.1064300| 3.208641|  0.0320860|9.461   |
|iform_order_3_strong |NaN    |0.979  |0.021  |NaN    |0.000  |1.000  |0.000  |1.000  |0.000  |1.000  |0.000  |1.000  |  3.294337| 0.0309730| 3.428952| -0.0339540|2.265   |
|glinternet           |NaN    |0.429  |0.571  |NaN    |1.000  |0.983  |0.017  |0.000  |0.000  |1.000  |0.000  |1.000  |  1.002329| 0.6991235| 1.444997|  0.5614912|145.078 |
|hierNet              |NaN    |0.147  |0.853  |NaN    |1.000  |0.955  |0.045  |0.000  |0.000  |1.000  |0.000  |1.000  |  0.672495| 0.8021830| 1.758235|  0.4666160|4.491   |
|Oracle               |NA     |NA     |NA     |NA     |NA     |NA     |NA     |NA     |NA     |NA     |NA     |NA     |  0.968410| 0.7132750| 1.021514|  0.6886360|NA      |
</small>


Application 1 (C Elegans) Picture of Celegans Genome
========================================================
type: incomplete

 - To genetically dissect the causes of different variability among C. elegans traits, transcript abundances of 20,000 gene transcripts were measured by microarray in developmentally synchronized young adult hermaphrodites of 208 recombinant inbred  advanced intercross lines from a cross between the laboratory strain, N2, and a wild isolate from Hawaii, CB485618

 - Microarray data preprocessed through a normal–exponential convolution background correction and were normalized using quantile standardization  Rockman et al 2010



Application 1 (C Elegans)
========================================================
type: incomplete
incremental: true

<div class="midcenter" style="margin-left:-300px; margin-top:-300px;">
<img src='images/OutputTable_HighDeQTL.png'> </img>
</div>


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
incremental: true


|Coefficient                           | Estimate|      SE| T.value|P.value      |
|:-------------------------------------|--------:|-------:|-------:|:------------|
|(Intercept)                           |  0.18285| 0.07613|   2.402|0.0174 *     |
|AATTC_nn_np_2517_a                    |  0.40013| 0.06509|   6.147|5.13e-09 *** |
|AATTC_nn_np_2815_a                    |  0.15792| 0.06837|   2.310|0.0221 *     |
|CATG_nn_np_3479_a                     |  0.23433| 0.05285|   4.434|1.63e-05 *** |
|CATG_nn_np_1284_a                     |  0.22200| 0.05313|   4.179|4.61e-05 *** |
|AATTC_nn_np_2815_a×AATTC_lm_ll_3034_a |  0.45783| 0.09244|   4.953|1.71e-06 *** |


Application 2 (Mei Trees)
========================================================
type: incomplete
incremental: true


|Coefficient                                              | Estimate|      SE| T.value|P.value      |
|:--------------------------------------------------------|--------:|-------:|-------:|:------------|
|(Intercept)                                              |  0.16859| 0.05801|   2.906|0.00415 **   |
|AATTC_nn_np_2517_a                                       |  0.27773| 0.04396|   6.318|2.27e-09 *** |
|AATTC_nn_np_2815_a                                       |  0.26382| 0.05295|   4.983|1.54e-06 *** |
|CATG_nn_np_3479_a                                        |  0.20767| 0.03467|   5.990|1.23e-08 *** |
|CATG_nn_np_1284_a                                        |  0.04522| 0.04265|   1.060|0.29055      |
|AATTC_nn_np_2815_a×AATTC_lm_ll_3034_a                    |  1.82572| 0.17925|  10.185|< 2e-16 ***  |
|AATTC_nn_np_2815_a×AATTC_hk_hk_278_a                     |  0.25935| 0.03888|   6.671|3.48e-10 *** |
|CATG_lm_ll_3153_a                                        |  0.14877| 0.03491|   4.262|3.36e-05 *** |
|CATG_nn_np_1284_a×AATTC_nn_np_554_a                      |  0.22994| 0.05104|   4.505|1.23e-05 *** |
|AATTC_nn_np_2815_a.AATTC_lm_ll_3034_a×AATTC_nn_np_1615_a | -1.51714| 0.19060|  -7.960|2.39e-13 *** |
|AATTC_nn_np_2815_a×AATTC_nn_np_929_a                     | -0.30805| 0.05477|  -5.624|7.57e-08 *** |
|AATTC_hk_hk_479_d                                        |  0.16044| 0.03443|   4.660|6.37e-06 *** |
|AATTC_nn_np_2517_a×CATG_hk_hk_648_a                      |  0.14537| 0.02840|   5.118|8.33e-07 *** |


Application 2 (Mei Trees)
========================================================
type: incomplete

Show results in graphical form of the different SNPs and how they impact Tree height growth over time learned from running the selection procedure.  
<small>
$\alpha_1 =$ what SNP in the table  
$\alpha_2 =$ what SNP in the table  
$\alpha_3 =$ what SNP in the table  
</small>

<div class="tooltip">
<img src='images/GrowthCurveComparison.png', width = "1000px"></img>
  <span class="tooltiptext">

$\mu_{111} = \mu(t) + \alpha_1(t) + \alpha_2(t) + \alpha_3(t) + i_{12}(t) + i_{13}(t) + i_{23}(t) + i_{123}(t)$
$\mu_{112} = \mu(t) + \alpha_1(t) + \alpha_2(t) - \alpha_3(t) + i_{12}(t) - i_{13}(t) - i_{23}(t) - i_{123}(t)$
$\mu_{121} = \mu(t) + \alpha_1(t) - \alpha_2(t) + \alpha_3(t) - i_{12}(t) + i_{13}(t) - i_{23}(t) - i_{123}(t)$
$\mu_{122} = \mu(t) + \alpha_1(t) - \alpha_2(t) - \alpha_3(t) - i_{12}(t) - i_{13}(t) + i_{23}(t) + i_{123}(t)$
$\mu_{211} = \mu(t) - \alpha_1(t) + \alpha_2(t) + \alpha_3(t) - i_{12}(t) - i_{13}(t) + i_{23}(t) - i_{123}(t)$
$\mu_{212} = \mu(t) - \alpha_1(t) + \alpha_2(t) - \alpha_3(t) - i_{12}(t) + i_{13}(t) - i_{23}(t) + i_{123}(t)$
$\mu_{221} = \mu(t) - \alpha_1(t) - \alpha_2(t) + \alpha_3(t) + i_{12}(t) - i_{13}(t) - i_{23}(t) + i_{123}(t)$
$\mu_{222} = \mu(t) - \alpha_1(t) - \alpha_2(t) - \alpha_3(t) + i_{12}(t) + i_{13}(t) + i_{23}(t) - i_{123}(t)$

</span>
</div>


Application 2 (Mei Trees)
========================================================
type: incomplete

<div class="tooltip">
<img src='images/EpistasisComparison.png', width = "1000px"></img>
  <span class="tooltiptext">

$$i_{12}(t) = [(\mu_{111}(t) + \mu_{112}(t) + \mu_{221}(t) + \mu_{222}(t)) - (\mu_{121}(t) + \mu_{122}(t) + \mu_{211}(t) + \mu_{212}(t))]$$
$$i_{13}(t) =  [(\mu_{111}(t) + \mu_{121}(t) + \mu_{212}(t) + \mu_{222}(t)) - (\mu_{112}(t) + \mu_{122}(t) + \mu_{211}(t) + \mu_{221}(t))]$$
$$i_{23}(t) =  [(\mu_{111}(t) + \mu_{122}(t) + \mu_{211}(t) + \mu_{222}(t)) - (\mu_{112}(t) + \mu_{121}(t) + \mu_{212}(t) + \mu_{221}(t))]$$
$$i_{123}(t) =  [(\mu_{111}(t) + \mu_{122}(t) + \mu_{212}(t) + \mu_{122}(t)) - (\mu_{112}(t) + \mu_{121}(t) + \mu_{211}(t) + \mu_{222}(t))]$$

</span>
</div>


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

<img src='images/GrowthCurveExample.png', width = "1000px"> 


Legendre Polynomials
========================================================
type: incomplete

<img src='images/ExampleDataGrowthCurve.png', width = "1000px"> 


Legendre Polynomials
========================================================
type: incomplete

<img src='images/LegendreFit.png', width = "1000px"> 



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

100% screening consistent but also fitted up to cutoff amount in algorithm  



Application 2 (Mei Trees)
========================================================
type: incomplete
incremental: true


|Coefficient                                                      |  Estimate|       SE|  T.value| P.value|
|:----------------------------------------------------------------|---------:|--------:|--------:|-------:|
|mu_t                                                             |   0.98465| 0.014211| 69.28900|  0.0000|
|CATG_lm_ll_2801_a_P0                                             | -18.34603| 2.330582| -7.87187|  0.0000|
|CATG_lm_ll_2801_a_by_AATTC_lm_ll_2056_a_P0                       |  15.09978| 2.736415|  5.51809|  0.0000|
|CATG_lm_ll_2801_a_by_AATTC_nn_np_2246_a_P0                       |  17.44734| 2.589045|  6.73891|  0.0000|
|CATG_nn_np_3635_a_P0                                             |  35.64628| 2.457669| 14.50410|  0.0000|
|CATG_nn_np_3635_a_by_CATG_lm_ll_2637_a_P0                        | -20.86321| 2.162166| -9.64922|  0.0000|
|AATTC_lm_ll_2056_a_by_CATG_lm_ll_2801_a_by_CATG_nn_np_4200_a_P0  |  11.63512| 2.993626|  3.88663|  0.0001|
|AATTC_nn_np_2246_a_by_AATTC_nn_np_4205_a_by_CATG_lm_ll_2801_a_P0 | -22.10740| 3.305815| -6.68743|  0.0000|
|CATG_nn_np_3635_a_by_AATTC_nn_np_2292_a_P0                       | -16.84440| 2.092784| -8.04880|  0.0000|
|AATTC_hk_hk_753_d_P0                                             | -13.80524| 1.800741| -7.66642|  0.0000|
|AATTC_hk_hk_753_d_by_AATTC_hk_hk_485_a_P0                        |  18.90166| 2.076084|  9.10448|  0.0000|
|CATG_nn_np_3635_a_by_CATG_nn_np_3234_a_P0                        |  -7.26070| 2.454102| -2.95860|  0.0031|
|AATTC_hk_hk_753_d_by_CATG_nn_np_2208_a_P0                        |   7.81921| 2.288519|  3.41671|  0.0006|
|AATTC_hk_hk_485_a_by_AATTC_hk_hk_753_d_by_AATTC_nn_np_3611_a_P0  | -15.91187| 3.036233| -5.24066|  0.0000|
|CATG_lm_ll_2801_a_by_AATTC_nn_np_3395_a_P0                       | -13.76959| 2.301988| -5.98161|  0.0000|
|CATG_lm_ll_1732_a_P0                                             |  -8.98704| 2.578819| -3.48494|  0.0005|
|CATG_lm_ll_1732_a_by_AATTC_nn_np_4608_a_P0                       |  23.97652| 2.767672|  8.66307|  0.0000|
|AATTC_nn_np_4608_a_by_CATG_lm_ll_1732_a_by_CATG_nn_np_3543_a_P0  | -18.61259| 3.349466| -5.55688|  0.0000|
|CATG_lm_ll_789_a_by_CATG_nn_np_3234_a_by_CATG_nn_np_3635_a_P0    | -11.79580| 2.885184| -4.08840|  0.0000|
|CATG_lm_ll_1732_a_by_AATTC_lm_ll_2212_a_P0                       |   8.16593| 2.177668|  3.74985|  0.0002|
|CATG_lm_ll_1732_a_by_CATG_nn_np_2948_a_P0                        |   9.56873| 2.321721|  4.12139|  0.0000|
|AATTC_hk_hk_495_a_by_CATG_lm_ll_1732_a_by_CATG_nn_np_2948_a_P0   |  -8.64983| 2.183079| -3.96222|  0.0001|
|AATTC_nn_np_4608_a_by_CATG_lm_ll_1732_a_by_CATG_lm_ll_2787_a_P1  |  14.27351| 2.135087|  6.68521|  0.0000|



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
type: exclaim

 1. Taylor and Ehrenreich 2014
 2. Lim and Hastie 2014
 3. Bien et al 2013
 4. Wang 2012
 5. Wnag and Leng 2015
 6. Fan and Lv 2008
 7. Rockman et al 2010

Acknowledgments
========================================================
type: incomplete

Thank you
========================================================
type: incomplete
