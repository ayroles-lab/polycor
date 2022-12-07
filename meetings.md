# 14 Feb, 2019

TODOs
- [ ] (Andrey) Update `relmer.R` by adding option to specify relationship matrices via `zmat`: `GRM = ZZ'`, 
  where `GRM` is a NxN matrix, while `Z` is a NxM matrix (N = sample size ~ 10K, M = the number of causal SNPs ~ 100).
  This update will assist in computation speed up for dennse GRM matrices.
- [ ] Write new functions for data simulation in `simFunctions.R`: move from family-based sample to unrealted + SNP-based GRM.
  - The original simulation functions are implemented like this. I'll (Simon) have a look at them and update if needed.
- [ ] Create diagrams with causal relationships for t1, t2 and SNPs. 
- An example R package for drawing DAGs: [ggdag](https://github.com/malcolmbarrett/ggdag). See also the issue about [plotting moderators](https://github.com/malcolmbarrett/ggdag/issues/6).
- Mediator vs moderator via [link](https://www.theanalysisfactor.com/five-common-relationships-among-three-variables-in-a-statistical-model/)
- [Article](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.169.4836&rep=rep1&type=pdf) on decomposition of the moderator model 
- [ ] (Andrey) Review the LDSC paper, where `rhog` is defined via SNP-based GWAS summary statistics.

List of hypotheses for testing on real datasets
- Can random-slope model capture the `G_slope` variance and distinguish it from `rhog`?
- Are there any pair of traits that (i) has no genetic correlation `rhog`; (ii) exibits the `G_slope` variance?
  - If we could find any example of this (or even just `rhog` < `G_slope`), it would be a very nice case to support the "conceptual" part of this paper. As I picture it, this is where we try to articulate how pleiotropy (`rhog`) and genetic control of individual correlation (`G_slope`) are two different phenomenon. 
- Perform association tests on SNPs (similar to longitudinal GWAS).
- Once the polygenic effect of `G_slop` is established, individuals are stratified based on magnitude of correlation. Is such stratification linked to a disease?
- Possibly try to perform causal inference by fitting two random slope models, one "in each direction". Relating this to the diagram below, the two models would estimate `G_slope1` and `G_slope2` respectively. We would need to explore this with simulations first.

Possible non-genetic applications
- For a given grouping variable, e.g. environment, try to find group-level correlations using the same-random slope model.

### Causality Diagram
![](https://github.com/hemostat/polycor/blob/master/illustrations/causalityDiagram.JPG)
This is an attempt to visualize what we are simulating, and what our models capture. `t1` and `t2` are our two traits of interest. `G_pleio` is the pleiotropic effect, quantified by `rhog`. `G_add1` and `G_add2` are the additive genetic effects acting independently on the two traits. `G_slope1` and `G_slope2` are the genetic effects on the individual correlation/causality. In our latest simulation, we simulate all of these effects simultaneously (except `G_slope2`) and show that the slope_add model combined with bivar can distinguish them. 

If we want to pursue the causal inference angle, that might involve comparing `G_slope1` and `G_slope2`, estimated from random slope models fitted in different "directions". 

Based on our discussion, I also added a potential third trait `t3`. Such a trait could be the cause of an observed correlation between `t1` and `t2`. In that case, a genetic effect `G_add3` on `t3` should look like pleiotropy when analyzing `t1` and `t2` right? In other words, `G_add3 ~ G_pleio`. So such a signal should not be confounded with `G_slope`. 

# 1 Feb, 2019
Slide 21 [here](https://docs.google.com/presentation/d/1i_GdvuWhyWT5TzZH279RPxKEL5azU54mxBdcuBW4RyY/edit#slide=id.g4e83fb7b4c_4_1) gives a very brief summary of where we are at the moment. Given these results, we should in theory be able to distinguish pleiotropy (genetic cor) from ind level cor (h2.cov) by comparing a slope+add model and a "simple" bi-variate model fitted to the same data. At least as long as we dont't have both phenomenon occuring at the same time. 
- if h2.cov > 0 (slope_add model) & rhog = 0 (bi-var model): 
   - individual-level corr. and NOT pleiotropy
- if h2.cov = 0 (slope_add model) & rhog > 0 (bi-var model): 
   - pleiotropy and NOT individual-level corr.  
- if h2.cov > 0 (slope_add model) & rhog > 0 (bi-var model): 
   - We haven't tried to simulate the scenario where have both pleiotropy and individual-level corr. So we don't yet know if this case implies that scenario (but it should)
   
It would be preferable if we had one extended bi-variate model that could distinguish the two. This would be more elegant, and it could be tested against a simpler model using a LRT. 

## To Do (in order or priority)
- Simulate the pleiotropy + ind level cor scenario (last row in slide 21 above) - Andrey
- Run simulations with parameters based on real data - Simon
    - Julien is going to provide the data
- Continue working on an extended bi-variate model, capable of simultaneously estimating h2.cov and rhog (last column in slide 21 above)
    - We might end up having to drop this idea because it's not practically feasible. But it's worth spending some more time trying.
  
  
# 24 Jan, 2019

- When data simulated under bi-variate model, is the random-slope data fitting immune to genetic corr.?
    - (Simon says): What I'm saying is that when data is simulated under the bi-variate model, there should be no within individual correlation (could be validated by a simple scatter plot per individual). Whether the random slope is immune to "variance leakage" I'm not sure. My simulations suggested that there was leakage from additive effects but not from genetic corr, but Andreys results suggest otherwise. I'll revisit this.
    - (Andrey says): Indeed, my simulations suggest similar conclusions; see two figures in [this post](https://github.com/hemostat/polycor/issues/5#issuecomment-455611571) on issue 5.
    - TODO (Andrey): increase the sample size (n) and check whether the "leakage" tends to zero for large n.
- When data simulated under random-slope model, can the bi-variate partially capture the variance of ind.-level corr.?
    - (Simon says): Not sure how much money I would bet here :) But it is possible I think.
    
A nice diagram:

 - horizontal arrows are causal correlations at ind. level
 - vertical arrows are causal pleiotropy
```
  /------\
y1        y2
  \------/

   \   /
     G
```

- [x] TODO (Andrey): do a simple experiment with repeated measurements to understand ind.-level cor.
   - data simulated under bi-variate model + rep.: plot ind.-level corr and observse a noise
   - data simulated under random-slope model + rep: plot ind.-level corr and observse variation in ind.-level corr.
   - See [02-learn-cor.m](https://github.com/hemostat/polycor/blob/master/projects/02-learn-cor/02-learn-cor.md)
- (Simon says): The goal here is to find out if our models can distinguish the different scenarios. The figures below illustrate the three scenarios we dicussed (like the diagram above)

### Figure 1
![](https://github.com/hemostat/polycor/blob/master/illustrations/pleio.png)

### Figure 2
![](https://github.com/hemostat/polycor/blob/master/illustrations/indCor.png)

### Figure 3
![](https://github.com/hemostat/polycor/blob/master/illustrations/indCor_plusAdd.png)

Possible applications beyond genetics
- Longitudinal/repeated medical records (ARIC, Milieu Interieur)

TODO:
- Experiment with data simulation under random-slope model + no genetics
   - compare the same pair of model fitting: bi-variate vs. random slope

# 19 Dec, 2018

A good example of the phenomenon: two traits, blood preasure and heart rate.
- the genetic corr. (`rhog`) estimated on population level is a single number
- `rhog` heterogeneity across conditions/environments: a person after running vs a sitting person
- `rhog` heterogeneity across individuals: younger/helthier/genetically different

Data sim. scenario 1 (random-slope)
- `y1 ~ rnorm(n)`, the first trait has nothing about genteics
- `g <- mvnorm(n, 0, G)`
- `y2 ~ g y1`, the second trait has a genetic component, which is individual-level genetic correlations.

Data sim. scenario 2 (genetically correlated traits)

|              | Model fitting 1 | Model fiting 2|
|--------------|-----------------|---------------|
| Data sim. 1  | done | ? |
| Data sim. 2  | done | ? |

Real data sources

- Medical records with bloop phenotypes
- Pylogenetic data 

TODO:

- Andrey: implement Model fiting 2 & perform experiments with Data sim. scenario 1/2
