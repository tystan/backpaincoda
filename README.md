# backpaincoda

Repository containing code and files to generate the supplementary material for the article:
"Are reallocations of time between physical activity, sedentary behaviour, and sleep associated with low back pain? A compositional data analysis." [under journal review]


## Repository file structure

* In the top level directory there is:
    + The R analysis script  [`back-pain.R`](https://github.com/tystan/backpaincoda/tree/main/back-pain.R) that performs all the pre-processing of the data, the analysis and output of results/figures
    + [`back-pain.pdf`](https://github.com/tystan/backpaincoda/blob/main/back-pain.pdf) is the supplementary material as a transcript of the R session (code and outputs) to perform the pre-processing of the data, the analysis and output of results/figures
    + The `.qmd` file [`back-pain.qmd`](https://github.com/tystan/backpaincoda/tree/main/back-pain.qmd) is the [Quarto](https://quarto.org/) (Rmarkdown/Jupiter Notebook-esque documents) that generates the corresponding `*.pdf` document [`back-pain.pdf`](https://github.com/tystan/backpaincoda/blob/main/back-pain.pdf)
* [dat/](https://github.com/tystan/backpaincoda/tree/main/dat) contains the analysis data
* [res/](https://github.com/tystan/backpaincoda/tree/main/res) contains the results used for the time-use reallocation plots/predictions with 95% bootstrapped confidence intervals
* [fig/](https://github.com/tystan/backpaincoda/tree/main/fig) contains stand alone, 600 dpi `.png` images of the time-use reallocation plots/predictions (see example embedded below)


## Figures

The generated figures can be found in the [fig/](https://github.com/tystan/backpaincoda/tree/main/fig) directory.

![](https://github.com/tystan/backpaincoda/blob/main/fig/lbp_intens_negbin_abs_v2.png)
Figure: *Associations of reallocating time between 24-hour movement behaviours with the intensity of low back pain (n = 1660). The analyses were adjusted for age, sex, body mass index, smoking, stress, education, and socio-economic status. SB, sedentary behaviour; LPA, light physical activity; MVPA, moderate-to-vigorous physical activity. Note that absolute difference estimates (% intensity scale) are based upon an 'average' participant in the sample and may will differ for different participant covariates.*






