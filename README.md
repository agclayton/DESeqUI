# DESeqUI
This shiny app is based on the package DESeq2 and provides a user interface for RNAseq data analysis. The app is designed for working with **Trypanosoma brucei** RNAseq data. Feel free to fork this repo and modify it as you wish. Pull requests are welcome.

## Dependencies
- shiny
- DESeq2
- ggplot2
- plotly

## How to use DESeqUI
### Setup
It is recommended to use [RStudio](https://www.rstudio.com) for the following:

- Download the package
- Fire up RStudio and set the package directory as working directory
- open: DESeqUI/app.r in Rstudio
- Click 'Run App' (top right corner)

The program will load its dependencies into the environment (observe RStudios console) and the user interface will appear.

- Click 'Browse' and select the Count-table
- Below the 'Browse'-Button, details about the file can be specified. If you used [TrypRNAseq](https://github.com/klprint/TrypRNAseq) for data preparation, the default can be used
- If the upload of the data worked properly, the first couple of line of the file will be shown in the main-window

Now set up the DESeq2 parameters:

**Column description**: Group the columns in the read file in their corresponding groups. For example, if the first three samples are treated and the following two control, write: treat,treat,treat,ctrl,ctrl

The descriptions can be anything, but be consistent. You could have also written: t,t,t,c,c

Please make sure not to use spacebar, just write down a comma-seperated list.


**Normalization**: Here, the contrast which should be studied can be defined. To keep our previous example, we are interested in significant differences between the treated and the control group. Therefore, we should enter: treat/ctrl

If the data-frame consists of multiple groups, pick the ones you want to study.


**Significance level**: Here you can define which level of significance you want to set. The default DESeq2 value is 0.1 which was kept in DESeqUI.


**Alternative Hypothesis**: As the name implies, this option defines the alternative hypothesis which should be tested for, Depending on the experimental design.

The class enrichment parameters can be left as they are for now.


### Data analysis
- Click on the PCA tab
In the RStudio console you should see the programm working. Depending on the data-size, this can take some minutes. When the calculations are done a PCA-plot (principle component analysis) will be shown. The plot is interactive and by hovering over each point, its identity is shown as a tooltip

- Click on the DESeq tab
A MA-plot will be generated according to the DESeq2 package description. For details please check out the vignette of DESeq2. In short: Each point represents one gene. Points which are colored red were identified as significantly (see above) different in both specified groups.

- Click on the Classes tab
Here, class enrichment analysis can be done. Classes are specified using a currated list of genes (Prof. Christine Clayton). Using the **Class enrichment Parameters** panel, you can setup the analysis. If you are interested in the classes of genes which are upregulated (at least 2-fold) set up as follows: log2FoldChange -> greater, Cutoff log2Fold-Change -> 1. Now, the programm counts the class appearance in the defined subset of genes and tests whether each class appears more often than in the background (the unique gene list), using fishers exact test and adjustment of the p-value using the Benjamini-Hochberg method. The uppermost barplot shows the result of this analysis. Only classes below the user-defined level of significance will be reported. The p-adjusted value is displayed below each bar-group. If a 0 is printed, the p-adjusted value is below 0.001. The exact value is printed to the RStudio console.

The next plot shows a boxplot which reports the distribution of each class within **all** significantly different genes in both groups. The red line indicates the log2FoldCutoff defined previously.

The thrid plot, another barplot, reports the peak time during cellcycle for each subset significant gene. This plot is affected by the log2FoldChange option and the cutoff.

The scatterplot below resembles in its coloring the MA-plot but shows the transcript length for each gene against their log2Fold change. Here, a length-bias could be identified. For example you did a RITseq experiment and the matrix you were using interacts with RNA itself. Therefore the longer the RNA the more it binds to the matrix. A bias like that would be visible in this plot.


- Click on Result Table
This is a preview of the DESeq2 output table. It is subset using the Class Enrichment Parameters defined by the user. The full table can be downloaded using the download button.

## References
Love MI, Huber W and Anders S (2014). “Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.” Genome Biology, 15, pp. 550. doi: [10.1186/s13059-014-0550-8](http://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8)

