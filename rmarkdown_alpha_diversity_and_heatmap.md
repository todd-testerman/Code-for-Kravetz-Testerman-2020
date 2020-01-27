Paper\_R\_Markdown
================
Todd Testerman
1/24/2020

  - [Data Import and Alpha Diversity
    Calculation](#data-import-and-alpha-diversity-calculation)
  - [Heatmap](#heatmap)

# Data Import and Alpha Diversity Calculation

*Load
packages*

``` r
pacman::p_load(phyloseq, BiocManager, qiime2R, tidyverse, viridis, pheatmap)
```

*Import data from Q2 artifacts into Phyloseq*

``` r
#Import qza files into phyloseq object
phy = qza_to_phyloseq("table-all-data-final.qza", "rooted_tree.qza", "taxonomy.qza", "All_data_with_NAFLD_PNPLA3txt.txt", tmp = "C:/tmp")
#Remove Greengenes prefixes
phy@tax_table = gsub("k__|p__|o__|c__|f__|g__|s__", "", phy@tax_table)
```

*Rarefaction of data and calculation of Shannon Index value using
estimate\_richness*

``` r
#Removing two samples without necessary metadata
phy_no_19 = subset_samples(phy, Stool_DNA_number != 19)
phy_no_19_87 = subset_samples(phy_no_19, Stool_DNA_number != 87)
phy = phy_no_19_87
#Adjust NAFLD status for one subject due to data entry error
phy@sam_data$NAFLD[73] = "yes"
#Rarefy to 10,000 reads
phy_rarefy = rarefy_even_depth(phy, 10000)
```

    ## You set `rngseed` to FALSE. Make sure you've set & recorded
    ##  the random seed of your session for reproducibility.
    ## See `?set.seed`

    ## ...

    ## 274OTUs were removed because they are no longer 
    ## present in any sample after random subsampling

    ## ...

``` r
#Estimate Shannon alpha diversity
alpha_diversity <- estimate_richness(phy_rarefy, measures = "Shannon")
#these values were then exported for plotting in Prism
```

# Heatmap

*Heatmap Generation Setup*

``` r
#Adjusting sample names
sample_names(phy) = phy@sam_data$Sample_ID
#Agglomerate to the genus level
phy_genus_glom = tax_glom(phy, "Genus", NArm = T)
#Filter low abundance taxa (if less than 250 reads summed across all samples, remove)
phy_tax_glom_gen_trimmed = filter_taxa(phy_genus_glom, function(x) sum(x) > 250, TRUE)
#Normalize all samples to 10,000 reads. *THIS IS DIFFERENT THAN RAREFACTION*
phy_tax_glom_gen_trimmed_normalize = transform_sample_counts(phy_tax_glom_gen_trimmed, function(x) (x / sum(x))*10000)
#Convert genus names to a factor
genfac = factor(tax_table(phy_tax_glom_gen_trimmed_normalize)[, "Genus"])
#Generate table using genfac as index and pulling values from the otu_table
gentab = apply(otu_table(phy_tax_glom_gen_trimmed_normalize), MARGIN = 2, function(x) {
     tapply(x, INDEX = genfac, FUN = sum, na.rm = TRUE, simplify = TRUE)
 })
```

*Heatmap experimentation*

``` r
#Generate basic heatmap (messy)
pheatmap(gentab)
```

![Figure1](link-to-image)

*Create row
annotations*

``` r
#Use psmelt function to convert phyloseq object into single, large R object
melted_phy_table = psmelt(phy_tax_glom_gen_trimmed_normalize)
#Reorder based on genus name
reordered_melted = melted_phy_table[order(melted_phy_table$Genus),]
#Remove duplicate genus names
reordered_melted_unique = reordered_melted %>% distinct(Genus, .keep_all = T)
#Create data frame with phylum and genus names paired 
phyla_frame = subset(reordered_melted_unique, select = c(Phylum, Genus))
#Set rownames to the genus names
rownames(phyla_frame) = phyla_frame$Genus
#Remove column we just used as rownames (no longer needed)
phyla_frame = subset(phyla_frame, select = -c(Genus))
```

*Create column annotations*

``` r
#Only keep unique samples within the melted R object created above
reordered_melted_unique_sample = reordered_melted %>% distinct(Sample_ID, .keep_all = T)
#Create data frame with Sample ID and NAFLD status
NAFLD_frame = subset(reordered_melted_unique_sample, select = c(NAFLD, Sample_ID))
#Set rownames to SampleID
rownames(NAFLD_frame) = NAFLD_frame$Sample_ID
#Remove column we just used as rownames (no longer needed)
NAFLD_frame = subset(NAFLD_frame, select = -c(Sample_ID))
```

*Produce heatmap with both rows and columns clustered hierarchically but
annotated based on phylum and NAFLD status,
respectively.*

``` r
pheatmap(log10(gentab+1), annotation_row = phyla_frame, annotation_col = NAFLD_frame, angle_col = 315, fontsize_row = 10, fontsize_col = 8, show_colnames = F, border_color = NA)
```

![](rmarkdown_for_github_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

*Alphabetize phyla and create annotation for alphabetical Phylum-Genus
pairs*

``` r
#Order phyla within melted object  alphabetically. The "order" function defaults to this automatically
reordered_melted_phyla = melted_phy_table[order(melted_phy_table$Phylum),]
#Remove duplicate genus entries
reordered_melted_phyla_no_genus_dups = reordered_melted_phyla %>% distinct(Genus, .keep_all = TRUE)
#Reorder gentab based on genera ordered by phylum, alphabetically
gentab_reorder = gentab[reordered_melted_phyla_no_genus_dups$Genus, order(colnames(gentab))]
```

*Regenerate phyla frame (for row annotations) using reordered objects
from
above*

``` r
phyla_frame = data.frame(Phylum = reordered_melted_phyla_no_genus_dups$Phylum)
rownames(phyla_frame) = reordered_melted_phyla_no_genus_dups$Genus
```

*Group Samples by NAFLD Status (rather than hierarchal)*

``` r
#Order rows and columns of melted object by NAFLD status and Sample ID
reordered_melted_NAFLD = melted_phy_table[order(melted_phy_table$NAFLD, melted_phy_table$Sample_ID),]
#Remove duplicate Sample IDs
reordered_melted_NAFLD_no_ID_dups = reordered_melted_NAFLD %>% distinct(Sample_ID, .keep_all = TRUE)
#Reorder gentab by Sample ID
#Received "subscript out of bounds" error with this method, attempting to troubleshoot
#*******FIXED********
#Had to treat column from melted phy table "as.character" to allow appropriate ordering to occur
gentab_reorder_ID_final = gentab_reorder[, as.character(reordered_melted_NAFLD_no_ID_dups$Sample_ID)]
```

*Generate Heatmap using pheatmap function*

``` r
#Custom colors for annotations
ann_colors = list(
    NAFLD = c(yes = "red", no = "green"))
#First argument is inputing finalized gentab object. This has undergone a log transformation (log10) and we have added a pseudocount (+1) to account for zero column entries. Cluster columns and rows is turned off as we are manually grouping taxa and samples. Phyla frame is used to annotate the genera by phylum. NAFLD frame is used to annotate samples by NAFLD status. Other arguments are aesthetic. 
pheatmap(log10(gentab_reorder_ID_final+1), cluster_rows = FALSE, cluster_cols = FALSE, annotation_row = phyla_frame, annotation_col = NAFLD_frame, angle_col = 315, fontsize_row = 10, fontsize_col = 8, show_colnames = F, gaps_col = 29, border_color = NA)
```

![](rmarkdown_for_github_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

``` r
sessionInfo()
```

    ## R version 3.5.2 (2018-12-20)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 10 x64 (build 18362)
    ## 
    ## Matrix products: default
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_United States.1252 
    ## [2] LC_CTYPE=English_United States.1252   
    ## [3] LC_MONETARY=English_United States.1252
    ## [4] LC_NUMERIC=C                          
    ## [5] LC_TIME=English_United States.1252    
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] pheatmap_1.0.12     viridis_0.5.1       viridisLite_0.3.0  
    ##  [4] forcats_0.4.0       stringr_1.4.0       dplyr_0.8.3        
    ##  [7] purrr_0.3.3         readr_1.3.1         tidyr_1.0.2        
    ## [10] tibble_2.1.3        ggplot2_3.2.1       tidyverse_1.3.0    
    ## [13] qiime2R_0.99.11     BiocManager_1.30.10 phyloseq_1.26.1    
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] nlme_3.1-137        fs_1.3.1            lubridate_1.7.4    
    ##  [4] RColorBrewer_1.1-2  httr_1.4.1          tools_3.5.2        
    ##  [7] backports_1.1.5     R6_2.4.1            vegan_2.5-6        
    ## [10] rpart_4.1-15        Hmisc_4.3-0         DBI_1.1.0          
    ## [13] lazyeval_0.2.2      BiocGenerics_0.28.0 mgcv_1.8-31        
    ## [16] colorspace_1.4-1    permute_0.9-5       ade4_1.7-13        
    ## [19] nnet_7.3-12         withr_2.1.2         tidyselect_0.2.5   
    ## [22] gridExtra_2.3       compiler_3.5.2      cli_2.0.1          
    ## [25] rvest_0.3.5         Biobase_2.42.0      pacman_0.5.1       
    ## [28] htmlTable_1.13.3    xml2_1.2.2          scales_1.1.0       
    ## [31] checkmate_1.9.4     digest_0.6.23       foreign_0.8-71     
    ## [34] rmarkdown_2.1       XVector_0.22.0      base64enc_0.1-3    
    ## [37] pkgconfig_2.0.3     htmltools_0.4.0     dbplyr_1.4.2       
    ## [40] htmlwidgets_1.5.1   rlang_0.4.3         readxl_1.3.1       
    ## [43] rstudioapi_0.10     farver_2.0.3        generics_0.0.2     
    ## [46] jsonlite_1.6        acepack_1.4.1       magrittr_1.5       
    ## [49] Formula_1.2-3       biomformat_1.10.1   Matrix_1.2-18      
    ## [52] fansi_0.4.1         Rcpp_1.0.3          munsell_0.5.0      
    ## [55] S4Vectors_0.20.1    Rhdf5lib_1.4.3      ape_5.3            
    ## [58] lifecycle_0.1.0     stringi_1.4.5       yaml_2.2.0         
    ## [61] MASS_7.3-51.5       zlibbioc_1.28.0     rhdf5_2.26.2       
    ## [64] plyr_1.8.5          grid_3.5.2          parallel_3.5.2     
    ## [67] crayon_1.3.4        lattice_0.20-38     Biostrings_2.50.2  
    ## [70] haven_2.2.0         splines_3.5.2       multtest_2.38.0    
    ## [73] hms_0.5.3           knitr_1.27          pillar_1.4.3       
    ## [76] igraph_1.2.4.2      reshape2_1.4.3      codetools_0.2-16   
    ## [79] stats4_3.5.2        reprex_0.3.0        glue_1.3.1         
    ## [82] evaluate_0.14       latticeExtra_0.6-28 data.table_1.12.8  
    ## [85] modelr_0.1.5        vctrs_0.2.2         foreach_1.4.7      
    ## [88] cellranger_1.1.0    gtable_0.3.0        assertthat_0.2.1   
    ## [91] xfun_0.12           broom_0.5.3         survival_3.1-8     
    ## [94] iterators_1.0.12    IRanges_2.16.0      cluster_2.1.0
