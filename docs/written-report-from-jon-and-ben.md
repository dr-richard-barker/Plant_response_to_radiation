# Report from Jon and Ben

**Introduction/background**

Exploration beyond Earth’s atmosphere introduces a host of novel challenges that threaten normal biological processes in living organisms. Specifically, space missions expose humans to a radiation environment of a composition and intensity not encountered within our biosphere (Horneck, 1999). Challenging the normal functioning within living organisms presents the potential for decreased survivability, particularly in the context of high intensity radiation. Relatively common on Earth, UV radiation is known to be a potent agent for the induction of programmed cell death in human skin (Martin et al., 1991) At an average frequency of 10^16 Hz, the intensity of UV radiation is roughly one hundred millionths that of typical cosmic rays; consequently, cosmic rays may amplify the effects observed in UV radiation experimentation as well as pose novel challenges for biological life. Therefore, the inherent presence of extreme high intensity radiation in outer space, like that of cosmic rays demands a fundamental understanding of its effects on biological systems.

Many studies have already examined how the transcriptome of Arabidopsis thaliana is altered in response to sublethal doses of ionizing radiation. Culligan et al., 2006 \[1], subjected five day old WS-2 seedlings to 100 Gy gamma radiation, and measured RNA expression levels 1.5 after exposure through microarray analysis. Missirian et al., 2014 \[2] subjected WS-2 five-day old seedlings to both 30 Gy high-energy Fe nuclei (HZE) and 100 Gy of gamma (GW) ionizing radiation, and assessed gene expression changes through a time-series analysis of microarray data (1.5 hours to 24 hours post ionizing radiation treatment). Differentially expressed genes were identified (p < 0.05) and genes predicted to be associated with regulating circadian rhythm were removed. However, very few studies have examined the transcriptomics of A. thaliana in the context of ionizing radiation in a spaceflight environment where variables prove more difficult to control. To identify potential radiation response candidates that respond to cosmic radiation, several spaceflight datasets were chosen to compare against the aforementioned literature. Fengler et al., 2014 \[3] gathered microarray data from callus Col-0 cultures grown within a 1.0g centrifuge aboard the Shenzhou 8 rocket for 5 days. By controlling for the factor of microgravity through use of a 1.0g centrifuge, differentially expressed genes are more likely to be involved in the radiation response than responses to other unique characteristics of spaceflight. Choi et al., 2019 \[4] report RNA-seq data from Col-0 seedlings grown within BRIC19 aboard the ISS. Ionizing radiation can directly lead to genomic instability in the form of double strand breaks and indirectly through oxidative damage induced by ROS production (Woodmine et al., 2011). In the context of ROS, Willems et al., 2016 \[5] defined “high light early” (HLE) and late (HLL) responding clades of genes in high frequency light conditions. Availability of these transcriptomic experiments and their results provides opportunity for a meta-analysis identifying putative radiation-response genes.

\


**Methods**

Meta-analysis

Transcriptomic gene expression data was gathered from the experiments mentioned above for a bioinformatics-based meta-analysis. All significantly differentially expressed (p<0.05) genes between experiments 1-5 were compared in order to identify shared differentially expressed genes across multiple studies. Two lists were generated in parallel using genes in either HZE bombardment or GW radiation datasets obtained from experiment 2. The original HZE list is composed of genes common across experiments 2 (HZE bombardment), 3, 4 and 5. In the case of experiment 5, differentially expressed genes from HLE and HLL were combined as temporal ROS response is not of particular concern, thereby reducing the unnecessary filtering of false negatives. An amended HZE list adds to its comparison genes from experiment 1 because the original HZE list yields an unmanageable amount of genes for further analysis. The GW list is composed of genes common across experiments 1, 2 (GW radiation), 3, 4 and 5. The same genes were included from experiment 5 as the HZE list. Finally, a final list was constructed by identifying common genes between the amended HZE list and GW list.&#x20;

Circadian rhythm-filtered time series microarray data from experiment 2 was then analyzed with DREM2.0 software in order to assign shared gene expression-paths with known transcription factor interactions. Transcription factors were identified in certain path clusters, and those transcription factors were used to select genes of interest; specifically, any gene with known interactions with selected transcription factors as defined by DREM2.0 was grouped to that transcription factor’s gene list. This analysis was done for gene expression in both HZE and gamma-irradiated environments. DREM2.0 settings used include: no normalization/add 0, maximum number of missing values is 4, minimum absolute expression change is 1.&#x20;

\


Cloning of HY5 promoter into pCAMBIA1302

One of the transcription factors identified by DREM2.0 for both HZE and GW treatments in the missirian et al. 2014 dataset was selected as a candidate for a transcriptional fusion. The 1 Kb upstream promoter sequence of HY5 was obtained using thalemine and was cloned into a modified pCAMBIA1302 vector (pUBQ10::mGFP) using the NEBuilder primer design software (Table 1). The 1 Kb upstream promoter sequence of HY5 was amplified using the Q5 Hotstart 2x Mastermix following the recommended guidelines using gDNA extracted from Col-0 (denaturing: 10 sec, annealing: 30 sec, extension 30 sec). PCR products were then purified using column-cleanup (QIAGEN) prior to ligation into the backbone.&#x20;

\


To generate 5’ overhangs, EcoRI and SalI were used to create a double-restriction digest in the backbone of the vector (Table 2).  The restriction enzymes were then heat inactivated at 65°C. 70 ng of digested vector were then combined with 20.5 ng of insert (1:2 ratio) and combined with the 10x NEBuilder assembly mastermix and was incubated for four hours at 50°C to improve ligation efficiency. Gibson assembly procedure is illustrated in Figure 3. Ligation products were then transformed into chemically competent DH5α E. coli and selected using Kanamycin (1000X) plates. The DNA was extracted from putative positive clones using DNA minipreps and were confirmed to be ligated correctly using Benchlings alignment software (Figure 4) and were sequenced with a manually designed sequencing primer (Table 1).&#x20;

\


Table 1. Primer design output obtained from NEBuilder (settings: min. overlap: 20, primer concentration: 500 nM, min primer length: 18). Gene-specific regions are designated with bolded characters

| Primer                    | Sequence                                         | Tm    | Ta    |
| ------------------------- | ------------------------------------------------ | ----- | ----- |
| pHY5 forward Primer       | ggaaacagctatgaccatgattacgGTTTGTCTTGACAGGGCTTATG  | 63.5° | 64.5° |
| pHY5 reverse Primer       | tcctttactagtcagatctaccatgTGGAACGTGGGTGGATTTTAG   | 64.7° | 64.5° |
| Forward Sequencing Primer | GCTTCCGGCTCGTATGTTGTG                            | NA    | 58.3° |

\
\


**Table 2.** Protocol for digesting pCAMBIA 1302

| Reagent           | Volume (uL) |
| ----------------- | ----------- |
| Vector (100ng/uL) | 30          |
| SalI              | 1           |
| EcoR1             | 1           |
| Cutsmart (10x)    | 5           |
| ddH2O             | 13          |
| Total             | 50          |
| Temperature       | 37°         |

\
**Nebuilder primers**



* &#x20;HY5\_1Kb\_Upstream\_promoter\_fwd ggaaacagctatgaccatgattacgGTTTGTCTTGACAGGGCTTATG 3'Tm=63.5 3'Ta(annealing temp)=64.5
* HY5\_1Kb\_Upstream\_promoter\_rev tcctttactagtcagatctaccatgTGGAACGTGGGTGGATTTTAG 3'Tm=64.7 3'Ta(annealing temp)=64.5
* PCR settings (Ta, hotstart polymerase recommended)
* 1kb upstream promoter obtained from thalemine hy5 webpage
* pCAMBIA1302: modified (GFP)
* Position relative to cutsite: Restriction enzyme&#x20;
* 5’ upstream: EcoR1
* 3’ downstream: Sal1

\
\
\


.&#x20;

\
\


**Results:**&#x20;

Both the HZE and GW timeseries microarray data was analyzed using DREM2.0 (no normalization, max. missing values: 4, default paramters).&#x20;

\
\
\
\
\
\
\
\
\
\


<figure><img src="https://lh7-us.googleusercontent.com/TUKmPHtGonor9ZBtpFvmRg40bo_8Olm1f75juferqy7CphC9T2O8HrBQx6ZZ0bICaC2XMgy1Lt6dDSaqaVobku6l7ZGolAsRWf1fTGNMPyH5rptuTJSc55b60dhRUvb6k3ds7LYCW1JOMd0U7CdJzw" alt=""><figcaption></figcaption></figure>

\
\


<figure><img src="https://lh7-us.googleusercontent.com/eqlP2mBvjOcwF8ACDLYd-QiUqH3mBAtv37VyS-owAAfFnkdKa0LwgdIW9XnWGqnmv_IeDpbq2w1XE0eUHO-3Nwekwquh5xZdaaDKVMyn4PqlflHdDlGwuRCMvwZHnUnPMT-W4XfsuyqQAL57Dqcvkg" alt=""><figcaption></figcaption></figure>

\


<figure><img src="https://lh7-us.googleusercontent.com/I0nTRSKjp-UKbc-e0iJdxEZYQO7iL3N53QTxNsLFGTauNvrK8wk8NokpTssjyHyevyU0uSN3e2wXD6CB8CzefX2G_Mx_fYuzflp7OfkC4JY81J3aC10oPhfrlvEwLskXJwQe-1a00j9oQlccWSEj2g" alt=""><figcaption></figcaption></figure>

Figure X: Heatmap depicting the eight radiation response candidates shared between all datasets in either HZE or GW based radiation. Genes and datasets were hierarchically clustered using a Pearson (1 - S) distance-based approach, and UPGMA. Genes of biological interest are outlined in black.&#x20;

\
\
\
\
\


**Discussion**

\


To our knowledge, this is the first meta-analysis that has attempted to identify genes in Arabidopsis thaliana. that are differentially expressed in response to prolonged exposure to ionizing radiation experienced during spaceflight. In total, five separate experimental studies were used to analyze differentially expressed shared between spaceflight and two distinct forms of ionizing radiation (HZE,GW). To account for differences between possible differences of HZE and GW radiation, two parallel analyses were conducted for each type of ionizing radiation. In total eight radiation response candidates were common between all five experimental studies. Of these eight, four genes were chosen for genes of biological interest (SIP4, SIP3, WRKY25, BCB) were selected due to similar patterns of gene expression across multiple datasets. Furthermore, three transcription factors (HY5, AGL15,AGL9) were enriched in the time series dataset for both HZE and GW data reported in Missirian et al. 2014. One such transcription factor HY5, was of biological interest to us, due to known associations with UV-B response activity, and its significant upregulation in response to GW reported in Missirian et al. 2014. For these reasons, HY5 was selected for a transcriptional fusion using a gibson-based assembly protocol.&#x20;

Overall, the gene ontology of the eight shared radiation response candidates yielded interesting and unexpected results. Cellular responses to temperature stimuli were highly represented within these eight candidates, suggesting that acute exposure to ionizing radiation may affect the composition of membrane lipids, either directly, or indirectly through oxidative damage induced by the generation of ROS. Likewise, several candidate genes have known associations to processes that are connected to pathways associated with the generation, regulation or response to reactive oxygen species. Additionally, a number of genes had connections to salt and osmotic stress.&#x20;

\
\


Justification of Pil5/Hy5 filtering as opposed to AGL15 Abridged:

\


**Pil5/Hy5**: Include the following in the justification:

\-Both Pil5 and Hy5 lack a clear and noticeable trend and fluctuate frequently unlike AGL15’s&#x20;

\-Possible missiran artefact with the NaN reported values&#x20;

\-Hard to detect this pattern using transcriptional fusions to GFP

\


Therefore only genes known to have associations with these transcription factors were selected to avoid generating false positives.

\
\
\


**AGL15:** The genes that pass through AGL15’s 24 hour node are noticeably suppressed at later time points. Considering the small standard deviation of the genes grouped at this node, the results suggest that AGL15 may play both a direct and indirect role. While some of these genes may not be directly bound by AGL15, the results suggest that the genes that pass through AGL15’s 24 hour node may be secondary downstream targets of AGL15 mediated transcriptional regulation.&#x20;

\


This pattern can be feasibly identified through transcriptional fusions and would provide good markers for detecting responses to radiation.

\
**Ignore all of the following:** Microarray data of Arabidopsis thaliana that is relevant in the context of ionizing radiation already exists in the literature. Fengler et al., 2014 \[1] gathered data from callus cultures grown within a 1.0g centrifuge aboard the Shenzhou 8 rocket for 5 days. Published data from seedlings in BRIC19 aboard the ISS also provides spaceflight gene expression data \[2]. Missirian et al., 2014 \[3] subjected seedlings to ionizing radiation in the form of either gamma radiation or HZE particles, and reported time-series data filtered for circadian rhythm. Culligan et al., 2006 \[4] subjected seedlings to gamma radiation as well. Ionizing radiation has well-documented effects on causing double stranded breaks in DNA. One mechanism by which these double stranded breaks may be caused is reactive oxygen species (ROS) generated in the cell under ionizing radiation conditions. In the context of ROS, Willems et al., 2016 \[5] defined “high light early” (HLE) and late (HLL) responding clusters of genes in high light environments. The availability of experiments like these provides opportunity for a meta-analysis identifying putative radiation-response genes. &#x20;

**JL Edit:** Many studies have already examined how the transcriptome of Arabidopsis thaliana is altered in response to sublethal doses of ionizing radiation.  Similarly, Culligan et al., 2006 \[2], subjected five day old WS-2 seedlings to 100 Gy gamma radiation, and measured RNA expression levels 1.5 after exposure. Missirian et al., 2014 \[1] subjected WS-2 five-day old seedlings to both 30 Gy high-energy Fe nuclei (HZE) and 100 Gy of gamma (GW) ionizing radiation, and assessed gene expression changes through a time-series analysis of microarray data (1.5 hours to 24 hours post ionizing radiation treatment). Differentially expressed genes were identified (p < 0.05) and genes predicted to be associated with regulating circadian rhythm were removed. However, very few studies have examined the effects of ionizing radiation in the context of genes that respond to the spaceflight environment. Yet it can be difficult to differentiate genes that are responding to certain aspects of the spaceflight environment, such as microgravity, to possible biological responses associated with ionizing radiation. To identify potential radiation response candidates that respond to cosmic radiation, several spaceflight datasets were chosen to compare against the aforementioned literature.&#x20;

\
\
\


Fengler et al., 2014 \[1] gathered data from callus cultures grown within a 1.0g centrifuge aboard the Shenzhou 8 rocket for 5 days. Published data from seedlings in BRIC19 aboard the ISS also provides spaceflight gene expression data \[2]. Missirian et al., 2014 \[3] subjected seedlings to ionizing radiation in the form of either gamma radiation or HZE particles, and reported time-series data filtered for circadian rhythm. Culligan et al., 2006 \[4] subjected seedlings to gamma radiation as well. Ionizing radiation has well-documented effects on causing double stranded breaks in DNA. One mechanism by which these double stranded breaks may be caused is reactive oxygen species (ROS) generated in the cell under ionizing radiation conditions. In the context of ROS, Willems et al., 2016 \[5] defined “high light early” (HLE) and late (HLL) responding clusters of genes in high light environments. Availability of these transcriptomic experiments and their results provides opportunity for a meta-analysis identifying putative radiation-response genes.

\


The gamma radiation identified candidates

\


&#x20;

\


the parallel comparisons were crossed with the spaceflight experiments (\_\_,\_\_) to generate two parallel lists shared between all datasets Figure 1. Finally, a third list was constructed by identifying the genes shared between the two parallel radiation analyses, resulting in eight radiation response candidates that are significantly differentially expressed across all studies and types of radiation.&#x20;

\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\


<figure><img src="https://lh7-us.googleusercontent.com/uViHAQ6gJaN2GiPqoggcMsXAKM--SGQwjqkip6U_61ZlfDU4sag3XToXYbFSV5W-S9ufg9cbBx61wtdbo6ml1vNdkC8koCtGy4Y1ElD4dzu6LwcE4UOR-zhwN0In9vG-SJhQm198vW-_v6oMu2oUmA" alt=""><figcaption></figcaption></figure>
