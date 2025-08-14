# Future work



Import counts and normalized counts into repo.&#x20;



DRB -> Radiation counts: Bourbousse et al., -> RawExpressionValue



### Research Plan: ROS Signature Identification Tool using Machine Learning

**1. Introduction**

Reactive Oxygen Species (ROS) play a crucial role in plant stress responses. Accurately identifying ROS signatures is vital for understanding plant adaptation to various environmental challenges, including those encountered in space. This research aims to develop a novel ROS signature identification tool using machine learning (ML) techniques.

**2. Background**

* Existing tools like ROS-o-meter rely on microarray data, limiting their scope.
* Recent research, such as the ROS wheel paper, utilizes a broader palette of ROS-related microarrays for ROS cluster identification.

**3. Research Objectives**

* Build on the ROS wheel concept and develop a ROS signature identification tool using RNA-seq data.
* Compile RNA-seq expression matrices and associated metadata from a variety of ROS-related studies focusing on leaves and roots.
* Employ machine learning algorithms to identify markers for specific ROS types from the compiled data.

**4. Methodology**

4.1 Data Acquisition:

* Leverage pre-processed RNA-seq datasets from NCBI's Sequence Read Archive (SRA) and resources like [http://bioinformatics.sdstate.edu/reads/](http://bioinformatics.sdstate.edu/reads/).
* Utilize search terms like "H2O2," "ROS," "Reactive oxygen species," and specific ROS types to identify relevant studies.
* Focus on acquiring at least 3 RNA-seq datasets for leaves and roots exposed to different ROS types.

4.2 Data Processing:

* Merge the acquired RNA-seq expression matrices and associated metadata (factors) matrices, potentially involving manual curation to ensure completeness.
* Factors of interest include:
  * Organ/tissue (leaf, root)
  * Age/developmental stage
  * Treatment (ROS type and intensity)

4.3 Machine Learning Approach:

* Employ supervised machine learning algorithms to identify gene expression patterns or markers indicative of specific ROS types based on the compiled data.
* Potential algorithms include:
  * Support Vector Machines (SVM)
  * Random Forests
  * Deep Learning techniques

**5. Expected Outcomes**

* Develop a novel ROS signature identification tool using an RNA-seq based approach.
* Identify robust gene expression markers associated with different ROS types in leaves and roots.
* Enhance our understanding of plant responses to ROS stress in space environments.





Link to table of studies selected for use in a ROS recognition system.&#x20;



{% embed url="https://docs.google.com/spreadsheets/d/10h0bL03afgNXbeM4uxJv9qzzwIlagTu1Mt3h-Fi7kmw/edit?usp=sharing" %}

