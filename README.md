# Astrocyte-metabolic-modeling-MCI-AD
This repository contains the code and data supporting the study “Transcriptome-Informed Metabolic Modeling Reveals Astrocyte-Specific Vulnerabilities in Mild Cognitive Impairment and Alzheimer’s Disease Progression.” The project focuses on the integration of transcriptomic data into genome-scale metabolic models (GEMs) to characterize astrocyte-specific metabolic alterations across disease progression.

In this work, gene expression data is used to construct and simulate condition-specific metabolic models of human astrocytes under four biological states: Control, Incipient Mild Cognitive Impairment (MCI), Moderate MCI, and Severe MCI/Alzheimer’s disease. Through this approach, we aim to capture the progressive metabolic rewiring that occurs during neurodegeneration and identify potential vulnerabilities associated with each stage.

The analytical framework follows a systems biology approach in which transcriptomic data is first integrated into a generic astrocyte metabolic network to generate context-specific models. These models are subsequently constrained and analyzed using flux balance analysis and related methods, enabling the comparison of metabolic flux distributions across conditions. The resulting outputs allow for the identification of pathway-level alterations and stage-dependent metabolic changes.

The repository includes RNA-seq datasets corresponding to each condition, four condition-specific genome-scale metabolic models, and the scripts required to perform data integration and metabolic simulations. The structure of the repository is organized to facilitate reproducibility and reuse, allowing users to follow the workflow from data input to final results.

To reproduce the analysis, users should clone the repository, open it in MATLAB, and ensure that the COBRA Toolbox is properly installed. The scripts provided in the repository can then be executed to reconstruct the models and perform the simulations described in the study.

It is important to note that genome-scale metabolic modeling and flux-based analyses rely on well-established computational frameworks that are predominantly implemented in MATLAB, particularly through the COBRA Toolbox. For this reason, MATLAB was used as the primary environment in this project, ensuring methodological robustness and compatibility with standard practices in systems biology.

From a biological perspective, this work highlights the central role of astrocytes in maintaining brain metabolic homeostasis and provides insights into how their metabolism is progressively altered during the transition from MCI to Alzheimer’s disease. The results suggest a dynamic reorganization of metabolic processes, reflecting early adaptive responses followed by increasing metabolic dysregulation at later stages.

All data and scripts are made available to ensure full reproducibility and to support further development by the research community. For questions or potential collaborations, please contact Maria Andrea Angarita Rodríguez (maria.angaritar@javeriana.edu.co
), Doctoral Candidate in Biological Sciences and MSc (c) in Bioinformatics.

Maria Andrea Angarita Rodríguez
Doctoral Candidate in Biological Sciences
MSc (c) in Bioinformatics
maria.angaritar@javeriana.edu.co

