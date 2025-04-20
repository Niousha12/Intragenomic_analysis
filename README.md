# Intragenomic_analysis

This repository contains the code for the experiments and the CGR-Diff software proposed in our paper.

![Experiments](Figures/Experiments.jpg)

## Dataset Summary

The dataset used in our paper includes genomes from various kingdoms, organized into three subsets based on the type of analysis performed: Reference (Subset 1), Intergenomic Analysis (Subset 2), and Intragenomic Analysis (Subset 3). Below is a summary of the datasets used in our experiments, along with links to their corresponding assemblies in the NCBI database.

| Dataset                          | Kingdom   | Species (Common Name)                    | Assembly Link | Length (Mbp) | % N   |
|----------------------------------|-----------|-------------------------------------------|----------------|--------------|--------|
| *Reference* (Subset 1)           | Animalia  | *Homo sapiens* (human)                   | [GCA_009914755.4](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_009914755.4/) | 3117         | 0      |
| *Intergenomic Analysis* (Subset 2) | Animalia  | *Pan troglodytes* (chimpanzee)           | [GCA_028858775.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_028858775.2/) | 3178         | 0.16   |
|                                  |           | *Mus musculus* (house mouse)             | [GCA_000001635.9](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000001635.9/) | 2723         | 2.7    |
|                                  |           | *Drosophila melanogaster* (fruit fly)    | [GCA_000001215.4](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000001215.4/) | 80           | 0.57   |
|                                  | Fungi     | *Saccharomyces cerevisiae* (yeast)       | [GCA_000146045.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000146045.2/) | 12           | 0      |
|                                  | Plantae   | *Arabidopsis thaliana* (thale cress)     | [GCA_000001735.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000001735.2/) | 119          | 0.16   |
|                                  | Protista  | *Paramecium caudatum*[^1]                | [GCA_000715435.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000715435.1/) | 30           | 2.16   |
|                                  | Archaea   | *Pyrococcus furiosus*                    | [GCA_008245085.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_008245085.1/) | 2            | 0      |
|                                  | Bacteria  | *Escherichia coli*                       | [GCA_000005845.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000005845.2/) | 5            | 0      |
| *Intragenomic Analysis* (Subset 3) | Fungi     | *Aspergillus nidulans*                   | [GCA_000011425.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000011425.1/) | 30           | 0.04   |
|                                  | Plantae   | *Zea mays* (maize)                       | [GCA_022117705.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_022117705.1/) | 2179         | 0      |
|                                  | Protista  | *Dictyostelium discoideum*               | [GCA_000004695.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000004695.1/) | 34           | 0.07   |

[^1]: Among all the species in this table, the assembly for *Paramecium caudatum* is at the scaffold level.


