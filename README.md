# Intragenomic_analysis

This repository contains the code for the experiments and the CGR-Diff software proposed in our paper.

![Experiments](Figures/Experiments.jpg)

## Dataset

The dataset used in our paper includes genomes from various kingdoms, organized into three subsets based on the type of
analysis performed: Reference (Subset 1), Intergenomic Analysis (Subset 2), and Intragenomic Analysis (Subset 3). Below
is a summary of the datasets used in our experiments, along with links to their corresponding assemblies in the NCBI
database.

<table>
  <thead>
    <tr>
      <th>Dataset</th>
      <th>Kingdom</th>
      <th>Species (Common Name)</th>
      <th>Assembly Link</th>
      <th>Length (Mbp)</th>
      <th>% N</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td rowspan="1"><em>Reference (Subset 1)</em></td>
      <td rowspan="1">Animalia</td>
      <td><em>Homo sapiens</em> (human)</td>
      <td><a href="https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_009914755.4/">GCA_009914755.4</a></td>
      <td>3117</td>
      <td>0</td>
    </tr>
    <tr>
      <td rowspan="8"><em>Intergenomic Analysis (Subset 2)</em></td>
      <td rowspan="3">Animalia</td>
      <td><em>Pan troglodytes</em> (chimpanzee)</td>
      <td><a href="https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_028858775.2/">GCA_028858775.2</a></td>
      <td>3178</td>
      <td>0.16</td>
    </tr>
    <tr>
      <td><em>Mus musculus</em> (house mouse)</td>
      <td><a href="https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000001635.9/">GCA_000001635.9</a></td>
      <td>2723</td>
      <td>2.7</td>
    </tr>
    <tr>
      <td><em>Drosophila melanogaster</em> (fruit fly)</td>
      <td><a href="https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000001215.4/">GCA_000001215.4</a></td>
      <td>80</td>
      <td>0.57</td>
    </tr>
    <tr>
      <td rowspan="1">Fungi</td>
      <td><em>Saccharomyces cerevisiae</em> (yeast)</td>
      <td><a href="https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000146045.2/">GCA_000146045.2</a></td>
      <td>12</td>
      <td>0</td>
    </tr>
    <tr>
      <td rowspan="1">Plantae</td>
      <td><em>Arabidopsis thaliana</em> (thale cress)</td>
      <td><a href="https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000001735.2/">GCA_000001735.2</a></td>
      <td>119</td>
      <td>0.16</td>
    </tr>
    <tr>
      <td rowspan="1">Protista</td>
      <td><em>Paramecium caudatum</em> 📌</td>
      <td><a href="https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000715435.1/">GCA_000715435.1</a></td>
      <td>30</td>
      <td>2.16</td>
    </tr>
    <tr>
      <td rowspan="1">Archaea</td>
      <td><em>Pyrococcus furiosus</em></td>
      <td><a href="https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_008245085.1/">GCA_008245085.1</a></td>
      <td>2</td>
      <td>0</td>
    </tr>
    <tr>
      <td rowspan="1">Bacteria</td>
      <td><em>Escherichia coli</em></td>
      <td><a href="https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000005845.2/">GCA_000005845.2</a></td>
      <td>5</td>
      <td>0</td>
    </tr>
    <tr>
      <td rowspan="3"><em>Intragenomic Analysis (Subset 3)</em></td>
      <td rowspan="1">Fungi</td>
      <td><em>Aspergillus nidulans</em></td>
      <td><a href="https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000011425.1/">GCA_000011425.1</a></td>
      <td>30</td>
      <td>0.04</td>
    </tr>
    <tr>
      <td rowspan="1">Plantae</td>
      <td><em>Zea mays</em> (maize)</td>
      <td><a href="https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_022117705.1/">GCA_022117705.1</a></td>
      <td>2179</td>
      <td>0</td>
    </tr>
    <tr>
      <td rowspan="1">Protista</td>
      <td><em>Dictyostelium discoideum</em></td>
      <td><a href="https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000004695.1/">GCA_000004695.1</a></td>
      <td>34</td>
      <td>0.07</td>
    </tr>
  </tbody>
</table>

<p>📌 Among all the species in this table, the assembly for <em>Paramecium caudatum</em> is at the scaffold level.</p>

## Run the experiments in our paper

1. Clone this repository and install the required libraries by running:
```bash
# Clone the repository
git clone https://github.com/Niousha12/Intragenomic_analysis.git
cd Intragenomic_analysis

# Create a virtual environment and install requirements
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```
2. Download the assemblies for each species from the NCBI website, and organize the data in the `Data` folder to match the following structure:
```
📂 Data/
├── 📂 Human/
│   ├── 📂 chromosomes/
│   │   ├── chr1.fna
│   │   └── ...
│   └── 📂 bedfiles/
│       ├── cytobands.bed
│       ├── telomere.bed
│       └── centromere.bed
├── 📂 Chimp/
│   └── 📂 chromosomes/
│       ├── chr1.fna
│       └── ...
├── 📂 Mouse/
│   └── 📂 chromosomes/
│       └── ....fna
├── 📂 Drosophila melanogaster/
│   └── 📂 chromosomes/
│       └── ....fna
├── 📂 Saccharomyces cerevisiae/
│   └── 📂 chromosomes/
│       └── ....fna
├── 📂 Arabidopsis thaliana/
│   └── 📂 chromosomes/
│       └── ....fna
├── 📂 Paramecium caudatum/
│   └── 📂 chromosomes/
│       └── ....fna (These are scaffolds)
├── 📂 Pyrococcus furiosus/
│   └── 📂 chromosomes/
│       └── ....fna
├── 📂 Escherichia coli/
│   └── 📂 chromosomes/
│       └── ....fna
├── 📂 Aspergillus nidulans/
│   └── 📂 chromosomes/
│       └── ....fna
├── 📂 Maize/
│   └── 📂 chromosomes/
│       └── ....fna
└── 📂 Dictyostelium discoideum/
    └── 📂 chromosomes/
        └── ....fna
```

Alternatively for `Human`, `Chimp`, `Mouse`, and `Maize` you can use the following command to download the chromosomes assemblies directly from the NCBI FTP server.

```bash
# Homo sapiens (human)
wget -P Data/Human/chromosomes/ -r -nH --cut-dirs=12 --no-parent \
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/914/755/GCA_009914755.4_T2T-CHM13v2.0/GCA_009914755.4_T2T-CHM13v2.0_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/
gunzip Data/Human/chromosomes/*.gz

# Pan troglodytes (chimpanzee)
wget -P Data/Chimp/chromosomes/ -r -nH --cut-dirs=12 --no-parent \
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/028/858/775/GCA_028858775.2_NHGRI_mPanTro3-v2.0_pri/GCA_028858775.2_NHGRI_mPanTro3-v2.0_pri_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/
gunzip Data/Chimp/chromosomes/*.gz

# Mus musculus (house mouse)
wget -P Data/Mouse/chromosomes/ -r -nH --cut-dirs=12 --no-parent \
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/635/GCA_000001635.9_GRCm39/GCA_000001635.9_GRCm39_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/
gunzip Data/Mouse/chromosomes/*.gz

# Zea mays (maize)
wget -P Data/Maize/chromosomes/ -r -nH --cut-dirs=12 --no-parent \
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/022/117/705/GCA_022117705.1_Zm-Mo17-REFERENCE-CAU-T2T-assembly/GCA_022117705.1_Zm-Mo17-REFERENCE-CAU-T2T-assembly_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/
gunzip Data/Maize/chromosomes/*.gz
```

You can also download the complete datasets used in the paper from the [Google Drive](https://drive.google.com/file/d/1q7fbymvlAd7XLA7D94QN575tON1qk1fR/view?usp=sharing).


