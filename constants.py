HETERO_HETERO_DICT = {
    "1": ["p31.1", "p21.1", "q31.1", "q31.3", "q41"],
    "2": ["p16.3", "p16.1", "p12", "q22.1", "q22.3", "q34", "q36.3"],
    "3": ["p24.3", "q24", "q26.1"],
    "4": ["p15.1", "q13.1", "q28.3", "q32.1", "q32.3", "q34.3"],
    "5": ["p14.3", "p14.1", "q14.3", "q21.1", "q21.3", "q23.1", "q23.3", "q34"],
    "6": ["p12.3", "p12.1", "q12", "q16.1", "q16.3", "q22.31"],
    "7": ["p21.3", "p21.1", "q21.11"],
    "8": ["p22", "q21.11", "q21.3", "q23.3"],
    "9": ["p21.3", "p21.1", "q31.1"],
    "10": ["q21.1", "q21.3", "q23.1", "q25.1"],
    "11": ["p14.3", "p12", "q14.1", "q14.3", "q22.1", "q22.3"],
    "12": ["p12.3", "p12.1", "q12", "q21.31", "q21.33"],
    "13": ["q21.1", "q21.33", "q31.1", "q31.3", "q33.1", "q33.3"],
    "14": ["q12", "q21.1", "q21.3", "q31.1", "q31.3"],
    "15": ["q14", "q21.1", "q21.3"],
    "16": ["q21", "q23.1"],
    "17": ["p12", "q22", "q23.2", "q24.3"],
    "18": ["q12.1", "q22.1"],
    "19": ["p13.2", "p13.12", "q13.12", "q13.2", "q13.32", "q13.41", "q13.43"],
    "20": ["p12.3", "p12.1", "q12", "q13.2"],
    "21": ["q21.1", "q21.3"],
    "22": ["q12.1", "q12.3", "q13.2", "q13.32"],
    "X": ["p21.3", "p21.1", "q21.1", "q21.31", "q25", "q27.3"],
    "Y": ["p11.31", "q11.221", "q11.223"]
}

HETERO_EU_DICT = {
    "1": ["p31.1", "p34.1", "p13.1", "q25.3", "q32.1"],
    "2": ["p16.1", "p11.2", "p21", "q31.1", "q33.1"],
    "3": ["p24.3", "p21.31", "p21.1", "q21.1", "q23"],
    "4": ["p15.1", "p14", "p16.1", "q21.3", "q25"],
    "5": ["p14.3", "p13.3", "p15.1", "q23.2", "q31.3"],
    "6": ["p12.3", "p21.31", "p21.1", "q21", "q23.1"],
    "7": ["p21.3", "p15.3", "p13", "q22.1", "q36.1"],
    "8": ["p22", "p21.1", "p23.1", "q11.23", "q22.1"],
    "9": ["p21.3", "p22.3", "p13.3", "q22.1", "q32"],
    "10": ["q21.1", "p13", "p11.23", "q22.1", "q26.13"],
    "11": ["p14.3", "p15.1", "p13", "q13.1", "q23.3"],
    "12": ["p12.3", "p13.31", "p11.23", "q13.13", "q24.31"],
    "13": ["q21.1", "q14.11", "q14.3", "q22.1", "q32.3"],
    "14": ["q12", "q11.2", "q13.3", "q24.1", "q24.3"],
    "15": ["q14", "q11.2", "q13.3", "q15.3", "q22.31"],
    "16": ["q21", "p13.13", "p11.2", "q12.1", "q22.1"],
    "17": ["p12", "p13.1", "p11.2", "q11.2", "q21.31"],
    "18": ["q12.1", "p11.23", "p11.21", "q21.1", "q12.2"],
    "19": ["p13.2", "p13.13", "p13.11", "q13.13", "q13.33"],
    "20": ["p12.3", "p11.23", "p11.21", "q13.11", "q13.13"],
    "21": ["q21.1", "q11.2", "q21.2", "q22.11", "q22.13"],
    "22": ["q12.1", "q11.21", "q12.2", "q13.1", "q13.31"],
    "X": ["p21.3", "p22.13", "p11.4", "q13.1", "q22.1"],
    "Y": ["p11.31", "p11.32", "p11.2", "q11.21", "q11.222"]
}

TANDEM_REPEAT_DICT = {
    "Y": ["q12"],
    "13": ["p13", "p12", "p11.2"],
    "14": ["p13", "p12", "p11.2"],
    "15": ["p13", "p12", "p11.2"],
    "21": ["p13", "p12", "p11.2"],
    "22": ["p13", "p12", "p11.2"]
}

DISTANCE_METRICS_LIST = ["Normalized Euclidean", "Cosine", "Manhattan",
                         "Descriptor", "DSSIM", "LPIPS",
                         "K-S", "Wasserstein"]

SCIENTIFIC_NAMES = {
    "Human": 'H. sapiens',
    "Chimp": 'P. troglodytes',
    "Mouse": 'M. musculus',
    "Drosophila melanogaster": 'D. melanogaster',
    "Saccharomyces cerevisiae": 'S. cerevisiae',
    "Arabidopsis thailana": 'A. thaliana',
    "Paramecium caudatum": 'P. caudatum',
    "Pyrococcus furiosus": 'P. furiosus',
    "Escherichia coli": 'E. coli',
    "Maize": 'Z. mays',
    "Dictyostelium discoideum": 'D. discoideum',
    "Aspergillus terreus": 'A. terreus',
    "Aspergillus nidulans": 'A. nidulans'
}

GENOME_LENGTH = {"Human": 3117,
                 "Chimp": 3178,
                 "Mouse": 2723,
                 "Drosophila melanogaster": 80,
                 "Saccharomyces cerevisiae": 12,
                 "Arabidopsis thailana": 119,
                 "Paramecium caudatum": 30,
                 "Pyrococcus furiosus": 2,
                 "Escherichia coli": 5,
                 "Maize": 2179,
                 "Dictyostelium discoideum": 34,
                 "Aspergillus terreus": 29,
                 "Aspergillus nidulans": 30}
REVERSE_SCIENTIFIC_NAMES = {v: k for k, v in SCIENTIFIC_NAMES.items()}

RESOLUTION_DICT = {2: 2, 3: 2, 4: 2, 5: 2, 6: 2, 7: 2, 8: 4, 9: 4, 10: 4, 11: 4, 12: 4}
BITS_DICT = {"Human": 8, "Chimp": 8, "Mouse": 8, "Drosophila melanogaster": 8, "Saccharomyces cerevisiae": 6,
             "Arabidopsis thailana": 8, "Paramecium caudatum": 8, "Pyrococcus furiosus": 6,
             "Escherichia coli": 5, "Maize": 6, "Dictyostelium discoideum": 9, "Aspergillus terreus": 8,
             "Aspergillus nidulans": 8, "Custom": 8}
