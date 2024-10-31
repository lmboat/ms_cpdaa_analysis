# MS CpDAA: Mass Spectrometry-based Chemoproteomics Detected Amino Acids Analysis Suite

The MS CpDAA Analysis Suite was created and developed by Lisa Boatner. 

About: Processing pipeline for aggregating and analyzing outputs from various mass-spectromtery experiments (ex. ABPP, Silac, TMT, Iodination, Label-Free, N-terminal, DIA), software packages (ex. FragPipe, MaxQuant, Skyline, ProLuCID), and  replicates. 

## Output
1. Sites of labeling mapped to amino acid positions in a reference protein sequence
2. Sites of labeling for every labeled peptides, residues and proteins in each experiment/replicate
3. Counts of distinct peptides, residues and proteins in each experiment/replicate
4. Matrices of distinct and shared peptides, residues and proteins across different experiments/replicates
5. Medians of raw ratios, log2 ratios, or ion intensities for every peptide, residue and protein in each experiment and replicate
6. Average of medians for raw ratios, log2 ratios, or ion intensities for every peptide, residue and protein in each experiment and replicate

## Set Up

1. Download python
2. Download latest analysis script
3. Create a folder called "data"
4. Move desired output folders into the "data" folder

- data
  - Experiment-1
    - Replicate-1
      - psm.tsv
      - peptide.tsv
      - peptide_label_quant.tsv
      - protein.fas 
    - Replicate-2   
  - Experiment-2
    - Replicate-1
    - Replicate-2
    - Replicate-3
  ... 
  - Experiment-50  

## Usage

```python
python3 231101_ms_cpdaa_process.py -h
```

## Test Usage
```python
python3 231101_ms_cpdaa_process.py -exp 'isotop' -lpm '521.3074' -hpm '527.3213' -dbv '20'
```

## Identification Experiments Usage

### Without a Probe Mass
```python
python3 231101_ms_cpdaa_process.py -pm '0'
```

### With a Probe Mass
```python
python3 231101_ms_cpdaa_process.py -pm 'your_probe_mass'
```

### Label Free with a Probe Mass
```python
python3 231101_ms_cpdaa_process.py -exp 'lfq' -pm 'your_probe_mass'
```

### N-terminal with a Probe Mass
```python
python3 231101_ms_cpdaa_process.py -aa 'N-term' -pm 'your_probe_mass'
```

## Quantitative Experiments Usage

### competitive-ABPP, isoTOP-ABPP, silac, iodiniation
```python
python3 231101_ms_cpdaa_process.py -aa 'C;K;H' -exp 'your_quant_experiment_type' -lpm 'your_light_probe_mass' -hpm 'your_heavy_probe_mass' 
```

## Additional Notes
* -aa : add any or all amino acids (ex. 'A;C;D;E;F;G;H;I:K;L;M;N;P;Q;R;S;T;V;W;Y')
* -rr : ratios listed are not Log2Ratios (default False)
* -mid : separate multiplexed peptides into separate enteries from ProteinID_AA#1_AA#2 -> ProteinID_AA#1 and ProteinID_AA#2
* -pepc : compare the intensities of modified and unmodified peptides (default False)

<p align="center">
  <img src="https://github.com/lmboat/cpdaadb/assets/35751646/68c3c416-b213-4a51-82c7-317a0df17af6">
</p>

## Works published using this suite:
1. Cao J, Boatner LM, Desai HS, Burton NR, Armenta E, Chan NJ, Castellón JO, Backus KM. Multiplexed CuAAC Suzuki-Miyaura Labeling for Tandem Activity-Based Chemoproteomic Profiling. Anal Chem. 2021 Feb 2;93(4):2610-2618. doi: 10.1021/acs.analchem.0c04726. Epub 2021 Jan 20. PMID: 33470097; PMCID: PMC8849040.
2. Yan T, Desai HS, Boatner LM, Yen SL, Cao J, Palafox MF, Jami-Alahmadi Y, Backus KM. SP3-FAIMS Chemoproteomics for High-Coverage Profiling of the Human Cysteinome*. Chembiochem. 2021 May 14;22(10):1841-1851. doi: 10.1002/cbic.202000870. Epub 2021 Feb 18. PMID: 33442901; PMCID: PMC8942723.
3. Tang KC, Cao J, Boatner LM, Li L, Farhi J, Houk KN, Spangle J, Backus KM, Raj M. Tunable Amine-Reactive Electrophiles for Selective Profiling of Lysine. Angew Chem Int Ed Engl. 2022 Jan 26;61(5):e202112107. doi: 10.1002/anie.202112107. Epub 2021 Dec 16. PMID: 34762358; PMCID: PMC10111338.
4. Yan T, Julio AR, Villanueva M, Jones AE, Ball AB, Boatner LM, Turmon AC, Nguyễn KB, Yen SL, Desai HS, Divakaruni AS, Backus KM. Proximity-labeling chemoproteomics defines the subcellular cysteinome and inflammation-responsive mitochondrial redoxome. Cell Chem Biol. 2023 Jul 20;30(7):811-827.e7. doi: 10.1016/j.chembiol.2023.06.008. Epub 2023 Jul 6. PMID: 37419112; PMCID: PMC10510412.
5. Castellón JO, Ofori S, Burton NR, Julio AR, Turmon AC, Armenta E, Sandoval C, Boatner LM, Takayoshi EE, Faragalla M, Taylor C, Zhou AL, Tran K, Shek J, Yan T, Desai HS, Fregoso OI, Damoiseaux R, Backus KM. Chemoproteomics Identifies State-Dependent and Proteoform-Selective Caspase-2 Inhibitors. J Am Chem Soc. 2024 Jun 5;146(22):14972-14988. doi: 10.1021/jacs.3c12240. Epub 2024 May 24. PMID: 38787738.
6. Ofori S, Desai HS, Shikwana F, Boatner LM, Dominguez Iii ER, Castellón JO, Backus KM. Generating cysteine-trypsin cleavage sites with 2-chloroacetamidine capping. Chem Commun (Camb). 2024 Aug 15;60(67):8856-8859. doi: 10.1039/d4cc01583e. PMID: 39081146.
7. Takechi S, Ngo C, Burton N, Villanueva M, Boatner L, Yu F, et al. Silyl Ether Enables High Coverage Chemoproteomic Interaction Site Mapping. ChemRxiv. 2024; doi:10.26434/chemrxiv-2024-21r7b This content is a preprint and has not been peer-reviewed.

## Contact
Lisa Boatner - lisamboatner@gmail.com
