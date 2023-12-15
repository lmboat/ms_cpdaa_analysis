# MS CpDAA: Mass Spectrometry-based Chemoproteomics Detected Amino Acids Analysis Suite

The MS CpDAA Analysis Suite was created and developed by Lisa Boatner. 

About: Processing pipeline for aggregating and analyzing outputs from various mass-spectromtery experiments (ex. ABPP, Silac, TMT, Iodination, Label-Free, N-terminal), software packages (ex. FragPipe, MaxQuant, Skyline, ProLuCID), and  replicates. 

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

### Identification Experiments without a Probe Mass
```python
python3 231101_ms_cpdaa_process.py -pm '0'
```

### Identification Experiments with a Probe Mass
```python
python3 231101_ms_cpdaa_process.py -pm 'your_probe_mass'
```

### Label Free Identification Experiments with a Probe Mass
```python
python3 231101_ms_cpdaa_process.py -exp 'lfq' -pm 'your_probe_mass'
```

### N-terminal Identification Experiments with a Probe Mass
```python
python3 231101_ms_cpdaa_process.py -aa 'N-term' -pm 'your_probe_mass'
```

### Quantitative Experiments (ex. isotop, silac, iodiniation) 
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

## Contact
Lisa Boatner - lisaboatner@g.ucla.edu
