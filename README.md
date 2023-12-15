# MS CpDAA Analysis: Mass Spectrometry-based Chemoproteomics Detected Amino Acids Analysis Suite

The MS CpDAA Analysis Suite was created and developed by Lisa Boatner. 

About: Processing pipeline for aggregating and analyzing outputs from various mass-spectromtery experiments (ex. ABPP, Silac, TMT, Iodination, Label-Free, N-terminal), software packages (ex. FragPipe, MaxQuant, Skyline, ProLuCID), and  replicates.

## Set Up

1. Download python
2. Download latest scripts
3. Create a "data" folder
4. Move any output folders into the data folder

<p align="center">
  <img width="460" height="250" src="https://github.com/lmboat/cpdaadb/assets/35751646/7c1c9c21-31b8-481c-a448-5651c6c9030c">
</p>


## Usage

```python
python3 231101_ms_cpdaa_process.py -h
```
<p align="center">
  <img width="460" height="250" src="https://github.com/lmboat/cpdaadb/assets/35751646/68c3c416-b213-4a51-82c7-317a0df17af6">
</p>

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

## Contact
Lisa Boatner - lisaboatner@g.ucla.edu
