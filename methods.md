## 06. How does mamp-ml compare to other structural models

Assessing how mamp-ml compares to other structural models is important for determining its potential advantages. There are many structural models that exist but we will test a couple to A) see how they perform and B) determine if the outputs of those models in combination with a simple regession model can outcompete mamp-ml performance. 

```
# we then wrote a script that inputs All_LRR_PRR_ligand_data.xlsx and outputs an updated excel sheet that has co-receptor information as well (we will want to compare structural models with and without co-receptor information)

# install packages and run script
conda install -y -c bioconda -c conda-forge parasail-python
python 00_scripts/01_serk_pipeline.py
```

Next, we need to build the multiple sequence alignments (MSAs) for each receptor sequence. This is an essential step for structural modeling, as high-quality MSAs will be used as input features for downstream structure prediction tools such as Boltz and will directly impact the accuracy of the resulting structural models.

```
# output for msas
mkdir -p 04_receptor_only_msas

# install colabfold
pip install uv
uv venv ~/envs/colabfold --python 3.11
source ~/envs/colabfold/bin/activate
uv pip install colabfold

# MMseqs2 static binary for macOS (Homebrew's version can lag behind what ColabFold needs)
colabfold_search 03_out_data/lrr_domain_sequences.fasta ../colabfold_db 04_receptor_only_msas \
  --mmseqs /opt/homebrew/bin/mmseqs --use-env 1 --use-templates 0 --threads 14
```


```
# run in root of mamp-ml-structure -- will generate serk_ectodomain_sequences.fasta in 03_out_data
python - <<'EOF'
import pandas as pd, parasail
df = pd.read_csv("SERK_coreceptor.tsv", sep="\t", usecols=["SERK_ID", "SERK_Sequence"])
df = df[df.SERK_ID.ne("NA") & df.SERK_Sequence.notna()].drop_duplicates("SERK_ID")
ref = df.loc[df.SERK_ID == "AT4G33430", "SERK_Sequence"].iat[0]
# AtBAK1 26-220: mature ectodomain (after signal peptide, before transmembrane helix)
START, END = 26, 220
with open("03_out_data/serk_ectodomain_sequences.fasta", "w") as fh:
    for gid, seq in zip(df.SERK_ID, df.SERK_Sequence):
        tr = parasail.nw_trace_scan_16(seq, ref, 11, 1, parasail.blosum62).traceback
        qi = ri = 0
        s = e = None
        for q, r in zip(tr.query, tr.ref):
            qi += q != "-"
            ri += r != "-"
            if q != "-" and r != "-":
                if s is None and ri >= START:
                    s = qi
                if ri <= END:
                    e = qi
        ecto = seq[s - 1:e]
        fh.write(f">{gid}\n{ecto}\n")
        print(f"{gid}\t{len(seq)} -> {len(ecto)} aa")
EOF


# output for msas - co-receptor
source ~/envs/colabfold/bin/activate
mkdir -p 05_coreceptor_msas

# search to build msas for co-recptor
colabfold_search 03_out_data/serk_ectodomain_sequences.fasta ../colabfold_db 05_coreceptor_msas \
  --mmseqs /opt/homebrew/bin/mmseqs --use-env 1 --use-templates 0 --threads 14

```

Next, we need to prepare the necessary input files in the specific formats required by each structural model we plan to use. This involves running scripts that convert our processed sequence and annotation data into model-specific input directories and file types to ensure compatibility for downstream structure prediction analyses.

```
  python 00_scripts/02_build_inputs.py --primary_only
  python 00_scripts/02_build_inputs.py --complex

# run the build input to run throguh AF/Boltz/RF3 
python 00_scripts/02_build_inputs.py \
  SERK_coreceptor.tsv \
  03_out_data/lrr_domain_sequences.fasta \
  04_receptor_only_msas \
  06_input_files_fold_primary_only

# outputs the following:
├── 06_input_files_fold_primary_only/
│   └── boltz/
│   └── rf3.json
│   └── rf3.json

```
First, we will install and run all receptor-epitope and receptor-epitope-co-receptor members from Model_Building section though Boltz-1. We will begin by running Boltz, a diffusion-based structural modeling tool, on our prepared input files. This will allow us to generate structural predictions for all receptor-epitope and receptor-epitope-co-receptor combinations, serving as a benchmark for our model comparison.

```
# --- Boltz-1 + Boltz-2 (same package, pick model with --model) -------------
# boltz-community is a maintained fork with native `--accelerator mps`.
uv venv envs/boltz --python 3.11
source envs/boltz/bin/activate
uv pip install boltz-community pandas pyyaml
#export PYTORCH_ENABLE_MPS_FALLBACK=1

# ------ run on A5000 card on-----

# make directory for outputs
mkdir -p 08_output_fold_primary_only

# --------------- Boltz-1 (primary only) ---------------
sbatch \
  --job-name="mamp_prediciton_ml2" \
  --partition=savio3_gpu \
  --qos=a40_gpu3_normal \
  --account=fc_kvkallow \
  --time=16:00:00 \
  --gres=gpu:A40:1 \
  --cpus-per-task=8 \
  --exclude=n0215.savio3 \
  --mail-user=dmstev@berkeley.edu \
  --mail-type=ALL \
  --wrap 'cd /global/scratch/projects/vector_kvklab/dani/mamp-ml-structure && source /global/scratch/users/dmstev/mamp_prediction_ml/envs/boltz/bin/activate && nvidia-smi && python -c "import torch; assert torch.cuda.is_available(), \"no CUDA\"" && boltz predict 06_input_files_fold_primary_only/boltz --model boltz1 --out_dir 08_output_fold_primary_only/boltz1 --accelerator gpu --devices 1 --recycling_steps 3 --sampling_steps 150 --max_parallel_samples 1 --diffusion_samples 1 --output_format mmcif'


# --------------- Boltz-2 (primary only) ---------------
sbatch \
  --job-name="mamp_prediciton_ml2_boltz2" \
  --partition=savio3_gpu \
  --qos=a40_gpu3_normal \
  --account=fc_kvkallow \
  --time=16:00:00 \
  --gres=gpu:A40:1 \
  --cpus-per-task=8 \
  --exclude=n0215.savio3 \
  --mail-user=dmstev@berkeley.edu \
  --mail-type=ALL \
   --wrap 'cd /global/scratch/projects/vector_kvklab/dani/mamp-ml-structure && source /global/scratch/users/dmstev/mamp_prediction_ml/envs/boltz/bin/activate && nvidia-smi && python -c "import torch; assert torch.cuda.is_available(), \"no CUDA\"" && boltz predict 06_input_files_fold_primary_only/boltz --model boltz2 --out_dir 08_output_fold_primary_only/boltz2 --accelerator gpu --devices 1 --recycling_steps 3 --sampling_steps 150 --max_parallel_samples 1 --diffusion_samples 1 --output_format mmcif'


  
# -------- run on A40 card -------------
uv venv envs/rf3 --python 3.12
source envs/rf3/bin/activate

# assume torch with CUDA support is already installed
uv pip install pandas
uv pip install "rc-foundry[all] @ git+https://github.com/fnachon/foundry.git"

# downloads latest RF3 checkpoint to ~/.foundry/checkpoints
foundry install rf3   


# --------------- RF3 (primary only) ---------------
sbatch \
  --job-name="rf3_mamp_prediction_ml" \
  --partition=savio3_gpu \
  --qos=a40_gpu3_normal \
  --account=fc_kvkallow \
  --time=20:00:00 \
  --gres=gpu:A40:1 \
  --cpus-per-task=8 \
  --exclude=n0215.savio3 \
  --mail-user=dmstev@berkeley.edu \
  --mail-type=ALL \
  --wrap 'cd /global/scratch/projects/vector_kvklab/dani/mamp-ml-structure && source /global/scratch/users/dmstev/mamp_prediction_ml/envs/rf3/bin/activate && nvidia-smi && python -c "import torch; assert torch.cuda.is_available(), \"no CUDA\"" && rf3 fold inputs="06_input_files_fold_primary_only/rf3.json" out_dir="08_output_fold_primary_only/rf3" early_stopping_plddt_threshold=0.0 diffusion_batch_size=3 raise_if_missing_msa_for_protein_of_length_n=50 skip_existing=true'

```
RF3’s skip_existing=true did not skip finished folds because our installed Foundry build checks for *_metrics.csv, while completed runs only write *_ranking_scores.csv. As a workaround, copy each *_ranking_scores.csv to a matching *_metrics.csv (or patch rf3/utils/inference.py to look for ranking_scores) before resubmitting so existing outputs are recognized and not overwritten.

```

# note, run the below comand inbetween runs:
# current version 
cd /global/scratch/projects/vector_kvklab/dani/mamp-ml-structure
find 08_output_fold_primary_only/rf3 -name '*_ranking_scores.csv' | while read f; do
  m="${f%_ranking_scores.csv}_metrics.csv"
  [ -e "$m" ] || cp "$f" "$m"
done

```

Finally, we will run the same files with alphafold2 and alphafold3.

```
# --------------- Colabfold - AlphaFold2-Multimer (primary only) ---------------

sbatch \
  --job-name="colabfold_af2_primary" \
  --partition=savio3_gpu \
  --qos=a40_gpu3_normal \
  --account=fc_kvkallow \
  --time=20:00:00 \
  --gres=gpu:A40:1 \
  --cpus-per-task=8 \
  --exclude=n0215.savio3 \
  --mail-user=dmstev@berkeley.edu \
  --mail-type=ALL \
  --wrap 'cd /global/scratch/projects/vector_kvklab/dani/mamp-ml-structure && module load anaconda3 && source activate localfold && module load gcc/10.5.0 && export PATH="/global/scratch/users/dmstev/localcolabfold/colabfold-conda/bin:$PATH" && nvidia-smi && which colabfold_batch && mkdir -p 08_output_fold_primary_only/colabfold_af2 && colabfold_batch --model-type alphafold2_multimer_v3 --num-recycle 3 --num-models 3 --pair-mode unpaired 06_input_files_fold_primary_only/colabfold_msa 08_output_fold_primary_only/colabfold_af2'


# --------------- AlphaFold3 (primary only) ---------------


```

Now that we've done all the primary receptor-epitope combinations, we need to do the complex as well. Therefore, similar as above, we will run the following commands for the complexes.

```
# make a directory that holds the outputs
mkdir -p 09_output_fold_complex

# --------------- Boltz-1 (receptor complex (primary-peptide-co-receptor)) ---------------
sbatch \
  --job-name="mamp_boltz1_complex" \
  --partition=savio3_gpu \
  --qos=a40_gpu3_normal \
  --account=fc_kvkallow \
  --time=16:00:00 \
  --gres=gpu:A40:1 \
  --cpus-per-task=8 \
  --exclude=n0215.savio3 \
  --mail-user=dmstev@berkeley.edu \
  --mail-type=ALL \
  --wrap 'cd /global/scratch/projects/vector_kvklab/dani/mamp-ml-structure && source /global/scratch/users/dmstev/mamp_prediction_ml/envs/boltz/bin/activate && nvidia-smi && python -c "import torch; assert torch.cuda.is_available(), \"no CUDA\"" && boltz predict 07_input_files_complex/boltz --model boltz1 --out_dir 09_output_fold_complex/boltz1 --accelerator gpu --devices 1 --recycling_steps 3 --sampling_steps 150 --max_parallel_samples 1 --diffusion_samples 1 --output_format mmcif'

# --------------- Boltz-2 (receptor complex (primary-peptide-co-receptor)) ---------------
sbatch \
  --job-name="mamp_boltz1_complex" \
  --partition=savio3_gpu \
  --qos=a40_gpu3_normal \
  --account=fc_kvkallow \
  --time=16:00:00 \
  --gres=gpu:A40:1 \
  --cpus-per-task=8 \
  --exclude=n0215.savio3 \
  --mail-user=dmstev@berkeley.edu \
  --mail-type=ALL \
  --wrap 'cd /global/scratch/projects/vector_kvklab/dani/mamp-ml-structure && source /global/scratch/users/dmstev/mamp_prediction_ml/envs/boltz/bin/activate && nvidia-smi && python -c "import torch; assert torch.cuda.is_available(), \"no CUDA\"" && boltz predict 07_input_files_complex/boltz --model boltz2 --out_dir 09_output_fold_complex/boltz2 --accelerator gpu --devices 1 --recycling_steps 3 --sampling_steps 150 --max_parallel_samples 1 --diffusion_samples 1 --output_format mmcif'

# --------------- RF3 (receptor complex (primary-peptide-co-receptor)) ---------------
sbatch \
  --job-name="rf3_mamp_complex" \
  --partition=savio3_gpu \
  --qos=a40_gpu3_normal \
  --account=fc_kvkallow \
  --time=20:00:00 \
  --gres=gpu:A40:1 \
  --cpus-per-task=8 \
  --exclude=n0215.savio3 \
  --mail-user=dmstev@berkeley.edu \
  --mail-type=ALL \
  --wrap 'cd /global/scratch/projects/vector_kvklab/dani/mamp-ml-structure && source /global/scratch/users/dmstev/mamp_prediction_ml/envs/rf3/bin/activate && nvidia-smi && python -c "import torch; assert torch.cuda.is_available(), \"no CUDA\"" && rf3 fold inputs="07_input_files_complex/rf3.json" out_dir="09_output_fold_complex/rf3" early_stopping_plddt_threshold=0.0 diffusion_batch_size=3 raise_if_missing_msa_for_protein_of_length_n=50 skip_existing=true'

find 09_output_fold_complex/rf3 -name '*_ranking_scores.csv' | while read f; do
  m="${f%_ranking_scores.csv}_metrics.csv"
  [ -e "$m" ] || cp "$f" "$m"
done

# --------------- Colabfold - AlphaFold2-Multimer (receptor complex (primary-peptide-co-receptor)) ---------------

# note need to run on H200 card....
sbatch \
  --job-name="colabfold_af2_complex" \
  --partition=savio3_gpu \
  --qos=a40_gpu3_normal \
  --account=fc_kvkallow \
  --time=20:00:00 \
  --gres=gpu:A40:1 \
  --cpus-per-task=8 \
  --exclude=n0215.savio3 \
  --mail-user=dmstev@berkeley.edu \
  --mail-type=ALL \
  --wrap 'cd /global/scratch/projects/vector_kvklab/dani/mamp-ml-structure && module load anaconda3 && source activate localfold && module load gcc/10.5.0 && export PATH="/global/scratch/users/dmstev/localcolabfold/colabfold-conda/bin:$PATH" && nvidia-smi && which colabfold_batch && mkdir -p 09_output_fold_complex/colabfold_af2 && colabfold_batch --model-type alphafold2_multimer_v3 --num-recycle 3 --num-models 3 --pair-mode unpaired 07_input_files_complex/colabfold_msa 09_output_fold_complex/colabfold_af2'

  sbatch \
  --job-name="colabfold_af2_complex" \
  --partition=savio4_gpu \
  --qos=savio_lowprio \
  --account=fc_kvkallow \
  --time=20:00:00 \
  --gres=gpu:H200:1 \
  --cpus-per-task=8 \
  --mail-user=dmstev@berkeley.edu \
  --mail-type=ALL \
  --wrap 'cd /global/scratch/projects/vector_kvklab/dani/mamp-ml-structure && module load anaconda3 && source activate localfold && module load gcc/10.5.0 && export PATH="/global/scratch/users/dmstev/localcolabfold/colabfold-conda/bin:$PATH" && nvidia-smi && which colabfold_batch && mkdir -p 09_output_fold_complex/colabfold_af2 && colabfold_batch --model-type alphafold2_multimer_v3 --num-recycle 3 --num-models 3 --pair-mode unpaired 07_input_files_complex/colabfold_msa 09_output_fold_complex/colabfold_af2'


# --------------- AlphaFold3 (receptor complex (primary-peptide-co-receptor)) ---------------

```

We predicted structures for each receptor–ligand pair (and, separately, receptor–ligand–SERK co-receptor complexes) with Boltz-1, Boltz-2, RoseTTAFold3, and AlphaFold2-Multimer (ColabFold), using precomputed monomer MSAs for the receptor and co-receptor chains. For each model we treated the reported interface predicted TM score (ipTM) as a confidence proxy for whether the predicted complex supports a productive protein–protein interaction, with higher ipTM indicating stronger predicted interface confidence. We will not plot recognition (immunogenic) vs. non-recognition (non-immunogenic) and save the plors in 10_model_plots.

```
# once the combinations have been completed
Rscript 00_scripts/03_structure_models_comparisons.R
```
AlphaFold’s ipTM is a whole-chain interface score, so large length imbalances, disordered regions, or non-contacting accessory domains can systematically inflate or depress the metric even when the predicted contacting domains are unchanged (Dunbrack et al.). This is especially problematic for domain–peptide complexes such as LRR receptors with short ligands, where chain lengths differ by orders of magnitude. We therefore also scored interfaces with ipSAE, which restricts the calculation to interchain residue pairs with good predicted aligned error (PAE) and rescales the TM length parameter to those well-aligned residues, providing a more reliable proxy for interaction confidence than ipTM alone.

```
# Interface scoring with DunbrackLab IPSAE (https://github.com/DunbrackLab/IPSAE)
# on top-ranked primary-only folds. Model-specific PAE / distance cutoffs:
#   AlphaFold2 (colabfold_af2): 15 Å / 15 Å
#   Boltz-1, Boltz-2, RF3:      10 Å / 10 Å
# Scores (Type == max ipSAE) and plots go under 11_ipsae_plots/; existing
# per-pair score files are skipped. RF3 is run through IPSAE's AF3 (json+cif)
# path with a patched summary_confidences.json that supplies chain_pair_iptm
# from the overall iptm (needed by the AF3 parser; plotted metric is ipSAE).

# this will calculate ipSAE scores (takes about 1-2 hours) and plot in 11_ipsae_plots
Rscript 00_scripts/04_model_comparisons_ipsae.R
```
Finally, we treated each structure model’s interface confidence (ipTM, and separately ipSAE) as a binary immunogenicity classifier by applying a fixed score threshold and comparing predictions to experimental labels on the shared mamp-ml test set. We summarized performance with accuracy and macro- and weighted-F1 scores and ranked the best-performing structural models against SeqOnly and mamp-ml baselines from their published classification reports.

```
# compares F1 metrics of structural models vs. pLMs including mamp-ml
Rscript 00_scripts/05_structure_model_accuracy.R \
  --seqonly ../model_results/02_immuno_stratify_esm2_with_receptor/classification_report.txt \
  --mamp-ml ../model_results/06_immuno_stratify_esm2_bfactor_weighted/classification_report.txt
```

