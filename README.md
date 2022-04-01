# Degradation-aware isoform quantification

daiso is a tool for isoform abundance estimation from ONT direct RNA-seq that models degradation bias. 

## Installation

Clone the package:
```
git clone https://github.com/jleechung/daiso
cd daiso
```

Set up a python virtual environment:
```
python3 -m venv myenv
source ./myenv/bin/activate
```

Install requirements:
```
pip3 install -r requirements.txt
```

Test the installation:
```
python3 daiso.py -h
```
This should print the help message.

## Basic usage

daiso requires only one input: alignment (bam files) of long ONT reads to the transcriptome with [minimap2](https://github.com/lh3/minimap2). 

From raw ONT reads (fasta or fastq files), you can execute
```
# reads.fa are your reads, ref.fa is a reference transcriptome
minimap2 -ax map-ont -N10 ref.fa reads.fa | samtools view -b > alnm.bam
samtools sort -O BAM alnm.bam > alnm_sorted.bam
samtools index alnm_sorted.bam
```

Now, run daiso:
```
python3 daiso.py -a alnm_sorted.bam
```

## Examples

We provide example data in `examples/` which comprises simulated degraded reads aligning to the [human transcriptome](https://www.gencodegenes.org/human/) (see transcript sequences). The expected output is included in the `output/`. To generate the output, execute
```
python3 daiso.py -a examples/subset.bam -o output/out
```

## Advanced usage

We provide customizable options for alignment filtering, degradation estimation and parameter inference. In addition, if one is only keen on estimating degradation, daiso can be run with the `--no_quant` flag. 

```
usage: daiso -a alignment.bam [options]

Degradation-aware isoform quantification

Required arguments:
  -a ALIGNMENT, --alignment ALIGNMENT
                        BAM file (sorted and indexed) containing reads aligned
                        to the transcriptome

Optional arguments:
  -o OUTPUT, --output OUTPUT
                        Output prefix. The count and degradation plot files
                        will be written out to <output>_counts.tsv and
                        <output>_survival.png respectively (default: out)
  --no_quant            Estimate degradation only. No quantification is run.
                        (default: False)
  --seed SEED           Set seed (default: 42)

Parameters for alignment filtering:
  --filter_distance FILTER_DISTANCE
                        Filter alignments where end position is further than
                        this distance away from the annotated 3 prime end
                        (default: 50 bp)
  --filter_score FILTER_SCORE
                        Filter alignments where the alignment score is lower
                        than this fraction of the best alignment score for
                        this read (default: 0.95)

Parameters for degradation estimation:
  --deg_rate DEG_RATE   Degradation rate in (0,1). If this argument is
                        supplied, degradation rate will not be estimated from
                        the data and deg_const will be True
  --deg_const           Uses exact constant degradation model for read length-
                        isoform agreement (default:False)
  --bin_size BIN_SIZE   Bin isoform lengths into bins of this size (default:
                        500 bp)
  --min_read_count MIN_READ_COUNT
                        Only estimate degradation with isoforms with at least
                        this many reads (default: 5)
  --min_iso_count MIN_ISO_COUNT
                        Only estimate degradation with length bins supported
                        by this many isoforms (default: 1)
  --full_len_tol FULL_LEN_TOL
                        Consider a read full length if it is within this many
                        bp of the annotated isoform length (default: 50 bp)
  --delta DELTA         Left and right shift for estimating exact
                        probabilities for ecdf (default: 50 bp)
  --return_survival     Return the survival function (1-ecdf). Writes values
                        to <output>_survival.tsv (default: True)

Parameters for inference:
  --inference {EM,VB}   Inference algorithm - expectation maximization or
                        variational Bayesian inference (default: EM)
  --max_iter MAX_ITER   Maximum number of iterations (default: 200)
  --return_loglik       Return the log likelihood over iterations for EM.
                        Writes values to <output>_loglik.tsv (default: True)
  --prior {symmetric,gamma_hyper}
                        Choice of prior - symmetric prior or a gamma
                        hyperprior (default: gamma_hyper)
  --alpha_zero ALPHA_ZERO
                        Value of the concentration parameter for the symmetric
                        prior (default: 1)
  --gamma_rate GAMMA_RATE
                        Rate parameter for the gamma hyperprior (default: 5)
  --gamma_scale GAMMA_SCALE
                        Scale parameter for the gamma hyperprior (default: 5)
  --return_cred         Return credible interals for parameters (default:
                        False)
  --cred_int CRED_INT   Width of credible interval (default: 0.95)

General help:
  -h, --help            Show this help message and exit
  -v, --version         Print version
```

## Software limitations

Not supported on windows due to dependencies on pysam. Consider WSL.


