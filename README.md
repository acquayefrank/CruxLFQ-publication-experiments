# CruxLFQ-publication-experiments

In order to reporduce the experiments in "Label-free quantification in the Crux toolkit" follow the steps below:

## Assumptions

- You have basic familiarity with the command line.
- You have access to a computer running Windows, macOS, or Linux. (Some experiments were run on Linux and others on Windows)
- You have sufficient permissions to install and run software.
- Required dependencies (such as Python, if needed) are installed.
- Code is found in `src`, hence each time an instruction is provided to run it is assumed the user is running from `src`.
- You can change scripts accordingly for the specific platform (OS) you are using.
- Scripts will be run in order in which it is presented here - from top to bottom.
- When using crux users would have a seemless pipeline, howver in order to compare tools that are not interoperable a bit more work needs to be done hence this many scripts.
- It is assumed that results folder will contain a single folder(instance) of any results type e.g. only one instance of sage search results.
- Note that all plots will be outputted to the src directory or the directory from which you run the code

# Database Searching:

## Tide:

As at the time of this writing, crux-lfq had not been merged crux main.
The PR can be viewed [here](https://github.com/crux-toolkit/crux-toolkit/pull/717) and the raw code [here](https://github.com/acquayefrank/crux-toolkit/tree/feat/LFQ). In order to run `lfq` in crux, download a custom build of crux: [crux-lfq]()

- Once done downloading, place the binary in `bin` folder

- Download the fasta file from [ProteomeXchange - PXD003881](https://ftp.pride.ebi.ac.uk/pride/data/archive/2018/05/PXD003881/Human_ecoli_trypsin_1501v_uniprot_sprot.fasta) and place the downloaded file in the folder `data`

- Download the spectrum files from [ProteomeXchange - PXD005590](https://www.ebi.ac.uk/pride/archive/projects/PXD005590)

- Convert the spectrum files to centroided `.mzML` files and and place them in the subfolder within `data` called `spectrum_files` crux-lfq works with both profile and centroided data.

- Run the script `run_tide.sh` :

  - We assume you are using a Linux distro if not, then read the shell script and run the equivalent commands on your operating system.

  - We assume you are using `.mzML` files if not change the `extension` on line 9 in `run_tide.sh`

### FDR for Tide:

- Run the script `run_make-pin.sh`

- Once the previous script is done run `run_percolator.sh`

# CruxLFQ

- Run crux-lfq with this shell script: `run_crux-lfq.sh`

# FlashLFQ

- Owing to the fact the FlashLFQ and Crux (percolator) are not interoperable, you need to download a CLI tool that converts from percolator to crux from [crux-to-flash]().

- Once it's downloaded place this file in the `bin` folder and run the command `run_crux_to_flash.sh`

- In order to perform benchmarking tests as well as use the search results from `Tide` it is easier to use FlashLFQ as a library within [mzlib](https://www.nuget.org/packages/mzLib/). You can find the code for our FlashLFQ wrapper in `src/FlashLFQWrapper`, you can either compile the code yourself or download a compiled version from [zenodo](). \
  If you choose to compile then run `dotnet publish -c Release -r linux-x64 --self-contained true /p:PublishSingleFile=true`

- Once you have the binary file, place it in `bin` folder

# Sage LFQ

- Download sage from [here](https://github.com/lazear/sage/releases/tag/v0.14.7)

- Untar and copy the sage binary into the `bin` folder.

- Run sage with this shell script: `run_sage.sh`

- Run merge for the sage data `python merge_sage_data.py` - We run sage for single spectrum files since it's automatically run in MBR mode.

# MaxQuant & Fragpipe

- Download MaxQuant from [here](https://maxquant.org/). Follow the instructions and install it. Once done, Run the experiments ensuring it conforms to Tides params and MBR is turned off.
  MaxQuant produces a few files but for our use case we are most interested in the one uploaded to [zenodo]()

- Download Fragpipe from [here](https://github.com/Nesvilab/FragPipe/releases). Follow the instructions and install it. Once done, Run the experiments ensuring it conforms to Tides params and the MBR is turned off. Fragpipe produces a few files but for our use case we are most interested in the one uploaded to [zenodo]()

# OpenMS

- Install OpenMS based on the instructions found [here](https://openms.readthedocs.io/en/latest/about/installation/installation-on-gnu-linux.html)

- Download comet from [here](https://github.com/UWPR/Comet/releases/tag/v2025.02.0)

- Place the binary in the `bin` folder, you may need to make it executable

- Run `run_proteomics-lfq.sh`

# Data cleaning

- For the sake of convenience, organize these sepecific files into the folder `organized_results` within `results` folder.

- Copy `crux-lfq-mod-pep.txt` from CruxLFQ results folder into `organized_results`

- Copy `FlashLFQ+mods+protein_id_modpep.txt` from FlashLFQ results folder into `organized_results`

- Copy and rename `lfq.tsv` from sage results folder into `organized_results` as `sage_lfq.tsv`

- Copy and rename `peptides.txt` from MaxQuant results frolder into `organized_results` as `maxquant_peptides.txt`

- Copy and rename `combined_peptide.tsv` from Fragpipe results folder into `organized_results` as `ionquant_combined_peptide.tsv`

- Copy and rename `lfq.mztab` from proteomics-lfq results folder into `organized_results` as `proteomicslfq.mzTab`

- Run the command `python clean_data.py`

You can also download the data from [zenodo]()

# Benchmarking tests

Owing to the fact that `crux` is a CLI tool we used FlashLFQ from the mzlib library. The binary for our FlashLFQ used can be found here [zenodo]()

Download the PSM file used from here [zenodo]()

You can also just download the results from our timing tests from here [zenodo]()

# Plots

- To plot figure 1 run `python plot_fig1.py`

- To plot figure 2 run `python plot_fig2.py`

- To plot the benchmarking results run `python plot_timing_to_measure_runs.py`

- To generate the pearson correlations run `python stats.py`
