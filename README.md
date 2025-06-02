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
