# CruxLFQ-publication-experiments

In order to reporduce the experiments in "Label-free quantification in the Crux toolkit" follow the steps below:

## Assumptions

- You have basic familiarity with the command line.
- You have access to a computer running Windows, macOS, or Linux. (Some experiments were run on Linux and others on Windows)
- You have sufficient permissions to install and run software.
- Required dependencies (such as Python, if needed) are installed.
- Code is found in `src`, hence each time an instruction is provided to run it is assumed the user is running from `src`.

# Database Searching:

## Tide:

- Download Crux version 4.3 for your specific operating system from the [releases page](https://github.com/crux-toolkit/crux-toolkit/releases/tag/crux-4.3).

- Extract the the downloaded zip file into the `bin` folder, if you are on linux, you can use the command:

        unzip crux-4.3.Linux.x86_64.zip -d bin

- Download the fasta file from [ProteomeXchange - PXD003881](https://ftp.pride.ebi.ac.uk/pride/data/archive/2018/05/PXD003881/Human_ecoli_trypsin_1501v_uniprot_sprot.fasta) and place the downloaded file in the folder `data`

- Download the spectrum files from [ProteomeXchange - PXD005590](https://www.ebi.ac.uk/pride/archive/projects/PXD005590)

- Convert the spectrum files to centroided `.mzML` files and and place them in the subfolder within `data` called `spectrum_files` crux-lfq works with both profile and centroided data.

- Run the script `run_tide.sh` :

  - We assume you are using a Linux distro if not, then read the shell script and run the equivalent commands on your operating system.

  - We assume you are using `.mzML` files if not change the `extension` on line 9 in `run_tide.sh`

### FDR for Tide:

As at the time of this publication crux's make-pin did not return retention time, hence download a custom compiled version from [zonodo]() for your experiments. Note that this will e added in future versions of crux.

- Download custom crux build : [crux-percolator]()

- Place the downloaded file in the `bin` folder

- Run the script `run_make-pin.sh`

- Once the previous script is done run `run_percolator.sh`

# CruxLFQ

- In order to run `lfq` in crux, download a custom build of crux: [crux-lfq](). This is because at the time of writing lfq had not yet been merged into crux main.

- Once done downloading, place the binary in `bin` folder

- Run crux-lfq with this shell script: `run_crux-lfq.sh`
