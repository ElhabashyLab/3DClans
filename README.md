# Clans-3D

Clans-3D is a bioinformatics tool designed to create structural clans files from protein structure data.
It supports multiple input formats

- FASTA
- A2M
- TSV

It uses various structural similarity tools to compute pairwise structural similarities:

- Foldseek
- TM-align
- US-align

The resulting clans files can be used for further analysis and visualization of protein structure relationships with the CLANS software.

## Installation

### Requiremets

- python >= 3.8
- Java Runtime Environment (JRE)

### External Tools

The following external tools need to be installed and accessible in your system's PATH:

- foldseek
- TMalign
- USalign

### Creating a Virtual Environment and Installing Dependencies

Clans-3D can be installed by cloning the repository and installing the required dependencies.

```bash
cd Clans-3D
```

create a virtual environment:

```bash
python -m venv venv
source venv/bin/activate  # On Windows use `venv\Scripts\activate`
```

install dependencies:

```bash
pip install -r requirements.txt
```

## Usage

Clans-3D can be run via the command line interface (CLI) provided in `main.py`.  
The resulting clans file will be saved in the `clans_files` directory.

### Command Line Arguments

1.  `--load/-l` : Path to the input file.

    Provide a file path to a FASTA, A2M, or TSV file you want processed.
    The CLI will copy this file into `input_file_storage` as `basename.<input_type>`.  
    Example: `--load path/to/my_input.fasta`

2.  `--input_type/-i` : Type of the input file.

    Accepted values: `fasta`, `tsv`, `a2m`.  
    This tells the program how to parse the input.

    - **FASTA**: expects a standard FASTA file with sequence identifiers and sequences.  
      The header should specify the region of interest using the format: `>identifier/start-end`  
      Example: `>sp|P49811|MYOD1_PIG/2-300`

    - **A2M**: expects an aligned FASTA file in A2M format.
      The header should specify the region of interest using the format: `>identifier/start-end`  
      Example:  
      `>C5BN68_TERTT/528-601`

    - **TSV**: expects a tab-separated values file with columns: `entry, region_start, region_end`  
      Example:

          ```
          entry	region_start	region_end
          P49811	2	300
          C5BN68	528	601
          ```

3.  `--tool/-t` : Structural similarity tool to use.

    Choose the tool you want: `foldseek`, `tmalign`, `usalign`

    - _foldseek_: fast
    - _TMalign/USalign_: slower, more accurate

    Example: `--tool foldseek`

4.  `--score` : Which score to use for Foldseek.

    Accepted values in the: `evalue` or `TM`.  
    Only meaningful if `--tool foldseek`. If you pass `--score` while `--tool` is not foldseek, the program raises an error. Default is `evalue`.  
    Example: `--score evalue`

### Examples:

```bash
python main.py --load <path/to/input.fasta> --input_type fasta --tool foldseek --score <evalue>
```

### Help

To get an overview of all command line options, run:

```bash
python main.py -h
```

### Further information for usage

If you pass `--score` but `--tool` ≠ foldseek, the code raises a ValueError stating The foldseek score (`-s --score`) can only be specified if foldseek is chosen as tool.

If no region is specified in the FASTA/A2M header or TSV file, the full length of the structure is used.  
If region_start/region_end region is specified, the region_end/region_start must be specified as well.  
A Fasta input with specified regions expects the sequences to be cut to the specified regions.

Saved input files go to `input_file_storage` as `basename.<input_type_value>` and a cleaned version of the input file is saved in the same directory as `basename_cleaned.<input_type_value>`(only contains structures which were found and downloaded).

The resulting CLANS files are saved into the `Scores_Evaluation/clans_files*` directory.

## Benchmark

The Benchmarking pipeline is implemented to evaluate how the different tools can process large datasets to create clans files.

This is done by running the code of `benchmark_tool_speed.py` on a predefined set of protein structures in a folder or by supplying a valid input-file. The Script measures and records the time taken for each tool.

## Dataset-Generator

The Dataset-Generator module is used to create datasets for benchmarking and evaluation of Clans-3D.
The datasets are designed to have known cluster structures, allowing for assessment of the clustering performance of Clans-3D.

The generated fasta files consist of `n` sequences. The sequences are grouped into `n_cl` clusters.
Each cluster is defined by a `seed` sequence and contains very similar sequences.
The remaining sequences of a cluster are obtained by blasting the `seed` sequence.

## Scores-Evaluation

The Scores-Evaluation module is used to evaluate the quality of the structural clans files created by Clans-3D.

For a given dataset, it computes clans files based on structure and sequence similarity scores.  
The files are then run by clans.jar in headless mode.
The resulting clans files with updated coordinates are then compared:

- nummerically
- graph wise
- cluster wise

## Recovered Clans (recovered_CLANS)

The Recovered Clans module defines methods in `utils_old_clans.py` to run the old clans.jar in headless-mode.  
This is useful for comparing the clans files generated by Clans-3D with those containing sequence-based similarity scores.
The headless-mode ind restricted to only clustering clans files for a given number of rounds.

Recovered Clans can be run by using `run_clans_headless`.  
If you want to run recovered clans on structural similarity based clans files, set `input_file_type` to `InputFileType.CLANS`.  
If `input_file_type` is set to `InputFileType.FASTA`, recovered clans will be run using sequence similarity scores computed by BLAST based on the input fasta files.
