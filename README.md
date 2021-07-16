## Robust earthquake spectra from envelope inversion

This repository contains the source code for the reproduction of results and figures of the following publication:

Tom Eulenfeld, Sebastian Heimann, Torsten Dahm and Ulrich Wegler (2021),
Fast and robust earthquake source spectra and moment magnitudes from envelope inversion,
submitted to *Bulletin of the Seismological Society of America*. <!---[[preprint]()]-->

Additionally, this repository contains result files and figures.

#### Preparation

1. Download or clone this repository.
2. Install the the relevant python packages: `python>=3.7 qopen>=4.1 grond>=1.4 obspy>=1.2 pyrocko>=2021.06.29 matplotlib>=3.4 utm mtspec`. The minimum version also specifies the versions I used.
   The following is working for me when using conda

```
conda --add channels conda-forge
conda create -n eqspec obspy statsmodels utm mtspec
conda activate eqspec
conda install -c pyrocko pyrocko
pip install qopen grond
```

3. Prepare data: The data is hosted on zenodo with doi:[10.5281/zenodo.3741464](https://www.doi.org/10.5281/zenodo.3741464). Copy the data into a new `data` folder.
   Finish preparing the data with the script `python prepare_data.py`. The `data/waveforms` folder can be deleted afterwards.

#### Running Qopen

Switch to the `qopen` folder with the Qopen configuration and results. Check out the `run_qopen.sh` script which runs the Qopen commands described in the article (figure 3).
Run these one by one or the script itself to recreate the results.

#### Running Grond

First download the [used Green's function database](https://greens-mill.pyrocko.org/vogtland_scatter_v4-f87b40) and insert the directory where it is stored in the file `$HOME/.pyrocko/config.pf` in the `gf_store_superdirs: [your_gf_path]` field.
Switch to the `grond` folder with the Grond configuration already in place. Run the commands in the script `run_grond.sh` or the script itself.
Optionally the `--parallel` option can be used, please check the help with `-h` flag.

#### Calculate spectra from direct onsets

Switch to the `scripts` directory and run

```
python calc_onset_spectra_and_plot_comparison.py
```

#### Reproduce figures

Switch to the `scripts` directory and run all the `plot*.py` python scripts. The figure with Grond fits is created by the bash script `montage_grond_fits.sh` (montage program needed).

#### Results

The above scripts repopulate the folders `figs`, `grond`, `onsets` and `qopen`.
The results for Qopen are located in the `qopen` folder. Grond results can be viewed [here](https://data.pyrocko.org/publications/grond-reports/west-bohemia-2018/) <!---([doi]())-->.
Results from the spectra of direct onsets are located in the `onsets` folder. All created figures are stored in the `figs` folder.
