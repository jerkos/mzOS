# Running mzOS as a script

## Perform analysis 

Installing mzOS made a script available in your current environnement. It is called `mzos_script`. The aim of this script is to perform an analysis given a XCMS peaklist:

* Deisotoping
* Feature annotation (adducts, fragments)
* Database search (only HMDB and LMSD are supported for now)
* Annotation confidence estimation
	* using network presence/missing metabolites information (Bayesian inference)
	* using observed/theoritical isotopic pattern

## First time running

Open a terminal in the folder where you want to perform the analysis.

Running `mzos_script` for the first time will create some configuration files you will have to modify to match your
experimental settings in the current directory:

* mzos.yml
* FRAGMENTS.csv
* ADDUCTS.csv

The first row of csv files describes their content. Some common adducts/fragments are predefined. You will have to 
adjust it. The yaml file contains the main parameters for internal algorithms such as clustering, deisotoping etc...


## Second running

Relaunching the script will launch the analysis using configuration present in the current directory.


## Parameters

	parser.add_argument("-x", "--xcms_pkl", type=str, help="path to the xcms peaklist", required=True)
    parser.add_argument("-p", "--polarity", default='negative', choices=['negative', 'positive'], help='experiment polarity', required=True)
    parser.add_argument("--mz_tol_ppm", type=float, default=10.0, help='mass over charge tolerance', required=False)
    parser.add_argument("--dims", default=False, action='store_true', help='direct infusion MS experiment', required=False)
    parser.add_argument('--db', default='hmdb', choices=['hmdb', 'lmsd', 'hmdb + lmsd'], required=False)
    parser.add_argument("--output", type=str, default="annotations.tsv", required=False)
    parser.add_argument("--bayes", default=True, required=False)
	
### --xcms_pkl

Simply the path to the XCMS peaklist

### --polarity

Polarity used to acquire spectra. Should be 'negative' or 'positive'

### --mz_tol_ppm

Tolerance in mass precision used in algorithm (in ppm). Defautlt 10.

### --dims

If the experiment is a direct infusion experiment True or False.

### --db

Database to search for. 'hmdb' or 'lmsd' or 'hmdb + lmsd'

### --ouput

Path to the result directory

### --bayes

Perform Bayesian algorithm or not (can be time consuming). True or False.

## Example

Be sure to activate the virtual environnement where you installed mzOS.

`mzos --xcms_pkl C:\Users\M\xcms_results.tsv --polarity 'negative' --db 'hmdb + lmsd' --output '.\annotations.csv'`

## Result file

Matrix with several columns:

* **id**: id of the feature.

* **mz**: mass over charge

* **time**: time elution

* **Main putative attribution**: the most probable tag that can be assigned to this feature. It can be:

	* a *monoistope + n isotopes + n adducts/fragments*
	
	* an isotope (C13 ...)
	
	* an adduct
	
	* a fragment
	
	* nothing


* **Main attribution pattern composition**: one of the most important column, showing all relations between features. For example:

	* feature detected as an isotope in *Main putative tag* shows *Isotope C13 of 2546 for charge=1*. Here you can find to which feature it has been detected as a C13 isotope. You have also a charge information.

	* feature detected as a monoisotope: 

* **Putative secondary attributions**: All possible other attributions (or tag) that algorithm found but set as secondary (main is considered as the most probable)

* **Putative annotation**: matching metabolite(s)

* **Putative formula**: chemical formula of the matching metabolite

* **Inchi**: inchi metabolite formula

* **Database_id**: Database Id such as HMDB_ID and KEGG_ID

* **Isotopic pattern matching score**: RMSD between observed and theoritical isotopic pattern intensities

* **Annotation assignment probability**: Using network analysis to infer a presence probability using a bayesian algorithm (see metsamp)

* **Annotation pattern composition**: Isotopic pattern detected to use the *Isotopic pattern matching score*


