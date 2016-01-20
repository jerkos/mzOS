# Welcome to mzOS project

**mzOS** is a small library to perform deisotoping and feature annotation from an extracted feature list (commonly designed by _peaklist_) in LC-MS (Liquid Chromatography coupled to mass spectrometry) datasets.

Basically, you run first XCMS (using the classical pipeline), giving as a result a peaklist file. This is this result file that mzOS uses. Other feature extraction softwares might be supported later (OpenMS). 

## Installation

mzOS install should be pretty straightforward. The last version at the moment is  0.1.3. Simply using

* `pip install mzos`

You can also download a zip archive from [github](http://github.com/jerkos/mzOS) for example. Unzip the archive then

* `cd mzOS`
* `python setup.py install`

The project is quite heavy (around 20Mo), actually it comes with LMSD and HMDB databases.

It is still in alpha stage and is not bug free.

## Project layout

    scripts/    	# Scripts used to rebuild database sqlite files.
    third_party/	# small utility to calculate efficiently theoritical isotopic pattern
    	emass/
    		... 	# emass files
    ressources/		# contains databases, isotopes and adducts list
    
    ...       		# mzOS python modules
