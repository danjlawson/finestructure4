# Finestructure Version 4

## About fineSTRUCTURE

fineSTRUCTURE is a fast and powerful algorithm for identifying population structure using dense sequencing data.  By using the output of ChromoPainter as a (nearly) sufficient summary statistic, it is able to perform model-based Bayesian clustering on large datasets, including full resequencing data, and can handle up to 1000s of individuals. Full assignment uncertainty is given.

finestructure works on Linux and Mac. You can also compile it for Windows if you set up the required toolchains, but you need to know what you are doing and be comfortable with the command line.

A Stochastic optimization routine is available for performing faster EDA and dealing with larger datasets - see FAQ under "What if my dataset is too big for MCMC".

Important Note: fs4 includes ChromoPainter, which have different licences and authors. Both are free for Academic use only, and explicitly exclude commercial applications. See the file COPYING in the download for details.

Please consider [Registering](https://forms.office.com/Pages/ResponsePage.aspx?id=MH_ksn3NTkql2rGM8aQVGwBpDrkt7zVLlSvbqowMvq1UQzNXV1hMQkE1QzBNTDJPS0JEWFE4TDFaVi4u), if you have not already done so.

## Technical details:

NB Installation instructions are below: Linux and Mac OS X should work fully.

### LICENCE

* The "fs" code is for non-commercial purposes only.
* It is free to use for Academic, personal and non-profit purposes.
* Attribution to software creator must be made:
    	* Acknowledgement in personal or other non-commercial work. 
		* For academic use, citation of the appropriate article(s), currently: [Lawson, Hellenthal, Myers & Falush 2012, PLoS Genetics e1002453 "Inference of population structure using dense haplotype data"](https://journals.plos.org/plosgenetics/article/info%3Adoi%2F10.1371%2Fjournal.pgen.1002453).
* For commercial licencing, contact the author dan.lawson <at> bristol.ac.uk.
* We are likely to agree to modification and distribution for non-commercial use, but please ask.
* fs makes use of the non-free software ChromoPainter version 1, which has been modified with permission of the author Garrett Hellenthal (ghellenthal <at> gmail.com). You are bound by the terms of that software SEPARATELY. It uses the same licence as fs but the copyright belongs to Garrett Hellenthal and permission to use the software outside of these terms must be arranged with both authors separately.

See LICENCE for further details.

### REQUIREMENTS

The requirements are a "recent" version of the GCC compiler, for the GLIBC libraries. These are typically present by default.

There are probably other requirements which are present on most systems by default.  If you encounter any problems let us know. 

### INSTALLATION

A binary for linux and mac is provided, which you may just be able to use. However some users have different versions of some core C++ libraries that make compilation desirable.

fs4 comes with an **installation script** `fs_install.sh`. This is optional but will help novices get the correct version for their computer as well as adding it to their PATH variable.

### COMPILATION 

If you have a similar enough build-chain to me you can use:

```{sh}
./configure
make
make install #optional
```

If you need to reconfigure the configure options, to change the 

If you are compiling on a Mac, the default C compiler (clang) does not support multi-threading.  If you would like to use another C compiler that is not the default, you can use:

```{sh}
./configure CXX=g++-11
```

(substitute your CXX version). This is implemented in ./configure.mac so that you don't have to remember each time.

### Dependencies

Dependencies for the command line version are (from version 4.1.0) nearly non-existent!

* Linux:
  1. GCC c++ compiler (package build-essential in ubuntu)
  2. Automake (optional)
* Mac OS: Either:
  1. Xcode Clang compiler `xcode-select –install` (See [e.g. This Howto](https://mac.install.guide/commandlinetools/4.html); but note that this does not enable parallel computation outside of HPC mode.)
  2. or GCC c++ compiler (package gcc in brew) for full functionality.
 
### FURTHER INFORMATION

You need to prepare your data in CHROMOPAINTERv2 format. This is nearly the same as PHASE format, but HAS BEEN UPDATED since chromopainter v 0.0.5. We provide some tools for this in the scripts directory.

Run "fs" for help. Examples are included in the "examples" directory; it is recommended to work through the examples to establish how to run this program on your own dataset.

IMPORTANT: If you have a small dataset, you can run 
"fs project -phasefile <data> -recombfile <recombination map file> -go"
and it will do everything with default settings that should work. BUT IF YOU HAVE A LARGE DATA SET THIS IS GOING TO TAKE A LONG LONG TIME. You will want to parallelise the work, which is done with the "-hpc" flag. You can then run commands on your own HPC or in parallel on a multi-core machine.

The examples show you how to do this!

[www.paintmychromosomes.com](www.paintmychromosomes.com) is the place to start for help.  Visit the FAQ page for standard issues.

fs4 (finestructure) is written by Daniel Lawson (dan.lawson@bristol.ac.uk) COPYRIGHT University of Bristol 2022.
