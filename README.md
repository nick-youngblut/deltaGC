# deltaGC -- scripts for simulating GC content and buoyant density of DNA fragments in a CsCl gradient

__Main application: DNA-SIP (stable isotope probing)__


## DESCRIPTION

The notebooks describe the analyzses conducted in the paper: __Youngblut,N.D. and Buckley,D.H.
Intra-genomic variation in G+C content and its implications for DNA Stable Isotope Probing (DNA-SIP).
*In prep*__

These notebooks are in the ./ipynb directory. html versions of each notebook
can found in the ./ipynb/html directory

The perl scripts & packages in this repo were used in the analyzses.
They are in the ./bin and ./lib directories (see PERL SCRIPTS below).

A list of the accession numbers for all bacterial and archaeal genomes
used in the analyzes can be found in the 'accessions' directory.


## PERL SCRIPTS 

These scripts are design to assess how GC content of DNA fragments governs the buoyant density distribution
of these fragments in a CsCl gradient. 

An important consideration is that only fragments containing the template (amplicon, such as 16S rRNA amplicons
or shotgun reads for (meta)genomics). Therefore, amplicons or shotgun reads are simulated with GRINDER, then
the template fragment that each read originate from is simulated and the GC content is calculated, which
can be used to determine the buoyant density of the molecule in a CsCl gradient (assuming equilibruim).


### PREREQUISTS

* [GRINDER](http://sourceforge.net/projects/biogrinder/ "GRINDER")

#### perl modules:

All of these (except possibly bioperl) should be easy to install with cpan or cpanm
(at least on *nix systems).

* Getopt::Euclid
* Parallel::ForkManager
* MCE
* Math::Random
* Math::Gauss
* List::MoreUtils
* Memory::Usage
* Statistics::Descriptive
* Regexp::Common
* Term::ProgressBar

__BioPerl modules__

* Bio::SeqIO 
* Bio::DB::Fasta


### INSTALLATION

To install this module, run the following commands:

	perl Build.PL
	./Build
	./Build test
	./Build install

### SUPPORT AND DOCUMENTATION

After installing, you can find documentation for each script (eg., SCRIPT.pl)
by typing: 

	SCRIPT.pl --man
	   OR
	perldoc SCRIPT.pl


### LICENSE AND COPYRIGHT

Copyright (C) 2014 Nick Youngblut

This program is free software; you can redistribute it and/or modify it
under the terms of either: the GNU General Public License as published
by the Free Software Foundation; or the Artistic License.

See http://dev.perl.org/licenses/ for more information.
