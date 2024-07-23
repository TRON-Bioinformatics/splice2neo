# splice2neo 0.6.11

- fix of liftover function 
- improve warning and add more tests
- update of combine_mut_junc() to accept empty tables as input
- update of format_cispliceai() to optionally keep gene_ids
- add more test

# splice2neo 0.6.9

- Simplification of function `annotate_mut_effect()` and ensuring that return of same column names in case of empty output


# splice2neo 0.6.8

- Add function `liftover_junc_id() to liftOver splice junctions to other 
 reference genomes or personalized genomes


# splice2neo 0.6.7

## Added:

 - Functionality from splice2neo_neoants:
    - CI-SpliceAI parsing and formatting functions with example data and tests
    - IRFinder parsing and filtering functions with example data and tests
    - StringTie parsing and formatting functions
    - SUPPA2 parsing and formatting functions
    - STAR parsing and formatting functions with example data and tests
    - RegTools parsing and formatting functions with example data and tests
 - add test for empty CI-SpliceAI input file
 - add more global variables to avoid notes in R CMD check

## Changed:

 - Update CI/CD pipeline for less verbosity during R pkg installations
 - Return single-row data.frame for empty input in `annotate_mut_effect()`
 - include 'number_of_supporting_reads` when parsing SplAdder output
 - generalize `generate_combined_dataset()` to multiple inputs
 - Update README with list of supported tools
 - add option consider_intron_retention to `annotate_mut_effect()`
 - adjust syntax for dplyr >= 1.1.1
 - remove rarely used pkgs from Imports
 - simplify code in parse_spliceai_thresh()
 - minor spell check fixes and doc updates

## Fixed: 

 - fix warnings due to NA's in CI-SpliceAI parsing function
 - fix use of `generate_combined_dataset()` in vignette

# splice2neo 0.6.6

* add vignette with example workflow
* fix parsing of pangolin results
* simplify README

# splice2neo 0.6.5

* Fixed bug while parsing VCF files containing only a single variant
* Improved error handling of spladder import

# splice2neo 0.6.4

* add new function to select more relevant junction-transcript combinations 
* adapt input and output of exon_in_intron()

# splice2neo 0.6.3

* revision of parse_pangolin and format_pangolin functions 
* annotate_mut_junc should remove predicted effects outside of transcript range
* add_peptide returns full protein sequence until first stop codon

# splice2neo 0.6.2

* ensure correct exon usage in MMSplice annotation
* two new functions to filter for unique junctions based on mut_id, junc_id and tx_id based on the highest effect score

# splice2neo 0.6.1

* annotate_mut_effect() now optionally considers only transcripts from genes provided by SpliceAI or Pangolin
* Column `class` was renamed to `event_type`
* combine_mut_junc() returns now rows that are unique based on `mut_id`, `junc_id`, `tx_id`, `event_type`
* update of README
* update of description of add_peptide()

# splice2neo 0.6.0

* major update of add_peptide() function + adding some tests
* minor bug fix in exon_in_intron() function

# splice2neo 0.5.6

* minor changes to support stringr >= 1.5.0

# splice2neo 0.5.5

* integrate splicing mutation tool *Pangolin*
* generalize mutation effect annotation 

# splice2neo 0.5.4

* adds a function to annotate if a exon of another transcript is located in an intron 
* update README  

# splice2neo 0.5.3

* adapts splice2neo to easyquant 0.4.0. (https://github.com/TRON-Bioinformatics/easyquant)
* prevents annotate_spliceai function from failing in case of empty input adds some more tests

# splice2neo 0.5.2

* support updated version of EasyQuant

# splice2neo 0.5.1

* update annotate_spliceai_junction() to ignore donor loss and acceptor loss when not overlapping with wild-type junctions + IR must be within transcript range
* fix minor bugs in add_context_seq() and  get_junc_pos()

# splice2neo 0.5.0

* intron retention events are now supported by `add_context_seq()`. the resulting context sequence covers the complete intron instead of the exon/intron boundary only. Instead of a the junction position in the cts_seq, the positions are given in form of an interval in the `cts_junc_pos`  column for intron retentions. (0,start_IR, end_IR, end_cts) 
* `add_peptide()` was adjusted for intron retentions
* more tests were added for several functions
* fix small bug in `annotate_spliceai_junction()` that led to annotation with same transcript 


# splice2neo 0.4.1

* Support donor gain and acceptor gain in introns
* export annotate_spliceai_junction and canonical_junctions()
* More stable checks for correct format in generate_junction_id()
* Additional unit tests

# splice2neo 0.4.0

* Leafcutter: The strand information while leafcutter parsing is now retrieved from the cluster id in the count table.
* Spladder: previous code changes missed to consider the event type in generation of the junc id for alternative splice sites. This is fixed now.
* More tests were added
* An old bug in spladder_transform_mutex_exon was fixed
* Until now a junc_id could appear multiple times if it could be part of multiple events. As we are not interested in re-constructing the events, now unique junc_ids will be returned.

# splice2neo 0.3.1

* minor update of doc and README

# splice2neo 0.3.0

* update of the junction id format: chr:start-end:strand
* new function:
  * `generate_junction_id()`
  * `junc2breakpoint()`
  * `breakpoint2junc()`
* The recent versions of leafcutter (v0.2.9) and spladder (v3.0.2) will be supported from now on.
* update and add unit tests
* update documentation and README
* visualization of modified junctions added 

# splice2neo 0.2.0

* new function:
  * `bed_to_junc()`
  * `is_canonical()`
  * `is_in_rnaseq()`
  * `parse_gtf()`
* updated functions:
  * `add_context_seq()`
  * `add_peptide()`
  * `add_peptide()`
  * `get_junc_pos()`
* faster implementation of `modify_tx()`
* update documentation and README
* update and add unit tests

# splice2neo 0.1.3

* This minor release adds a pseudo example to the README + simplifies the `transform_for_requant` function

# splice2neo 0.1.2

* CI config for github action
  * R CMD check
  * test coverage with Codecov
  * pkgdown website on github.io
* update README

# splice2neo 0.1.1

* clean up code and doc

# splice2neo 0.1.0

* New implementation of transcript modifications by junction
* New user interface functions that mainly work with data.frames: 
  - `add_tx()`
  - `modify_tx()`
  - `add_context_seq()`
  - `add_peptide()`
* Updated toy example data sets
* First example in README

# splice2neo 0.0.2

* inclusion of functions from previous repositories

# splice2neo 0.0.1

* Added basic package structure
* Added functions `parse_spliceai()` and `format_spliceai()`
* Automatic build of pkgdown website
