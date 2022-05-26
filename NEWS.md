# splice2neo 0.5.2

* support updated version of EasyQuant

# splice2neo 0.5.1

* update annotate_spliceai_junction() to ignore donor loss and acceptor loss when not overlapping with wild-type junctions + IR must be within transcript range
* fix minor bugs in add_context_seq() and  get_junc_pos()

# splice2neo 0.5.0

* intron rentetion events are now supported by `add_context_seq()`. the resulting context sequence covers the complete intron instead of the exon/intron boundary only. Instead of a the junction position in the cts_seq, the positions are given in form of an interval in the `cts_junc_pos`  column for intron retentions. (0,start_IR, end_IR, end_cts) 
* `add_peptide()` was adjusted for intron rententions
* more tests were added for several functions
* fix small bug in `annotate_spliceai_junction()` that led to annotation with same transcript 

# splice2neo 0.4.1

* Support donor gain and acceptor gain in introns
* export annotate_spliceai_junction and canonical_junctions()
* More stabel checks for correct format in generate_junction_id()
* Additional unit tests

# splice2neo 0.4.0

* Leafcutterr: The strand information while leafcutter parsing is now retrieved from the cluster id in the count table.
* Spladder: previous code changes missed to consider the event type in generation of the junc id for alternative splice sites. This is fixed now.
* More tests were added
* An old bug in spladder_transform_mutex_exon was fixed
* Until now a junc_id could appear multiple times if it could be part of multiple events. As we are not interested in re-constructing the events, now unique junc_ids will be returned.

# splice2neo 0.3.1

* minor update of doc and README

# splice2neo 0.3.0

* new function:
  * `generate_junction_id()`
  * `junc2breakpoint()`
  * `breakpoint2junc()`
* implementation of new junction id format
* update documentation and README
* update and add unit tests

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

# splice2neo 0.0.0.9000

* Added a `NEWS.md` file to track changes to the package.

# splice2neo 0.0.1

* Added basic package structure
* Added functions `parse_spliceai()` and `format_spliceai()`
* Automatic build of pkgdown website
