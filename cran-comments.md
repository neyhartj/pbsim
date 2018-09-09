## Test environment
* local Windows install, R 3.5.1
* ubuntu 14.04.5 (on travis-ci), R 3.5.__
* win-builder (devel and release)


## R CMD Check comments
R CMD check results
Status: 1 ERROR, 3 WARNINGs, 5 NOTEs

See
  'C:/Users/jln54/Documents/Programming/Repositories/pbsim.Rcheck/00check.log'
for details.

checking examples ... ERROR
Running examples in 'pbsim-Ex.R' failed
The error most likely occurred in:

> base::assign(".ptime", proc.time(), pos = "CheckExEnv")
> ### Name: find_proxmarkers
> ### Title: Find the position of markers near a locus
> ### Aliases: find_proxmarkers
> 
> ### ** Examples
> 
> 
> # Simulate the genome
> n.mar  <- c(505, 505, 505)
> len <- c(120, 130, 140)
> 
> genome <- sim_genome(len, n.mar)
> 
> # Sample marker names to lookup
> sample_markers <- sample(markernames(genome), size = 3)
Error in qtlnames(genome) : 
  No genetic model has been declared for the genome
Calls: sample -> markernames -> qtlnames
Execution halted

checking S3 generic/method consistency ... WARNING
print:
  function(x, ...)
print.genome:
  function(x)

summary:
  function(object, ...)
summary.genome:
  function(x)

See section 'Generic functions and methods' in the 'Writing R
Extensions' manual.

Found the following apparent S3 methods exported but not registered:
  summary.genome
See section 'Registering S3 methods' in the 'Writing R Extensions'
manual.

checking Rd cross-references ... WARNING
Missing link or links in documentation object 'calc_exp_genvar.Rd':
  'sim_crossing.block'

Missing link or links in documentation object 'sim_family.Rd':
  'subset.pop'

Missing link or links in documentation object 'sim_family_cb.Rd':
  'sim_crossing.block'

See section 'Cross-references' in the 'Writing R Extensions' manual.


checking Rd \usage sections ... WARNING
Undocumented arguments in documentation object 'qtlnames'
  'chr'

Undocumented arguments in documentation object 'sim_family'
  '...'

Undocumented arguments in documentation object 'sim_founders'
  'ignore.gen.model'

Undocumented arguments in documentation object 'sim_multi_gen_model'
  '...'

Functions with \usage entries need to have the appropriate \alias
entries, and all their arguments documented.
The \usage entries must correspond to syntactically valid R code.
See chapter 'Writing R documentation files' in the 'Writing R
Extensions' manual.
S3 methods shown with full name in documentation object 'summary.genome':
  'summary.genome'

The \usage entries for S3 methods should use the \method markup and not
their full name.
See chapter 'Writing R documentation files' in the 'Writing R
Extensions' manual.

checking top-level files ... NOTE
Non-standard files/directories found at top level:
  'README.Rmd' 'cran-comments.md' 'sandbox.R'

checking dependencies in R code ... NOTE
':::' call which should be '::': 'qtl:::mf.h'
  See the note in ?`:::` about the use of this operator.
There are ::: calls to the package's namespace in its code. A package
  almost never needs to use ::: for its own objects:
  'recombine_hypred'

checking R code for possible problems ... NOTE
calc_exp_genvar : <anonymous> : <anonymous>: no visible global function
  definition for 'dist'
calc_genoval: no visible binding for global variable 'ind'
find_proxmarkers: no visible binding for global variable '.'
find_proxmarkers: no visible binding for global variable 'chr_len'
find_proxmarkers: no visible binding for global variable 'chr'
find_proxmarkers: no visible binding for global variable 'pos'
induce_dh: no visible binding for global variable 'gen'
map_to_popvar: no visible binding for global variable '.'
... 53 lines ...
subset_pop: no visible binding for global variable 'ind'
summary.genome: no visible binding for global variable 'add_eff'
summary.genome: no visible binding for global variable 'dom_eff'
Undefined global functions or variables:
  . add_eff check_pop chr chr_len dist dom_eff effect env gen head ind
  key marker na.omit obs parent1 parent2 pheno_mean phenomean phenoval
  pos qtl_name rnorm runif tail trait var
Consider adding
  importFrom("stats", "dist", "na.omit", "rnorm", "runif", "var")
  importFrom("utils", "head", "tail")
to your NAMESPACE file.

checking Rd files ... NOTE
prepare_Rd: sim_pedigree.Rd:33-35: Dropping empty section \details

checking Rd line widths ... NOTE
Rd file 'select_pop.Rd':
  \examples lines wider than 100 characters:
     map <- lapply(split(R CMD check results
1 error  | 3 warnings | 5 notes
s2_snp_info, s2_snp_info$chrom), function(chr) structure(chr$cM_pos, names = chr$rs) )

Rd file 'sim_family_cb.Rd':
  \examples lines wider than 100 characters:
     fam_cb <- sim_family_cb(genome = genome, pedigree = ped, founder.pop = founder.pop, crossing.block = cb)
     fam_cb <- sim_family_cb(genome = genome, pedigree = ped, founder.pop = founder.pop, crossing.block = cb)

Rd file 'sim_genome.Rd':
  \examples lines wider than 100 characters:
     map <- lapply(split(s2_snp_info, s2_snp_info$chrom), function(chr) structure(chr$cM_pos, names = chr$rs) )

Rd file 'sim_phenoval.Rd':
  \examples lines wider than 100 characters:
     map <- lapply(split(s2_snp_info, s2_snp_info$chrom), function(chr) structure(chr$cM_pos, names = chr$rs) )

Rd file 'subset_pop.Rd':
  \examples lines wider than 100 characters:
     map <- lapply(split(s2_snp_info, s2_snp_info$chrom), function(chr) structure(chr$cM_pos, names = chr$rs) )

These lines will be truncated in the PDF manual.



## Downstream dependencies
There are currently no downstream dependencies for this package.