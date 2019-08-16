# ANOR 0.2-1
1. Fixed phylogram when Map has multiple columns (Issue https://github.com/surh/AMOR/issues/28).

# AMOR 0.2-0
1. Added support in subset.Dataset for when only one sample remains.
2. Added travis CI build testing (Issue https://github.com/surh/AMOR/issues/6)
3. Fixed bug with plotgg_rankabun caused by update of pool_samples. Added
option to force return of a dataset object from pool samples.
4. Fixed bug with plotgg_rankabun2 caused because some columns of the
aggregated dataset became tables instead of vectors.
5. Moved functions out of AMOR.r into their own files and removed
unecessary functions. Might break backwards compatibility.
6. Removed functions that identify samples for rarefaction.
This has been deprecated since rarefaction is not recomended
to deal with uneven sampling. Might break backwards compatibility.
7. Eliminated summarizeOTUDistribution. Never implemented.
8. Removed normalizeSample function. Might break backwards compatibility.
9. Roxygen documentation and examples for get_tax_level.
10. Removed read.pool function. Might break backwards compatibility
11. Documentation and examples for write.qiime.
12. Roxygen2 documentation for bootglm and associated methods.
13. Homgenizing method variable names for PCA and PCO. Might break
backwards compatibility.
14. Moving documentation to roxygen2 for functions beta_diverstity,
clean, create_dataset, find_good_samples, findGOODSamples, matrix_glm,
matrix_glmNB, plotgg_rankabun2, plotgg_rankabun, phylogram,
plotgg_taxon, plotgg_var, rarefaction, remove_samples,
remove_taxons, site_diversity, compare_site_diversity and
total_richness. (Issue https://github.com/surh/AMOR/issues/4)
15. Switched theme_blackbox to a function and added roxygen documentation.
Might break backwards compatibility.
16. Added codecov support (Issue https://github.com/surh/AMOR/issues/7)

# AMOR 0.1-2
1. Fixed bug introduced on previous release regarding in remove_samples.

# AMOR 0.1-1
1. Fixed bug in remove_taxons, it now works correctly when only one
sample or taxon remains.
2. Fixed bug in remove_samples, it now works correclty when only
one sample or taxon remains.

# AMOR 0.1-0
1. Fixed error in summary.PCO. It required an object called Dat.pco
2. Added fill option to plotgg. Included in documentation
3. Added option to collapse_by_taxonomy with a custom factor
4. Fixed bug in collapse_by_taxonomy.Dataset so now it can take no group
5. collapse_matrix now uses plyr instead of apply functions. This allows for
consistent return of matrix with dimension names when the grouping factor
contains one level (Issue https://github.com/surh/AMOR/issues/10). This breaks
backwards compatibility. Roxygen2 documentation for collapse_matrix
7. Updated pool_samples. Changed naming conventions of parameters,
and removed the option to return a Dataset object from default method. These
changes break backwards compatibility. Pool samples can also now take a vector
or grouping factor. Updated and moved documentation to
roxygen2 and updated examples (Issue https://github.com/surh/AMOR/issues/11).

# AMOR 0.0-15
1. Adding samples & taxa functions
2. Fixed glitch in clean when only one variable in map
3. Mainteinance for heatgg function. Now there is an appropriate class and print method for when
clusering is used. Also documentation was moved to roxygen and updated. Old code was removed.
4. plotgg mainteinance. Now plots proportion of variance in the indicated axis, and value of PCO
function now has dimension names in the points variable.
5. Created variables function, and updated documentation for samples & taxa.
6. Mainteinance to the basic Dataset functions. Moved some functions from AMOR.r to create_dataset.r,
and added a Dataset print method.
7. Added simplify option to read.am
8. Replaced findGoodOTUs with measurable_taxa function.

# AMOR 0.0-14
1. Added space option to phylogram
2. Set width in geom_bar of phylogram to 1.
3. Added summary and print methods for PCA and PCO
4. Continued switching documentation to roxygen
5. Added README

# AMOR 0.0-13
1. Added subset.Dataset method.
2. Changed clean function, now it cheks that all the entries in a row/column are equal to zero to remove. Slighlty slower but more general.

# AMOR 0.0-12
1. read.am always returns a Dataset object. Documented
2. Added Dataset method for write.qiime and homogenized variable names. Documented

# AMOR 0.0-11
1. fixed glitched where MASS package was not being loaded and one couldn't use the matrix_glmNB function.
2. remove_samples now applies droplevels by default, and can be turned of by the user.
3. Added functions plotgg_var and plotgg_taxon which make boxplots of a variable in the mapping file or a taxon in the abundance matrix.
4. Now the pool_samples function return Dataset objects
5. Added plotgg_rankabun function and Dataset methods.
6. Also added experimental function plotgg_rankabun2
7. Improved phylogram function to return only plot and now it can also take an argument for showing only the top n taxa, and collapsing the rest. Also added Datast method for it.
8. Fixed glitch in PCA function. Totvar was not being calculated
9. Fixed glitch in remove_samples functiin, now it can handle Dataset onjects with NULL map.
10. Added ggdendro dependency which is required for the heatgg function
11. Addedd heatgg function with support for sample faceting and 2D clustering.
12. Added site_diversity functions and a plotgg method for it.

# AMOR 0.0-10
1. fixed bug in pool_samples. The test set was always being used.
2. read.am returns a dataset object if taxonomy is passed.
3. Fixed glitch in collapse_by_taxonomy help

# AMOR 0.0-9
1. Changed the PCA function, now it does not call pca from labdsv, but prcomp directly.
2. Changed PCO function, now it does not call pco, but cmdscale directly.
3. Fixed plotgg.PCA so that it returns all the data (including biplot) in the main object
4. Created PCO method for plotgg
5. Created clean, remove_samples and remove_taxons functions and docummented them.
6. Added Rhizo.tax as data distributed with the package
7. Added collapse_matrix functiona and documented it
8. labdsv now is suggested but not required. vegan is suggested as well.
9. Added collapse_by_taxonomy and documentation.
10. Added pool_samples and documented
11. Updated create_dataset to make both Map and Tab optional

# AMOR 0.0-8
1. Added matrix_glmNB function
2. Change pca and pco functions to PCA and PCO so they don't interfere with labdsv.
They are just wrappers for the labdsv functions. Documented both.
3. Function PCA was modified so that the pca object returned by its Dataset method contains
Map and Tax slots.
4. plotgg generic function and plotgg.pca method created
5. Changed rrarefy function to rarefaction, and also altered so the columns are samples
now. Created S3 method for Dataset objects and documented.

# AMOR 0.0-7
1. Important mainteinance to matrix_glm and bootstrap_glm.
2. Docummentation improved (read.am, write.qiime)
3. Added data Rhizo and Rhizo.map
4. Added PCA and PCoA methods for Dataset objects. Relies ehavily in the labdsv package
5. Updated matrix_glm to take ... argument and pass it to glm. bootstrap_glm was modified
accordingly.
6. Updated matrix_glm to take return design matrix X.
7. Updated docummentation

# AMOR 0.0-6
1. Added create_dataset() function which defines the "Dataset" class
2. Added the matrix_glm methods, and documentation. Includes a summary
method and a print method for the summary.
3. Added the bootstrap_glm methods, and documentation. Includes a summary
method and a print method for the summary.
4. Expanded documentation

# AMOR 0.0-5
1. Fixed bug for specifying the number of rows in the legend.

# AMOR 0.0-4
1. Reorganized functions into files according to general purpose.
2. Added an option to specify the maximum number of rows in the legend
of the function phylogram.
3. Added function beta_diversity, which calculates Bray-Curtis dissimilarity
so we don't need the vegan function anymore.
4. Removed vegan from the list of dependencies.