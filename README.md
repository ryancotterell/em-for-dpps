# Code Installation #

Start MATLAB from the trunk/code directory and compile the one java
file: "javac EqualByValueIntArray.java".


# Code Usage #

Start MATLAB from the trunk/code directory, or run startup.m to ensure
that the current folder is on the java classpath.

## Baby registry experiments: ##

To re-produce the baby registry experiments from the paper, run
registries\_test.m with the appropriate inputs.  For example:
```
  registries_test('../data/1_100_100_100_safety_regs.csv', ...
                         '../data/1_100_100_100_safety_item_names.txt', ...
                         '../safety_first.mat', 0.7, 10, 'wi', -1)'
```
will run ten trials testing the safety product category, with 70% of
the data allocated for training, and a Wishart type initialization.
The results will be saved in 'safety\_first.mat'.  Run
registries\_eval.m to print text describing the results:
```
  registries_eval('safety_first.mat');
```
When examining these results, note that the "true" log likelihood here
is just the KA one; we have aliased the two terms for the registries
test for simplicity, as no "true" likelihood exists.

# Baby Registry Data #

The folder trunk/data contains the results of filtering data from
close to 30,000 of Amazon's baby registries pages.  Each file refers
to a specific category of baby item (feeding, diaper, safety, etc.).
The file naming reflects the filtering criteria as follows:

  * Each registry obeys:
    * first # = minimum allowed number of items within category x
    * second # = maximum allowed number of items within category x; if the registry size is not within this range, it is discarded

  * Each category obeys:
    * third # = maximum allowed number of products per category (across all registries); if this limit is reached, the least frequent items are discarded in order to maintain the limit
    * fourth # = minimum allowed number of occurrences of a product (across all registries); products with lesser frequency are discarded

Notice that some of these criteria interact with each other (e.g. a
product may occur in enough registries to pass the minimum occurrences
cutoff, but only prior to the registries are filtered by size).  Thus,
we iterated the filtering process until all criteria were
simultaneously satisfied.

There are two types of files for each setting: item\_names and regs.
Item names lists the names of all N products.  Each line gives the
line number, then the item name:

  * 1 safety: Graco Sweet Slumber Sound Machine, White : Electronic Infant Sleep Aids : Baby
  * 2 safety: Snuza Baby Monitor, Hero : Baby
  * 3 safety: Summer Infant Baby Touch Digital Color Video Monitor : Baby
> ...

Note that while some items have colors (e.g. the "White" in item 1
above) or other very specific details associated with them, this does
not necessarily mean that every registry listed as containing that
item actually contained an item of that precise specification.
Rather, they are guaranteed to have had an item similar enough that
its Amazon product ID is identical.

Each line of a regs file constitutes one registry (restricted to the
category indicated by the filename).  The numbers index into the
item\_names.txt file.  For example, if the first line in the csv file
is "1,2" then this indicates the first two products from the
corresponding item\_names.txt file.