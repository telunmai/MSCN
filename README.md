# Minimum Span Clustering Network (MSCN)

This repository provides an R implementation of **Minimum Span Clustering Network (MSCN)**, an unsupervised framework for constructing multilevel co-expression modules from a square distance matrix or a similarity/correlation matrix.

MSCN initializes each node as an individual Level 0 module, denoted as `M(0)_i`, and recursively merges modules using nearest-neighbor distance relationships. At each retained clustering level, modules are named in the manuscript format `M(level)_module_index`; for example, the 317th module at Level 1 is written as `M(1)_317`.

The final global cluster, where all nodes belong to a single module, is not exported because it contains no informative module structure for downstream analyses. The final merging edges that produce this single global cluster are also excluded from the output.

## Method overview

Given an `N x N` matrix, MSCN proceeds as follows:

1. **Initialization**  
   Each node is initialized as an individual module: `M(0)_1`, `M(0)_2`, ..., `M(0)_N`.

2. **Nearest-neighbor edge search**  
   For each current connected component, the shortest edge(s) to nodes outside the component are identified. If multiple edges share the same minimum distance, all tied shortest edges are retained.

3. **Module merging**  
   Components connected by selected shortest edges are merged to form the next MSCN level.

4. **Module labeling**  
   Modules are labeled as `M(level)_module_index`, such as `M(1)_1`, `M(2)_3`, or `M(4)_1`.

5. **Stopping rule**  
   The algorithm continues until all nodes are merged into one global cluster. The final global cluster and its corresponding final merging edges are not returned.

## Input format

The input should be a square numeric matrix.

### Distance matrix input

Use `dist_matrix = TRUE` when the input is already a distance matrix.

```r
results <- minimum_span_clustering_network(my_distance_matrix, dist_matrix = TRUE)
```

Recommended distance definition for gene co-expression analysis:

```text
d_ij = 1 - |r_ij|
```

where `r_ij` is the Pearson correlation coefficient between gene `i` and gene `j`. Smaller distances indicate stronger absolute co-expression.

Expected properties:

- The matrix must be square.
- The matrix must not contain `NA` values.
- Diagonal entries should be `0`.
- Off-diagonal entries should be positive distances in the current implementation. Zero-valued off-diagonal entries are ignored by the line that filters candidate edges using `Distance > 0`.

### Similarity or correlation matrix input

Use `dist_matrix = FALSE` when the input is a similarity or correlation matrix. The function converts it internally as:

```r
adj_matrix <- 1 - abs(adj_matrix)
```

For a Pearson correlation matrix, this corresponds to `d_ij = 1 - |r_ij|`.

```r
results <- minimum_span_clustering_network(my_correlation_matrix, dist_matrix = FALSE)
```

When `dist_matrix = FALSE`, the output edge table also includes a `Similarity` column.

## Toy example

The repository includes `Distance.txt`, a 12 x 12 toy distance matrix. This matrix is intended to reproduce the schematic MSCN example.

Run the example from a terminal:

```bash
Rscript MSCN.R
```

Alternatively, run it interactively in R:

```r
sample_matrix <- as.matrix(read.table("Distance.txt", header = FALSE))
results <- minimum_span_clustering_network(sample_matrix)

results$cluster
results$edges
```

## Expected output

The function returns a list with two elements:

```r
results$cluster
results$edges
```

### `results$cluster`

`results$cluster` reports each original node and its module assignment at each retained MSCN level.

For the toy matrix, the output should be:

```text
Node      Level:1   Level:2
M(0)_1    M(1)_1    M(2)_1
M(0)_2    M(1)_1    M(2)_1
M(0)_3    M(1)_1    M(2)_1
M(0)_4    M(1)_1    M(2)_1
M(0)_5    M(1)_2    M(2)_2
M(0)_6    M(1)_2    M(2)_2
M(0)_7    M(1)_3    M(2)_2
M(0)_8    M(1)_3    M(2)_2
M(0)_9    M(1)_1    M(2)_1
M(0)_10   M(1)_1    M(2)_1
M(0)_11   M(1)_4    M(2)_1
M(0)_12   M(1)_4    M(2)_1
```

### `results$edges`

`results$edges` reports the retained node-node edges used during MSCN merging.

Columns:

- `Node1`: first node in manuscript notation
- `Node2`: second node in manuscript notation
- `Distance`: edge distance
- `Level`: numeric MSCN level at which the edge was retained

For the toy matrix, the output should be:

```text
Node1     Node2     Distance   Level
M(0)_1    M(0)_4    14         1
M(0)_2    M(0)_4     6         1
M(0)_3    M(0)_4     2         1
M(0)_5    M(0)_6    22         1
M(0)_7    M(0)_8    18         1
M(0)_2    M(0)_9     8         1
M(0)_9    M(0)_10   14         1
M(0)_11   M(0)_12   10         1
M(0)_10   M(0)_11   22         2
M(0)_6    M(0)_7    26         2
```

## Notes on the toy distance values

The values in the toy example are display-scaled distances used for illustrating the MSCN merging procedure. For example, if distances are shown as `100 x d_ij`, then a displayed value of `99` corresponds to a raw distance of `0.99`, and a displayed value of `2` corresponds to a raw distance of `0.02`.

This scaling is used only for readability in the toy example. Multiplying all distances by a positive constant preserves the rank order of distances and therefore does not change nearest-neighbor relationships or MSCN module assignments.

For empirical gene co-expression analysis, users can either:

1. compute the raw distance matrix using `d_ij = 1 - |r_ij|` and set `dist_matrix = TRUE`; or
2. provide a Pearson correlation matrix directly and set `dist_matrix = FALSE`.

## Using MSCN on your own data

Example using a distance matrix:

```r
source("MSCN.R")

my_distance_matrix <- as.matrix(read.table("my_distance_matrix.txt", header = FALSE))
results <- minimum_span_clustering_network(my_distance_matrix, dist_matrix = TRUE)

write.table(results$cluster, "MSCN_cluster_assignments.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(results$edges, "MSCN_edges.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
```

Example using a correlation matrix:

```r
source("MSCN.R")

my_correlation_matrix <- as.matrix(read.table("my_correlation_matrix.txt", header = FALSE))
results <- minimum_span_clustering_network(my_correlation_matrix, dist_matrix = FALSE)

write.table(results$cluster, "MSCN_cluster_assignments.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(results$edges, "MSCN_edges.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
```

## Environment

The MSCN function uses base R only. No additional R package is required.

## Interpretation notes

- MSCN produces a mutually exclusive module assignment at each retained level.
- The current implementation does not include a WGCNA-like grey or unassigned module.
- Each node is assigned to exactly one module at each retained MSCN level.
- The final single global cluster is not exported.

## Citation

Please cite the associated manuscript when using this code:

```text
Chen-Ling Lee, Geng-Ming Hu, Yi-Pei Li, and Te-Lun Mai*. "Mapping lncRNAs onto multilevel mRNA co-expression modules in autism."
```