sample_matrix <- as.matrix(read.table("Distance.txt", header = FALSE))

# Minimum Span Clustering Network (MSCN)
minimum_span_clustering_network <- function(adj_matrix, dist_matrix = TRUE) {
  
  # Check the input matrix.
  if (!is.matrix(adj_matrix) || nrow(adj_matrix) != ncol(adj_matrix)) {
    stop("adj_matrix must be a square matrix.")
  }
  if (any(is.na(adj_matrix))) {
    stop("adj_matrix contains NA values. Please clean the data.")
  }
  
  # Convert a non-distance matrix to a distance matrix if needed.
  if (!dist_matrix) {
    adj_matrix <- 1 - abs(adj_matrix)
  }
  
  n <- nrow(adj_matrix)  # Number of nodes.
  
  result <- data.frame(
    Node1 = integer(),
    Node2 = integer(),
    Distance = numeric(),
    Level = integer(),
    stringsAsFactors = FALSE
  )
  
  result_cluster <- data.frame(
    Node = paste0("M(0)_", seq_len(n)),
    stringsAsFactors = FALSE
  )
  
  result_weight <- 0  
  
  # Initialize connected components; each node starts as an individual component.
  components <- seq_len(n)
  cycle <- 0
  
  # Continue until all nodes are merged into one connected network.
  while (length(unique(components)) > 1) {
    
    cycle <- cycle + 1
    result_before_cycle <- result
    result_weight_before_cycle <- result_weight
    
    min_edges_all <- data.frame(
      u = integer(),
      v = integer(),
      Distance = numeric(),
      stringsAsFactors = FALSE
    )
    
    # Find the shortest edge(s) for each connected component.
    for (comp in unique(components)) {
      comp_nodes <- which(components == comp)
      
      # Identify all candidate edges from the current component to other components.
      candidate_edges <- expand.grid(
        u = comp_nodes,
        v = setdiff(seq_len(n), comp_nodes),
        KEEP.OUT.ATTRS = FALSE
      )
      candidate_edges$Distance <- adj_matrix[cbind(candidate_edges$u, candidate_edges$v)]
      candidate_edges <- candidate_edges[candidate_edges$Distance > 0, , drop = FALSE]  # Remove invalid edges.
      
      if (nrow(candidate_edges) > 0) {
        # Select all edges that share the minimum distance.
        min_distance <- min(candidate_edges$Distance)
        min_edges <- candidate_edges[candidate_edges$Distance == min_distance, , drop = FALSE]
        min_edges_all <- rbind(min_edges_all, min_edges)
      }
    }
    
    if (nrow(min_edges_all) == 0) {
      stop("No valid candidate edges remain. The input graph may be disconnected or contain no positive distances between components.")
    }
    
    # Add all selected minimum edges to the result.
    for (i in seq_len(nrow(min_edges_all))) {
      u <- min_edges_all$u[i]
      v <- min_edges_all$v[i]
      weight <- min_edges_all$Distance[i]
      edge_u <- min(u, v)
      edge_v <- max(u, v)
      
      # If the two nodes belong to different components, merge the components and record the edge.
      if (components[u] != components[v]) {
        result <- rbind(
          result,
          data.frame(
            Node1 = edge_u,
            Node2 = edge_v,
            Distance = weight,
            Level = cycle,
            stringsAsFactors = FALSE
          )
        )
        result_weight <- result_weight + weight
        
        # Merge the two connected components.
        components[components == components[v]] <- components[u]
      } else {
        # Check whether the edge has already been recorded.
        edge_index <- which(result$Node1 == edge_u & result$Node2 == edge_v)
        
        if (length(edge_index) == 0) {
          # If the edge is not in the result, add it as a new edge.
          result <- rbind(
            result,
            data.frame(
              Node1 = edge_u,
              Node2 = edge_v,
              Distance = weight,
              Level = cycle,
              stringsAsFactors = FALSE
            )
          )
        } else {
          # If the edge already exists, do not add a duplicate row.
          next
        }
      }
    }
    
    # Record module labels and edges for intermediate levels only.
    # The final global module (#cluster = 1) and its corresponding edges are not exported.
    if (length(unique(components)) > 1) {
      freq_table <- sort(table(components), decreasing = TRUE)
      label_map <- setNames(seq_along(freq_table), names(freq_table))
      result_cluster[[paste0("Level:", cycle)]] <- paste0(
        "M(", cycle, ")_",
        as.integer(label_map[as.character(components)])
      )
    } else {
      result <- result_before_cycle
      result_weight <- result_weight_before_cycle
    }
  }
  
  if (nrow(result) > 0) {
    result$Node1 <- paste0("M(0)_", result$Node1)
    result$Node2 <- paste0("M(0)_", result$Node2)
  }
  
  # Convert distances back to similarity values if the input was not a distance matrix.
  if (!dist_matrix) {
    result$Similarity <- 1 - abs(result$Distance)
  }
  
  # Return the MSCN edge list, and module assignment table.
  return(list(edges = result, cluster = result_cluster))
}


results <- minimum_span_clustering_network(sample_matrix)

results$cluster
results$edges
