library(tibble)
library(dplyr)

# ── Helpers ───────────────────────────────────────────────────────────────────

# Binary entropy H(p) in bits; vectorized, 0 at boundaries.
h_bin <- function(p) {
  ifelse(p <= 0 | p >= 1, 0, -p * log2(p) - (1 - p) * log2(1 - p))
}

# I(X ; Y) where Y is binary, X has arbitrary support.
# prior : probability vector over X values
# p_y1  : P(Y=1 | X=x_i) for each x_i
mi_binary_y <- function(prior, p_y1) {
  hy   <- h_bin(sum(prior * p_y1))
  hy_x <- sum(prior * h_bin(p_y1))
  hy - hy_x
}

# Internal: compute MI for one p_y1 vector. Returns named numeric vector.
# cell_prior = "sum" : P(cell) = sum of element priors (proper marginalization)
# cell_prior = "flat": P(cell) = 1/K uniform over cells; P(Y=1|cell) = mean of p_y1 in cell
mi_one <- function(p_y1, partitions, prior, cell_prior = "sum") {
  mi_full <- mi_binary_y(prior, p_y1)

  mi_parts <- vapply(partitions, function(part) {
    K <- length(part)
    if (cell_prior == "sum") {
      prior_cells <- sapply(part, function(g) sum(prior[g]))
      p_y1_cells  <- sapply(seq_along(part), function(ci) {
        g  <- part[[ci]]
        pc <- prior_cells[ci]
        if (pc == 0) 0 else sum(prior[g] * p_y1[g]) / pc
      })
    } else {
      prior_cells <- rep(1 / K, K)
      p_y1_cells  <- sapply(part, function(g) mean(p_y1[g]))
    }
    mi_binary_y(prior_cells, p_y1_cells)
  }, numeric(1))

  part_labels <- if (!is.null(names(partitions))) {
    names(partitions)
  } else {
    sapply(partitions, function(part) {
      cells <- sapply(part, function(g) paste0("{", paste(g, collapse = ","), "}"))
      paste(cells, collapse = " | ")
    })
  }

  c(`I(X;Y)` = mi_full, setNames(mi_parts, part_labels))
}

# ── Main function ─────────────────────────────────────────────────────────────

# Compute MI for X (4 values) and binary Y across a list of probability
# vectors, returning a tibble with one row per vector.
#
# prob_list  : list of length-4 numeric vectors, each giving P(Y=1 | X=x_i);
#              a single vector is also accepted and wrapped automatically
# partitions : list of partitions; each partition is a list of integer vectors
#              giving the x-indices in each cell, e.g.
#              list(c(1,2), c(3,4))  groups {x1,x2} vs {x3,x4}
#              Names on `partitions` are used as column labels.
# prior      : length-4 probability vector over X; defaults to uniform
# name       : optional experiment label, stored in the `name` column
# cell_prior : method for assigning probabilities to partition cells; either
#              "sum" (default) or "flat"; see `mi_one()` for details
#
# Columns returned:
#   name     — experiment label
#   p_y1     — list column holding each input probability vector
#   I(X;Y)   — standard mutual information (bits)
#   <label>  — MI after coarsening X by each named partition (bits)
#   loss-<label> — I(X;Y) minus the corresponding partition MI (bits)
explore_mi <- function(prob_list,
                       partitions = list(),
                       prior      = rep(0.25, 4),
                       name       = NA_character_,
                       cell_prior = "sum") {

  if (is.numeric(prob_list)) prob_list <- list(prob_list)

  stopifnot(
    all(sapply(prob_list, length) == 4),
    all(sapply(prob_list, function(v) all(v >= 0) && all(v <= 1))),
    length(prior) == 4,
    all(prior >= 0),
    abs(sum(prior) - 1) < 1e-9,
    cell_prior %in% c("sum", "flat")
  )

  rows <- lapply(prob_list, function(p_y1) {
    mis      <- mi_one(p_y1, partitions, prior, cell_prior)
    mi_std   <- mis[["I(X;Y)"]]
    part_nms <- setdiff(names(mis), "I(X;Y)")
    diffs    <- setNames(mi_std - mis[part_nms], paste0("loss-", part_nms))

    as_tibble(c(
      list(name = name, p_y1 = list(p_y1)),
      as.list(mis),
      as.list(diffs)
    ))
  })

  bind_rows(rows)
}


# ── Examples ──────────────────────────────────────────────────────────────────

parts <- list(
  "compressed"  = list(c(1, 2), c(3, 4)),
  "high" = list(c(1),    c(2, 3, 4)),
  "low"  = list(c(2),    c(1, 3, 4))
)

# Experiment 1: varying P(Y=1 | X=x1) while keeping others fixed at 0.70, 0.01, 0.01.
results_exp1 <- explore_mi(
  prob_list = list(
    c(0.70, 0.70, 0.01, 0.01),
    c(0.85, 0.70, 0.01, 0.01),
    c(0.98, 0.70, 0.01, 0.01)
  ),
  partitions = parts,
  name = "experiment 1",
  cell_prior = "flat"
) |> 
  mutate(
    compVShigh = `compressed` - `high`
  )

results_exp1


results_exp2 <- explore_mi(
  prob_list = list(
    c(0.55, 0.55, 0.01, 0.01),
    c(0.85, 0.55, 0.01, 0.01),
    c(0.98, 0.55, 0.01, 0.01)
  ),
  partitions = parts,
  name = "experiment 2",
  cell_prior = "flat"
) |>
  mutate(
    compVShigh = `compressed` - `high`
  )

results_exp2
