library(tidyverse)

# ==============================================================================
# Causal Bayes Net (CBN) — lightweight S3 class
# ==============================================================================
#
# A CBN encodes a joint probability distribution over discrete variables via:
#   exogenous  — root variables with unconditional prior distributions
#   endogenous — non-root variables with conditional probability tables (CPTs)
#
# The DAG topology is inferred from CPT column names:
#   every column in a CPT other than the variable's own name and "prob" is a parent.
#
# Example:
#
#   net <- CBN(
#     exogenous = list(
#       Rain    = c(yes = 0.3, no = 0.7),
#       Sprink  = c(yes = 0.4, no = 0.6)
#     ),
#     endogenous = list(
#       Wet = tribble(
#         ~Rain,  ~Sprink,  ~Wet,   ~prob,
#         "yes",  "yes",    "yes",   1.00,
#         "yes",  "yes",    "no",    0.00,
#         "yes",  "no",     "yes",   0.90,
#         "yes",  "no",     "no",    0.10,
#         "no",   "yes",    "yes",   0.80,
#         "no",   "yes",    "no",    0.20,
#         "no",   "no",     "yes",   0.10,
#         "no",   "no",     "no",    0.90
#       )
#     )
#   )
# ==============================================================================

CBN <- function(exogenous = list(), endogenous = list()) {

  # --- validate exogenous priors -----------------------------------------------
  for (var in names(exogenous)) {
    p <- exogenous[[var]]
    if (!is.numeric(p) || is.null(names(p)))
      stop(sprintf("Exogenous '%s': must be a named numeric vector", var))
    if (any(p < 0))
      stop(sprintf("Exogenous '%s': probabilities must be non-negative", var))
    exogenous[[var]] <- p / sum(p)
  }

  all_vars <- c(names(exogenous), names(endogenous))

  # --- validate CPTs and infer parents -----------------------------------------
  parents <- setNames(vector("list", length(endogenous)), names(endogenous))

  for (var in names(endogenous)) {
    cpt <- endogenous[[var]]

    if (!is.data.frame(cpt))
      stop(sprintf("CPT for '%s' must be a data frame", var))
    if (!var %in% names(cpt))
      stop(sprintf("CPT for '%s' must contain a column named '%s'", var, var))
    if (!"prob" %in% names(cpt))
      stop(sprintf("CPT for '%s' must contain a column named 'prob'", var))

    pa <- setdiff(names(cpt), c(var, "prob"))
    parents[[var]] <- pa

    unknown <- setdiff(pa, all_vars)
    if (length(unknown) > 0)
      stop(sprintf("CPT for '%s' references unknown variables: %s",
                   var, paste(unknown, collapse = ", ")))

    # check non-negativity and normalize within each conditioning context
    if (any(cpt$prob < 0))
      stop(sprintf("CPT for '%s': probabilities must be non-negative", var))
    if (length(pa) == 0) {
      endogenous[[var]]$prob <- cpt$prob / sum(cpt$prob)
    } else {
      contexts <- unique(cpt[, pa, drop = FALSE])
      for (i in seq_len(nrow(contexts))) {
        mask <- Reduce(`&`, lapply(pa, \(col) cpt[[col]] == contexts[[col]][i]))
        endogenous[[var]]$prob[mask] <- cpt$prob[mask] / sum(cpt$prob[mask])
      }
    }
  }

  # --- check acyclicity --------------------------------------------------------
  topo <- .cbn_topo_sort(names(exogenous), names(endogenous), parents)

  structure(
    list(
      exogenous  = exogenous,
      endogenous = endogenous,
      parents    = parents,
      topo_order = topo
    ),
    class = "cbn"
  )
}

# Kahn's algorithm — returns topological order or errors on cycle
.cbn_topo_sort <- function(exo_vars, endo_vars, parents) {
  all_vars  <- c(exo_vars, endo_vars)
  in_degree <- setNames(integer(length(all_vars)), all_vars)
  children  <- setNames(vector("list",  length(all_vars)), all_vars)
  for (v in all_vars) children[[v]] <- character(0)

  for (var in endo_vars) {
    for (pa in parents[[var]]) {
      in_degree[var]  <- in_degree[var] + 1L
      children[[pa]]  <- c(children[[pa]], var)
    }
  }

  queue <- names(in_degree)[in_degree == 0]
  order <- character(0)

  while (length(queue) > 0) {
    v     <- queue[1]
    queue <- queue[-1]
    order <- c(order, v)
    for (ch in children[[v]]) {
      in_degree[ch] <- in_degree[ch] - 1L
      if (in_degree[ch] == 0) queue <- c(queue, ch)
    }
  }

  if (length(order) < length(all_vars))
    stop("CBN contains a cycle — the graph must be acyclic")

  order
}

print.cbn <- function(x, ...) {
  has_partition <- "partition" %in% names(x$endogenous)
  has_compression <- "compression" %in% names(x$intervened)
  ctx_driven      <- if (has_compression) x$intervened[["compression"]]$var_names else character(0)

  # partition and compression-driven variables are displayed under "Interventions"
  endo_display <- setdiff(names(x$endogenous), c("partition", ctx_driven))

  n_exo  <- length(x$exogenous)
  n_endo <- length(endo_display)
  n_intv <- length(setdiff(names(x$intervened), "compression"))  # compression shown separately
  counts <- paste(Filter(nchar, c(
    if (n_exo  > 0) sprintf("%d exogenous",  n_exo),
    if (n_endo > 0) sprintf("%d endogenous", n_endo),
    if (n_intv > 0) sprintf("%d intervened", n_intv)
  )), collapse = ", ")
  cat(sprintf("Causal Bayes Net  [%s]\n", counts))

  # edges: exclude partition (internal); compression-driven vars ARE shown (compression->var->...)
  edge_vars <- setdiff(names(x$parents), "partition")
  edges <- Filter(nchar, vapply(edge_vars, function(var) {
    pa <- x$parents[[var]]
    if (length(pa) == 0) "" else sprintf("%s -> %s", paste(pa, collapse = ", "), var)
  }, character(1)))
  if (length(edges) > 0)
    cat(sprintf("  Edges: %s\n", paste(edges, collapse = ";  ")))
  cat("\n")

  if (n_exo > 0) {
    cat("Exogenous priors:\n")
    for (var in names(x$exogenous)) {
      p    <- x$exogenous[[var]]
      vals <- paste(names(p), sprintf("%.3g", p), sep = " = ", collapse = ",  ")
      cat(sprintf("  %-12s  %s\n", var, vals))
    }
    cat("\n")
  }

  if (n_endo > 0) {
    cat("Endogenous variables:\n")
    for (var in endo_display) {
      pa     <- x$parents[[var]]
      pa_str <- if (length(pa) == 0) "(none)" else paste(pa, collapse = ", ")
      cat(sprintf("  %-12s  parents: %s\n", var, pa_str))
    }
  }

  regular_intv <- setdiff(names(x$intervened), "compression")
  if (length(regular_intv) > 0 || has_compression || has_partition) {
    cat("\nInterventions  [do(...)]:\n")

    for (var in regular_intv) {
      p    <- x$intervened[[var]]$prior
      vals <- paste(names(p), sprintf("%.3g", p), sep = " = ", collapse = ",  ")
      cat(sprintf("  %-12s  %s\n", var, vals))
    }

    if (has_compression) {
      ctx <- x$intervened[["compression"]]
      focal_str <- paste(
        mapply(\(v, f) sprintf("%s \u2208 {%s}", v, paste(f, collapse = ", ")),
               ctx$var_names, ctx$focal_groups),
        collapse = ";  "
      )
      cat(sprintf("  Set-level intervention:  %s\n", focal_str))
      cat("  Context prior:\n")
      p <- ctx$prior
      for (i in seq_along(p))
        cat(sprintf("    %-14s  prob = %.3g\n", names(p)[i], p[i]))
    }

    if (has_partition) {
      cat("  Partition:\n")
      cell_probs <- joint(x, vars = "partition")
      for (i in seq_len(nrow(cell_probs)))
        cat(sprintf("    cell %-4s  prob = %.3g\n",
                    cell_probs$partition[i], cell_probs$prob[i]))
    }
  }

  invisible(x)
}

# ------------------------------------------------------------------------------
# joint(x, vars = NULL, cond = NULL)
#
# Computes the (conditional) marginal joint probability table for a CBN.
#   vars — character vector of variable names to include; NULL returns the full
#          joint over all variables.
#   cond — named list or named character vector specifying conditioning values.
#          Each entry gives the allowed values for one variable, e.g.:
#            cond = list(X = c("x1","x2"))   — P(... | X ∈ {x1, x2})
#            cond = c(X = "x1")              — P(... | X = x1)
#          Rows not matching the condition are dropped and probs renormalised.
#
# The full joint is built by multiplying exogenous priors and endogenous CPTs
# in topological order, conditioning (if cond supplied), then marginalizing.
# ------------------------------------------------------------------------------

joint <- function(x, ...) UseMethod("joint")

joint.cbn <- function(x, vars = NULL, cond = NULL) {

  # --- normalise cond to a named list -----------------------------------------
  # accepts: named list (multiple values per var) or named character vector
  if (!is.null(cond)) {
    if (!is.list(cond)) {
      if (is.null(names(cond)) || any(!nzchar(names(cond))))
        stop("'cond' must be a named list or named character vector")
      cond <- as.list(cond)
    }
    if (is.null(names(cond)) || any(!nzchar(names(cond))))
      stop("'cond' must be a named list with variable names as keys")
  }
  cond_vars <- names(cond)

  # --- validate vars and cond_vars --------------------------------------------
  all_vars <- c(names(x$exogenous), names(x$intervened), names(x$endogenous))
  if (!is.null(vars)) {
    unknown <- setdiff(vars, all_vars)
    if (length(unknown) > 0)
      stop(sprintf("Unknown variables: %s", paste(unknown, collapse = ", ")))
  }
  if (length(cond_vars) > 0) {
    unknown <- setdiff(cond_vars, all_vars)
    if (length(unknown) > 0)
      stop(sprintf("Unknown conditioning variables: %s", paste(unknown, collapse = ", ")))
  }

  # --- collect all possible values for each variable (in topo order) ----------
  all_values <- setNames(
    lapply(x$topo_order, \(var) {
      if (var %in% names(x$exogenous))   names(x$exogenous[[var]])
      else if (var %in% names(x$intervened)) names(x$intervened[[var]]$prior)
      else unique(x$endogenous[[var]][[var]])
    }),
    x$topo_order
  )

  # --- build full combination grid --------------------------------------------
  tbl <- do.call(expand.grid, c(all_values, list(stringsAsFactors = FALSE))) |>
    as_tibble() |>
    mutate(prob = 1.0)

  # --- multiply in exogenous priors -------------------------------------------
  for (var in names(x$exogenous)) {
    prior    <- x$exogenous[[var]]
    tbl$prob <- tbl$prob * prior[tbl[[var]]]
  }

  # --- multiply in intervened (degenerate) priors -----------------------------
  for (var in names(x$intervened)) {
    prior    <- x$intervened[[var]]$prior
    tbl$prob <- tbl$prob * prior[tbl[[var]]]
  }

  # --- multiply in endogenous CPTs (topological order) ------------------------
  for (var in x$topo_order[x$topo_order %in% names(x$endogenous)]) {
    cpt       <- x$endogenous[[var]]
    join_cols <- c(x$parents[[var]], var)
    tbl <- tbl |>
      left_join(
        cpt |> select(all_of(join_cols), prob) |> rename(prob_cpt = prob),
        by = join_cols
      ) |>
      mutate(prob = prob * prob_cpt) |>
      select(-prob_cpt)
  }

  # --- condition: filter to matching rows and renormalise ---------------------
  if (length(cond_vars) > 0) {
    mask <- Reduce(`&`, lapply(cond_vars, \(v) tbl[[v]] %in% cond[[v]]))
    tbl  <- tbl[mask, ]
    tbl$prob <- tbl$prob / sum(tbl$prob)
  }

  # --- marginalize if vars specified ------------------------------------------
  keep_vars <- if (!is.null(vars)) vars else setdiff(names(tbl), "prob")
  tbl <- tbl |>
    group_by(across(all_of(keep_vars))) |>
    summarise(prob = sum(prob), .groups = "drop")

  arrange(tbl, across(-prob))
}

# ------------------------------------------------------------------------------
# MI(x, x_var, y_var)
#
# Computes the mutual information I(X; Y) between two variables in a CBN,
# using their marginal joint probability table.
#
#   MI(net, "X", "Y")
#
# Returns a single numeric value. Logarithm base is controlled by 'base'
# (default exp(1) = nats; use base = 2 for bits).
# Cells with prob = 0 contribute 0 to the sum (0 log 0 := 0 by convention).
# ------------------------------------------------------------------------------

MI <- function(x, ...) UseMethod("MI")

MI.cbn <- function(x, x_var, y_var, base = exp(1), ...) {
  tbl <- joint(x, vars = c(x_var, y_var))
  names(tbl)[names(tbl) == x_var] <- ".X"
  names(tbl)[names(tbl) == y_var] <- ".Y"

  p_x <- aggregate(prob ~ .X, data = tbl, sum)
  p_y <- aggregate(prob ~ .Y, data = tbl, sum)
  tbl <- merge(tbl, p_x, by = ".X", suffixes = c("", "_x"))
  tbl <- merge(tbl, p_y, by = ".Y", suffixes = c("", "_y"))

  sum(ifelse(tbl$prob == 0, 0,
             tbl$prob * log(tbl$prob / (tbl$prob_x * tbl$prob_y), base = base)))
}

# ------------------------------------------------------------------------------
# cMI(x, intervention, y_var, prior = "marginal", base = exp(1))
#
# Causal mutual information: applies a do-intervention and computes MI between
# the intervened cause and an effect variable.
#
# The cause variable used for MI is determined automatically:
#   - List intervention (use case 3)          → "partition"
#   - Named vector, single var (use case 1)   → that variable
#   - Unnamed vector, single var (use case 2) → that variable
#   For multi-variable use cases 1/2, pass a list to obtain a partition instead.
#
# Examples:
#   cMI(net, list(X = c("x1","x2")), "Y")            # partition vs Y
#   cMI(net, c("X"), "Y", prior = "flat")             # X vs Y after soft do(X)
#   cMI(net, c(X = "x1"), "Y")                        # X vs Y after hard do(X=x1)
# ------------------------------------------------------------------------------

cMI <- function(x, intervention, y_var, prior = "marginal", base = exp(1)) {

  # --- determine cause variable ------------------------------------------------
  if (is.list(intervention)) {
    x_var <- "partition"
  } else {
    var_names <- if (!is.null(names(intervention))) names(intervention) else intervention
    if (length(var_names) > 1)
      stop("cMI: for multi-variable interventions use a list (use case 3) to obtain a partition")
    x_var <- var_names
  }

  intervened_net <- intervene(x, intervention, prior = prior)
  MI(intervened_net, x_var, y_var, base = base)
}

# ------------------------------------------------------------------------------
# intervene(x, intervention, prior = "marginal")
#
# Applies do-interventions to a CBN, returning a new (mutilated) CBN.
# Three use cases, distinguished by the type of 'intervention':
#
#   Use case 1 — hard (degenerate) intervention; named character vector:
#     intervene(net, c("X" = "x2", "Z" = "z37"))
#     Each variable is pinned to its stated value (degenerate prior).
#     'prior' must NOT be supplied; a warning is issued if it is.
#
#   Use case 2 — soft intervention; unnamed character vector:
#     intervene(net, c("X", "Z"))
#     intervene(net, c("X", "Z"), prior = "flat")
#     The 'prior' argument sets the distribution over each intervened variable:
#       "marginal" (default) — P(var) from the original (pre-intervention) model
#       "flat"               — maximum-entropy (uniform) distribution
#
#   Use case 3 — set-level / partitioned intervention; named list of value groups:
#     intervene(net, list(X = c("x1","x2"), Z = c("z3","z4")))
#     intervene(net, list(X = c("x1","x2"), Z = c("z3","z4")), prior = "flat")
#     Each entry specifies the "focal group" of values for that variable.
#     Internally: (a) performs a use-case-2 soft intervention on the listed
#     variables, then (b) adds a deterministic 'partition' variable whose CPT
#     maps each combination of variable values to a cell:
#       "1" — ALL intervened variables are in their focal group
#       "0" — otherwise
#     The 'prior' argument has the same meaning as in use case 2, applied to the
#     individual variable values (not directly to partition cells).
#     The 'partition' variable is included in joint() output for inspection.
#
# The original model is stored as attr(result, "original_model").
# ------------------------------------------------------------------------------

intervene <- function(x, ...) UseMethod("intervene")

intervene.cbn <- function(x, intervention, prior = "marginal") {

  # use case 3: list of focal value groups
  if (is.list(intervention))
    return(.intervene_partition(x, intervention, prior))

  # --- detect use case ---------------------------------------------------------
  has_values <- !is.null(names(intervention))

  if (has_values && !missing(prior))
    warning("'prior' is ignored when intervention values are specified (use case 1)")

  if (!has_values && !prior %in% c("marginal", "flat"))
    stop("'prior' must be \"marginal\" or \"flat\"")

  var_names <- if (has_values) names(intervention) else intervention
  all_vars  <- c(names(x$exogenous), names(x$endogenous))

  # --- validate variable names -------------------------------------------------
  unknown <- setdiff(var_names, all_vars)
  if (length(unknown) > 0)
    stop(sprintf("Unknown variables: %s", paste(unknown, collapse = ", ")))

  # --- validate values (use case 1 only) ---------------------------------------
  if (has_values) {
    for (var in var_names) {
      val        <- intervention[[var]]
      valid_vals <- if (var %in% names(x$exogenous))
                      names(x$exogenous[[var]])
                    else
                      unique(x$endogenous[[var]][[var]])
      if (!val %in% valid_vals)
        stop(sprintf("Variable '%s': value '%s' not in domain {%s}",
                     var, val, paste(valid_vals, collapse = ", ")))
    }
  }

  # --- helper: compute the prior for one variable ------------------------------
  .make_prior <- function(var) {
    valid_vals <- if (var %in% names(x$exogenous))
                    names(x$exogenous[[var]])
                  else
                    unique(x$endogenous[[var]][[var]])
    if (has_values) {
      # degenerate at the specified value
      p      <- setNames(rep(0, length(valid_vals)), valid_vals)
      p[intervention[[var]]] <- 1
      p
    } else if (prior == "flat") {
      setNames(rep(1 / length(valid_vals), length(valid_vals)), valid_vals)
    } else {  # "marginal"
      marg <- joint(x, vars = var)
      setNames(marg$prob, marg[[var]])
    }
  }

  # --- build mutilated network -------------------------------------------------
  new_exogenous  <- x$exogenous
  new_endogenous <- x$endogenous
  new_intervened <- if (is.null(x$intervened)) list() else x$intervened

  for (var in var_names) {
    p   <- .make_prior(var)
    val <- if (has_values) intervention[[var]] else NA_character_

    if (var %in% names(x$endogenous)) {
      new_intervened[[var]] <- list(
        cpt   = x$endogenous[[var]],
        prior = p,
        value = val
      )
      new_endogenous[[var]] <- NULL
    } else {
      new_intervened[[var]] <- list(
        original_prior = x$exogenous[[var]],
        prior          = p,
        value          = val
      )
      new_exogenous[[var]] <- NULL
    }
  }

  # recompute parents (only for remaining endogenous variables)
  new_parents <- x$parents[names(new_endogenous)]

  # recompute topo order: intervened variables are root nodes like exogenous ones
  new_topo <- .cbn_topo_sort(
    c(names(new_exogenous), names(new_intervened)),
    names(new_endogenous),
    new_parents
  )

  result <- structure(
    list(
      exogenous  = new_exogenous,
      endogenous = new_endogenous,
      intervened = new_intervened,
      parents    = new_parents,
      topo_order = new_topo
    ),
    class = "cbn"
  )

  attr(result, "original_model") <- x
  result
}

# Internal helper for use case 3 — set-level / partitioned intervention
#
# Architecture: introduces a root 'compression' variable whose values enumerate
# every atomic (X, Y, ...) combination. The intervened variables and 'partition'
# become deterministic endogenous children of 'compression'. This gives exact
# flat-over-cells priors for any partition shape, including non-rectangular ones,
# because the prior is placed directly on combinations rather than on independent
# variable marginals.
.intervene_partition <- function(x, focal_groups, prior) {

  var_names <- names(focal_groups)
  all_vars  <- c(names(x$exogenous), names(x$endogenous))

  if (is.null(var_names) || any(!nzchar(var_names)))
    stop("Use case 3: list entries must all be named (names are variable names)")

  # --- validate variable names -------------------------------------------------
  unknown <- setdiff(var_names, all_vars)
  if (length(unknown) > 0)
    stop(sprintf("Unknown variables: %s", paste(unknown, collapse = ", ")))

  # --- validate focal values ---------------------------------------------------
  for (var in var_names) {
    valid_vals <- if (var %in% names(x$exogenous))
                    names(x$exogenous[[var]])
                  else
                    unique(x$endogenous[[var]][[var]])
    bad <- setdiff(focal_groups[[var]], valid_vals)
    if (length(bad) > 0)
      stop(sprintf("Variable '%s': unknown focal values {%s}, domain is {%s}",
                   var, paste(bad, collapse = ", "), paste(valid_vals, collapse = ", ")))
  }

  # --- name conflict check -----------------------------------------------------
  reserved <- c("compression", "partition")
  conflicts <- intersect(reserved, all_vars)
  if (length(conflicts) > 0)
    stop(sprintf("Network already has variable(s) named {%s}; cannot build compression/partition",
                 paste(conflicts, collapse = ", ")))

  # --- get domains -------------------------------------------------------------
  var_vals <- lapply(var_names, function(var) {
    if (var %in% names(x$exogenous)) names(x$exogenous[[var]])
    else unique(x$endogenous[[var]][[var]])
  })
  names(var_vals) <- var_names

  # --- enumerate atomic combinations (one compression value per combo) ---------
  combos   <- do.call(expand.grid, c(var_vals, list(stringsAsFactors = FALSE)))
  ctx_vals <- if (length(var_names) == 1)
                combos[[1]]
              else
                apply(combos[, var_names, drop = FALSE], 1, paste, collapse = ":")

  # --- determine partition cell for each combo ---------------------------------
  in_focal <- sapply(var_names, \(v) combos[[v]] %in% focal_groups[[v]])
  if (!is.matrix(in_focal)) in_focal <- matrix(in_focal, ncol = 1)
  cell <- ifelse(rowSums(!in_focal) == 0L, "1", "0")

  # --- compute compression prior -----------------------------------------------
  if (prior == "flat") {
    # exact: P(compression = c) = 1/n_cells / |cell of c|
    n_cells    <- length(unique(cell))
    cell_sizes <- table(cell)
    ctx_prior  <- setNames((1 / n_cells) / as.numeric(cell_sizes[cell]), ctx_vals)
  } else {  # "marginal": P(compression = c) = P(X=xi, Y=yj, ...) in original model
    marg      <- joint(x, vars = var_names)
    combo_key <- if (length(var_names) == 1)
                   marg[[var_names]]
                 else
                   apply(marg[, var_names, drop = FALSE], 1, paste, collapse = ":")
    ctx_prior <- setNames(marg$prob[match(ctx_vals, combo_key)], ctx_vals)
  }

  # --- build mutilated network -------------------------------------------------
  new_exogenous  <- x$exogenous
  new_endogenous <- x$endogenous
  new_intervened <- if (is.null(x$intervened)) list() else x$intervened
  new_parents    <- x$parents

  # remove intervened vars from their original home
  for (var in var_names) {
    new_exogenous[[var]]  <- NULL
    new_endogenous[[var]] <- NULL
    new_parents[[var]]    <- NULL
  }

  # compression: intervened root node
  new_intervened[["compression"]] <- list(
    prior        = ctx_prior,
    value        = NA_character_,
    var_names    = var_names,
    focal_groups = focal_groups
  )

  # intervened variables: deterministic children of compression
  for (var in var_names) {
    vals        <- var_vals[[var]]
    var_of_ctx  <- combos[[var]]          # X-component of each compression value
    new_endogenous[[var]] <- bind_rows(lapply(seq_along(ctx_vals), \(i) {
      tibble(compression = ctx_vals[i], !!var := vals,
             prob = as.numeric(vals == var_of_ctx[i]))
    }))
    new_parents[[var]] <- "compression"
  }

  # partition: deterministic child of compression
  new_endogenous[["partition"]] <- bind_rows(lapply(seq_along(ctx_vals), \(i) {
    tibble(compression = ctx_vals[i], partition = c("0", "1"),
           prob = as.numeric(c("0", "1") == cell[i]))
  }))
  new_parents[["partition"]] <- "compression"

  # --- recompute topo order ----------------------------------------------------
  new_topo <- .cbn_topo_sort(
    c(names(new_exogenous), names(new_intervened)),
    names(new_endogenous),
    new_parents
  )

  result <- structure(
    list(
      exogenous  = new_exogenous,
      endogenous = new_endogenous,
      intervened = new_intervened,
      parents    = new_parents,
      topo_order = new_topo
    ),
    class = "cbn"
  )

  attr(result, "original_model") <- x
  result
}


# ------------------------------------------------------------------------------
# changePrior(x, vars)
#
# Returns a new CBN with updated priors for exogenous or intervened variables.
# Two calling modes:
#
#   Standard — supply full distributions:
#     changePrior(net, list(Rain = c(yes = 0.9, no = 0.1)))
#
#   Propensity — mix a degenerate prior on a focal value with the existing prior:
#     changePrior(net, list(Rain = "yes", Sprink = "yes"), propensity = 0.6)
#     New prior = propensity * delta(focal) + (1 - propensity) * old_prior
#
# Inputs are validated and normalized the same way as in CBN().
# Only exogenous and intervened variables may have their priors changed.
# ------------------------------------------------------------------------------

changePrior <- function(x, ...) UseMethod("changePrior")

changePrior.cbn <- function(x, vars, propensity = NULL) {

  valid_root_vars <- c(names(x$exogenous), names(x$intervened))

  # --- validate variable names -------------------------------------------------
  unknown <- setdiff(names(vars), valid_root_vars)
  if (length(unknown) > 0)
    stop(sprintf(
      "Variable(s) not found among exogenous or intervened: %s",
      paste(unknown, collapse = ", ")
    ))

  # --- propensity mode: vars entries are focal values (single strings) ---------
  if (!is.null(propensity)) {
    if (!is.numeric(propensity) || length(propensity) != 1 ||
        propensity < 0 || propensity > 1)
      stop("'propensity' must be a single number in [0, 1]")

    for (var in names(vars)) {
      focal <- vars[[var]]
      if (!is.character(focal) || length(focal) != 1)
        stop(sprintf("'%s': in propensity mode each entry must be a single focal value (character)", var))
      old_prior <- if (var %in% names(x$exogenous))
                     x$exogenous[[var]]
                   else
                     x$intervened[[var]]$prior
      if (!focal %in% names(old_prior))
        stop(sprintf("'%s': focal value '%s' not in domain {%s}",
                     var, focal, paste(names(old_prior), collapse = ", ")))
      degenerate     <- setNames(rep(0, length(old_prior)), names(old_prior))
      degenerate[focal] <- 1
      vars[[var]] <- propensity * degenerate + (1 - propensity) * old_prior
    }
  }

  # --- standard mode: validate and normalize full distributions ----------------
  for (var in names(vars)) {
    p <- vars[[var]]
    if (!is.numeric(p) || is.null(names(p)))
      stop(sprintf("'%s': must be a named numeric vector", var))
    if (any(p < 0))
      stop(sprintf("'%s': probabilities must be non-negative", var))
    valid_vals <- if (var %in% names(x$exogenous))
                    names(x$exogenous[[var]])
                  else
                    names(x$intervened[[var]]$prior)
    unknown_vals <- setdiff(names(p), valid_vals)
    if (length(unknown_vals) > 0)
      stop(sprintf("'%s': unknown values {%s}, domain is {%s}",
                   var, paste(unknown_vals, collapse = ", "),
                        paste(valid_vals,   collapse = ", ")))
    vars[[var]] <- p / sum(p)
  }

  # --- apply updates -----------------------------------------------------------
  result <- x
  for (var in names(vars)) {
    if (var %in% names(result$exogenous))
      result$exogenous[[var]] <- vars[[var]]
    else
      result$intervened[[var]]$prior <- vars[[var]]
  }

  attr(result, "original_model") <- x
  result
}

# ==============================================================================
# Speaker policy: soft-max on causal MIs
# ==============================================================================

# speaker policy
speaker <- function(cMIs, alpha = 1) {
  softmax <- exp(alpha * (cMIs - max(cMIs)))
  softmax <- softmax / sum(softmax)
  return(softmax)
}

# ==============================================================================
# Example: wet-grass network  (Rain, Sprinkler -> Wet)
# ==============================================================================

net <- CBN(
  exogenous = list(
    Rain   = c(yes = 0.3, no = 0.7),
    Sprink = c(yes = 0.4, no = 0.6)
  ),
  endogenous = list(
    Wet = tribble(
      ~Rain,  ~Sprink,  ~Wet,   ~prob,
      "yes",  "yes",    "yes",   1.00,
      "yes",  "yes",    "no",    0.00,
      "yes",  "no",     "yes",   0.90,
      "yes",  "no",     "no",    0.10,
      "no",   "yes",    "yes",   0.80,
      "no",   "yes",    "no",    0.20,
      "no",   "no",     "yes",   0.10,
      "no",   "no",     "no",    0.90
    )
  )
)

# # --- assessing joint probabilities (for subsets of variables) ----------------
# joint(net)                        # full joint over Rain, Sprink, Wet
# joint(net, vars = c("Wet", "Rain"))  # marginal over Rain and Wet
# joint(net, vars = "Wet")          # marginal over Wet only


# ==============================================================================
# Example: simple causal chain
# ==============================================================================

net_chain <- CBN(
  exogenous = list(
    "A" = c("1" = 0.7, "0" = 0.3)
  ),
  endogenous = list(
    "B" = tribble(
      ~A  , ~B  , ~prob ,
      "1" , "1" , 0.8   ,
      "1" , "0" , 0.2   ,
      "0" , "1" , 0.2   ,
      "0" , "0" , 0.8
    ),
    "C" = tribble(
      ~B  , ~C  , ~prob ,
      "1" , "1" , 0.6   ,
      "1" , "0" , 0.4   ,
      "0" , "1" , 0.1   ,
      "0" , "0" , 0.9
    )
  )
)

# print(net_chain)
# joint(net_chain)
# joint(net_chain, vars = c("B"))
# joint(net_chain, vars = c("C"))
# joint(net_chain, vars = c("A", "B"))

# net_chain_doB <- intervene(net_chain, c("B"), prior = "marginal")

# intervene(net_chain, c("B" = "1")) |> joint("C")
# intervene(net_chain, c("B" = "0")) |> joint("C")
# intervene(net_chain, c("B" = "1")) |> joint(c("B", "C"))

# changePrior(net_chain_doB, list("B" = "1"), propensity = .5)
# print(net_chain_doB)
# joint(net_chain_doB)
# net_chain_doB_priorTweak <- changePrior(net_chain_doB, vars = list(
#   B   = c("1" = 0.9, "0" = 0.1)
# ))
# joint(net_chain_doB_priorTweak, vars = "B")


# ==============================================================================
# Example: Sophie's pecking
# ==============================================================================

net_peck <- CBN(
  exogenous = list(
    "C" = c("R1" = 0.25, "R2" = 0.25, "B1" = 0.25, "B2" = 0.25)
  ),
  endogenous = list(
    "P" = tribble(
      ~C   ,  ~P   , ~prob,
      "R1" , "yes" ,  1   ,
      "R2" , "yes" ,  1   ,
      "B1" , "yes" ,  0   ,
      "B2" , "yes" ,  0   ,
      "R1" , "no"  ,  0   ,
      "R2" , "no"  ,  0   ,
      "B1" , "no"  ,  1   ,
      "B2" , "no"  ,  1
    )
  )
)

# net_peck_intervened <- intervene(
#   net_peck, 
#   list("C" = c("R1", "R2")), 
#   prior = "flat"
# )

# joint(net_peck)

# joint(net_peck_intervened)

# joint(net_peck_intervened, "P", cond = list(partition = 0))

# MI(net_peck_intervened, "partition", "P")

# cMIs_peck <- c(
#   "red"     = cMI(net_peck, list("C" = c("R1", "R2")), "P", prior = "flat"),
#   "scarlet" = cMI(net_peck, list("C" = c("R1")), "P", prior = "flat")
# )

# speaker(cMIs_peck, alpha=5)

# ==============================================================================
# Example: Conjunctive urns
# ==============================================================================

net_conj <- CBN(
  exogenous = list(
    "U1" = c("red" = 0.8, "blue" = 0.2),
    "U2" = c("red" = 0.2, "blue" = 0.8)
  ),
  endogenous = list(
    "win" = tribble(
      ~U1    , ~U2    , ~win , ~prob,
      "red"  , "red"  , 1    , 1 ,
      "red"  , "blue" , 1    , 0 ,
      "blue" , "red"  , 1    , 0 ,
      "blue" , "blue" , 1    , 0 ,
      "red"  , "red"  , 0    , 0 ,
      "red"  , "blue" , 0    , 1 ,
      "blue" , "red"  , 0    , 1 ,
      "blue" , "blue" , 0    , 1 
    )
  )
)

# update priors based on true world state
net_conj_updated <- changePrior(
  net_conj, 
  list(    U1 = "red", U2 = "red"),
  propensity = 0.5
)

cMIS <- c(
  "U1-red" = cMI(net_conj_updated, list("U1" = "red"), "win", prior = "marginal"),
  "U2-red" = cMI(net_conj_updated, list("U2" = "red"), "win", prior = "marginal"),
  "both"   = cMI(net_conj_updated, list("U1" = "red","U2" = "red"), "win", prior = "marginal")
)

# speaker(cMIS, alpha = 5)


# ==============================================================================
# Example: Disjunctive urns
# ==============================================================================

net_disj<- CBN(
  exogenous = list(
    "U1" = c("red" = 0.8, "blue" = 0.2),
    "U2" = c("red" = 0.2, "blue" = 0.8)
  ),
  endogenous = list(
    "win" = tribble(
      ~U1    , ~U2    , ~win , ~prob ,
      "red"  , "red"  ,    1 ,     1 ,
      "red"  , "blue" ,    1 ,     1 ,
      "blue" , "red"  ,    1 ,     1 ,
      "blue" , "blue" ,    1 ,     0 ,
      "red"  , "red"  ,    0 ,     0 ,
      "red"  , "blue" ,    0 ,     0 ,
      "blue" , "red"  ,    0 ,     0 ,
      "blue" , "blue" ,    0 ,     1
    )
  )
)

# update priors based on true world state
net_disj_updated <- changePrior(
  net_disj,
  list(U1 = "red", U2 = "red"),
  propensity = 0.5
)

cMIS <- c(
  "U1-red" = cMI(
    net_disj_updated,
    list("U1" = "red"),
    "win",
    prior = "marginal"
  ),
  "U2-red" = cMI(
    net_disj_updated,
    list("U2" = "red"),
    "win",
    prior = "marginal"
  ),
  "both" = cMI(
    net_disj_updated,
    list("U1" = "red", "U2" = "red"),
    "win",
    prior = "marginal"
  )
)

# speaker(cMIS, alpha = 5)
