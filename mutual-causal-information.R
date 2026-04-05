library(tidyverse)

# manual computation
U1 <- .54 * log(.54 / (.6 * .54), base = exp(1)) +
      .06 * log(.06 / (.6 * .46), base = exp(1)) +
       .4 * log(.4  / (.4 * .46), base = exp(1))

U2 <-
  .54 * log(.54 / (.9 * .54), base = exp(1)) +
  .36 * log(.36 / (.9 * .46), base = exp(1)) +
   .1 * log(.1  / (.1 * .46), base = exp(1))

# function to compute mutual information for two discrete probability distributions
# input is in long-table format (X, Y, joint probability)
mutual_information <- function(data, x_var = "X", y_var = "Y") {
  names(data)[names(data) == x_var] <- "X"
  names(data)[names(data) == y_var] <- "Y"
  # aggregate duplicate (X, Y) pairs by summing their probabilities
  data <- aggregate(prob ~ X + Y, data = data, sum)
  
  # complete missing (X, Y) combinations with prob = 0
  all_combinations <- expand.grid(
    X = unique(data$X),
    Y = unique(data$Y),
    stringsAsFactors = FALSE
  )
  data <- merge(all_combinations, data, by = c("X", "Y"), all.x = TRUE)
  data$prob[is.na(data$prob)] <- 0
  
  # compute marginal probabilities
  p_x <- aggregate(prob ~ X, data = data, sum)
  p_y <- aggregate(prob ~ Y, data = data, sum)
  
  # merge marginal probabilities with joint probabilities
  data <- merge(data, p_x, by = "X", suffixes = c("", "_x"))
  data <- merge(data, p_y, by = "Y", suffixes = c("", "_y"))
  
  # compute mutual information (0 * log(0) is treated as 0 per convention)
  data$mi_contrib <- ifelse(
    data$prob == 0,
    0,
    data$prob * log(data$prob / (data$prob_x * data$prob_y), base = exp(1))
  )
  mi <- sum(data$mi_contrib, na.rm = TRUE)
  return(mi)
}

# conjunctive case:
#    Urn 1 has 0.2/0.8 black/white balls
#    Urn 2 has 0.8/0.2 black/white balls
#    => we expect Urn 1 to have higher mutual information with the result than Urn 2
prob_table_conj <-
  tibble(
    U1 = c("B", "B", "W", "W", "B", "B", "W", "W"),
    U2 = c("B", "W", "B", "W", "B", "W", "B", "W"),
    R = c("Win", "Win", "Win", "Win", "Lose", "Lose", "Lose", "Lose"),
    prob = c(0.54, 0, 0, 0, 0, 0.06, 0.36, 0.04)
  )
mutual_information(
  prob_table_conj |> select(U1,R, prob) |> rename(X = U1, Y = R))
mutual_information(
  prob_table_conj |> select(U2, R, prob) |> rename(X = U2, Y = R)
)

# disjunctive case
#    Urn 1 has 0.2/0.8 black/white balls
#    Urn 2 has 0.8/0.2 black/white balls
#    => we expect Urn 2 to have higher mutual information with the result than Urn 1
prob_table_dis <-
  tibble(
    U1 = c("B", "B", "W", "W", "B", "B", "W", "W"),
    U2 = c("B", "W", "B", "W", "B", "W", "B", "W"),
    P1 = c("P", "q", "q", "q", "q", "q", "q", "q"),
    R = c("Win", "Win", "Win", "Win", "Lose", "Lose", "Lose", "Lose"),
    prob = c(0.54, 0.06, 0.36, 0, 0, 0, 0, 0.04)
  )
mutual_information(
  prob_table_dis |> select(U1,R, prob) |> rename(X = U1, Y = R))
mutual_information(
  prob_table_dis |> select(U2, R, prob) |> rename(X = U2, Y = R)
)
mutual_information(prob_table_dis, x_var = "P1", y_var = "R")



# playground
# Bricofly ::: causal relevance of "moderate temperature"
#  === "COMPRESSED" condition
tibble(
  X = c("moderate", "moderate", "extreme", "extreme"),
  Y = c("blue wings", "red wings", "blue wings", "red wings"),
  prob = c(0.84, 0.16, 0.01, 0.99) * 0.5
) |> 
  mutual_information()

# Bricofly ::: causal relevance of "moderate warm temperature" (under Laplace lumping)
#  === "HIGH" condition
tibble(
  X = c("moderate", "moderate", "extreme", "extreme"),
  Y = c("blue wings", "red wings", "blue wings", "red wings"),
  prob = c(0.24, 0.76, 0.98, 0.02) * 0.5
) |> 
  mutual_information()

# Bricofly ::: causal relevance of "moderate warm temperature" (under Bayes lumping)
#  === "HIGH" condition (altnerative)
tibble(
  X = c("moderate", "moderate", "extreme", "extreme"),
  Y = c("blue wings", "red wings", "blue wings", "red wings"),
  prob = c(c(0.24, 0.76) * 0.75, c(0.98, 0.02) * 0.25)
) |> 
  mutual_information()


# Bricofly ::: causal mutual information of uncompressed model (scenario 1)
blue_wing_probs <- c(0.01, 0.7, 0.7, 0.01)
tibble(
  X = rep(c("cold-extreme", "cold-moderate", "warm-moderate", "warm-extreme"),
          times = 2),
  Y = rep(c("blue wings", "red wings"), each = 4),
  prob = c(blue_wing_probs, 1 - blue_wing_probs) * 0.25
) |>
  mutual_information()


# Bricofly ::: causal mutual information of uncompressed model (scenario 2)
blue_wing_probs <- c(0.01, 0.7, 0.98, 0.01)
tibble(
  X = rep(c("cold-extreme", "cold-moderate", "warm-moderate", "warm-extreme"),
          times = 2),
  Y = rep(c("blue wings", "red wings"), each = 4),
  prob = c(blue_wing_probs, 1 - blue_wing_probs) * 0.25
) |>
  mutual_information()

# Bricofly ::: causal mutual information of uncompressed model (scenario 3)
blue_wing_probs <- c(0.43, 0.7, 0.97, 0.43)
tibble(
  X = rep(c("cold-extreme", "cold-moderate", "warm-moderate", "warm-extreme"),
          times = 2),
  Y = rep(c("blue wings", "red wings"), each = 4),
  prob = c(blue_wing_probs, 1 - blue_wing_probs) * 0.25
) |>
  mutual_information()



# Bricofly ::: causal mutual information of uncompressed model (scenario 3)

get_predictions <- function(blue_wing_probs, scenario_name){

  probability_table <- tibble(
    C = rep(c("cold-extreme", "cold-moderate", "warm-moderate", "warm-extreme"),
              times = 2),
    P1 = case_when(C == "warm-moderate" ~ 1, TRUE ~ 0),
    P2 = case_when(C %in% c("cold-moderate", "warm-moderate") ~ 1, TRUE ~ 0),B
    Y = rep(c("blue wings", "red wings"), each = 4),
    prob = c(blue_wing_probs, 1 - blue_wing_probs) * 0.25
  ) 
  tibble(
    scenario = scenario_name,
    KL_compressed = probability_table |> mutual_information(x_var = "P2"),
    KL_high = probability_table |> mutual_information(x_var = "C"),
    MF_compressed = probability_table |> mutual_information(x_var = "P2"),
    MF_high = probability_table |> mutual_information(x_var = "P1")
  )
}

results <- bind_rows(
  get_predictions(c(0.01, 0.7, 0.7, 0.01),  "Scenario 1"),
  get_predictions(c(0.01, 0.7, 0.98, 0.01), "Scenario 2"),
  get_predictions(c(0.43, 0.7, 0.7, 0.43), "Scenario 3")
)


# Sophie's pecking ::: Yablo's scenario

sophie_pecking <- tibble(
  C = rep(c("r1", "r2", "b1", "b2"), each = 2),
  E = rep(c("peck", "abstain"), times = 4),
  P1 = case_when(C %in% c("r1") ~ 1, TRUE ~ 0),
  P2 = case_when(C %in% c("r1", "r2") ~ 1, TRUE ~0),
  prob = c(0.25, 0, 0.25, 0, 0, 0.25, 0, 0.25)
)

mutual_information(sophie_pecking, x_var = "P1", y_var = "E")
mutual_information(sophie_pecking, x_var = "P2", y_var = "E")
