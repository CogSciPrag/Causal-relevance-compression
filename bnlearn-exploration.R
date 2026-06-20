library(bnlearn)
library(gRain)

# Define the DAG structure: X1 -> Y, X2 -> Y
dag <- model2network("[X1][X2][Y|X1:X2]")

# Define conditional probability tables (CPTs)

# Prior for X1: P(X1 = TRUE) = 0.8
cpt_X1 <- matrix(c(0.2, 0.8), ncol = 2, dimnames = list(NULL, c("FALSE", "TRUE")))

# Prior for X2: P(X2 = TRUE) = 0.6
cpt_X2 <- matrix(c(0.4, 0.6), ncol = 2, dimnames = list(NULL, c("FALSE", "TRUE")))

# Conditional probability table for Y | X1, X2
# Conjunctive: Y = TRUE only when both X1 and X2 are TRUE
cpt_Y <- array(
  c(
    1.0, 0.0,  # X1=FALSE, X2=FALSE -> Y=FALSE with prob 1
    1.0, 0.0,  # X1=FALSE, X2=TRUE  -> Y=FALSE with prob 1
    1.0, 0.0,  # X1=TRUE,  X2=FALSE -> Y=FALSE with prob 1
    0.0, 1.0   # X1=TRUE,  X2=TRUE  -> Y=TRUE  with prob 1
  ),
  dim = c(2, 2, 2),
  dimnames = list(
    Y  = c("FALSE", "TRUE"),
    X1 = c("FALSE", "TRUE"),
    X2 = c("FALSE", "TRUE")
  )
)

# Fit the network with these CPTs
fitted <- custom.fit(dag, dist = list(X1 = cpt_X1, X2 = cpt_X2, Y = cpt_Y))

# Inspect
fitted
cpquery(fitted, event = (Y == "TRUE"), evidence = (X1 == "TRUE"))
cpquery(fitted, event = (Y == "TRUE"), evidence = (X2 == "TRUE"))
cpquery(fitted, event = (Y == "TRUE"), evidence = (X1 == "TRUE" & X2 == "FALSE"))


grain_net <- as.grain(fitted)
grain_net <- compile(grain_net)

querygrain(setEvidence(grain_net, nodes = c("X1", "X2"), states = c("TRUE", "FALSE")), nodes = "Y")
querygrain(grain_net, nodes = "Y")


# General approach via graph mutilation (do-operator)
dag_mutilated <- mutilated(fitted, evidence = list(X1 = "TRUE"))
grain_do_x1_true <- compile(as.grain(dag_mutilated))
P_Y_do_x1_true <- querygrain(grain_do_x1_true, nodes = "Y")$Y["TRUE"]

dag_mutilated_2 <- mutilated(fitted, evidence = list(X1 = "FALSE"))
grain_do_x1_false <- compile(as.grain(dag_mutilated_2))
P_Y_do_x1_false <- querygrain(grain_do_x1_false, nodes = "Y")$Y["TRUE"]

# causal effect of X1 on Y
causal_effect_x1 <- P_Y_do_x1_true - P_Y_do_x1_false
causal_effect_x1


# example usage of grain

yn <- c("yes", "no")
a <- cpt(~asia, values = c(1, 99), levels = yn)
t.a <- cpt(~ tub | asia, values = c(5, 95, 1, 99), levels = yn)
s <- cpt(~smoke, values = c(5, 5), levels = yn)
l.s <- cpt(~ lung | smoke, values = c(1, 9, 1, 99), levels = yn)
b.s <- cpt(~ bronc | smoke, values = c(6, 4, 3, 7), levels = yn)
e.lt <- cpt(
  ~ either | lung:tub,
  values = c(1, 0, 1, 0, 1, 0, 0, 1),
  levels = yn
)
x.e <- cpt(~ xray | either, values = c(98, 2, 5, 95), levels = yn)
d.be <- cpt(
  ~ dysp | bronc:either,
  values = c(9, 1, 7, 3, 8, 2, 1, 9),
  levels = yn
)

chest_cpt <- compile_cpt(a, t.a, s, l.s, b.s, e.lt, x.e, d.be)
chest_cpt

chest_cpt$either

# my own mutual information calculations for the conjunctive case

X1 <- cpt(~X1, values = c(0.6, 0.4), levels = c("B", "W"))
X2 <- cpt(~X2, values = c(0.9, 0.1), levels = c("B", "W"))
Y  <- cpt(~Y | X1:X2, values = c(1, 0, 0, 1,0, 1, 0, 1), levels = c("Win", "Lose"))

conj_network <- grain(compile_cpt(X1, X2, Y))
# query marginal probabilty of Y
querygrain(conj_network, nodes = "Y")$Y
conj_network |> querygrain(nodes = c("X1", "X2", "Y"), type = "joint", simplify = TRUE) |> 
  as_tibble() |> 
  arrange(X1) |> 
  select(X1, X2, Y, Freq)

# compute mutual information between X1 and Y
joint <- querygrain(conj_network, nodes = c("X1", "Y"), type = "joint")
p_x1 <- margin.table(joint, 1)
p_y <- margin.table(joint, 2)
outer <- outer(p_x1, p_y)

mi_X1_Y <- sum(ifelse(joint > 0, joint * log(joint / outer), 0))
mi_X1_Y

# compute mutual information between X' and Y
# X' = 1 iff X1=B and X2=B (conjunctive cause)
joint_X1_X2_Y <- querygrain(conj_network, nodes = c("X1", "X2", "Y"), type = "joint")

joint_Xp_Y <- matrix(0, nrow = 2, ncol = 2,
  dimnames = list(Xp = c("W", "B"), Y = c("Win", "Lose")))

joint_Xp_Y["B", ] <- joint_X1_X2_Y["B", , "B"]
joint_Xp_Y["W", ] <- joint_X1_X2_Y["W", , "W"] +
  joint_X1_X2_Y["W", , "B"] +
  joint_X1_X2_Y["B", , "W"]

p_xp <- rowSums(joint_Xp_Y)
p_y  <- colSums(joint_Xp_Y)

mi_Xp_Y <- sum(ifelse(joint_Xp_Y > 0,
  joint_Xp_Y * log(joint_Xp_Y / outer(p_xp, p_y)), 0))
mi_Xp_Y
