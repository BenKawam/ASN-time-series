# Empirical data

This directory contains the empirical datasets used in [`01.2_core_model_empirical`](https://htmlpreview.github.io/?https://github.com/BenKawam/ASN-time-series/blob/main/01.2_core_model_empirical.html) and [`02.2_extended_model_empirical`](https://htmlpreview.github.io/?https://github.com/BenKawam/ASN-time-series/blob/main/02.2_extended_model_empirical.html) (see source directory).

## Data object `states_tibble.rds`

Each row of this data frame encodes one observation $j$—*i.e.*, a holding time in a dyadic behavioural state.
The features of these observations are encoded by the following columns:

| Variable | Description | Type |
|---|---|---|
| `id_iff` | Unique number identifying a focal-animal sampling protocol, or "individual focal-follow" (`iff`). `id_iff` is shared across several observations $j$ that are recorded during the same `iff`.| integer | / |
| `start_iff` | Start time of the `iff`. | datetime |
| `t` | Time of the action; _i.e._ start of the sojourn. | datetime |
| `focal_animal` |  Identifier of the focal animal. | integer |
| `actor` | Identifier of the behaviour's actor. | integer |
| `receiver` | Identifier of the behaviour's receiver.| integer |
| `gr` | Social group identifier. | integer |
| `focal_dyad` | Dyad identifier. | integer |
| `X` | Holding time. | numeric (minutes) |
| `S` | Dyadic state. | integer |
| `C` | Censoring. | integer |
| `row_num` | Row number. | integer |


## Data object `d_stan_format.rds`

This list contains both scalar values and one-dimensional containers.
Below, we first describe the scalar values, followed by four categories of containers, each corresponding to a different data structure.
The four data structures are combined in the original object but are presented separately here for clarity.

### Scalars

| Variable | Description | Type |
|---|---|---|
| `N_dyad` | Number of dyads. | integer |
| `N_ind` | Number of individuals. | integer |
| `N_group` | Number of social groups. | integer |
| `J` | Number of observations, _i.e._, all dyadic behavioural state sojourns $j$. | integer |
| `N_trans` | Number of observed state transitions. | integer |
| `N_row_D` | `N_dyad * 4`: number of rows of data set `D`. | integer |

### Data set A

Each element corresponds to one dyad; they have length `N_dyad`.

| Variable | Description | Type |
|---|---|---|
| `A_grp` | Social group in which the dyad is nested. | integer |
| `A_ind_A` | Identifier of individual $a$ in the dyad. | integer |
| `A_ind_B` | Identifier of individual $b$ in the dyad. | integer |
| `A_dyad` | Dyad identifier. | integer |
| `A_z_A` | Sex of individual $a$. | integer |
| `A_z_B` | Sex of individual $b$. | integer |

### Data set B

Each element corresponds to one dyadic state sojourn; they have length `J`.

| Variable | Description | Type |
|---|---|---|
| `B_grp` | Social group. | integer |
| `B_x` | Holding time of the dyadic state. | numeric |
| `B_s` | Dyadic behavioural state. | integer |
| `B_c` | Censoring indicator. | integer |
| `B_dyad` | Dyad identifier. | integer |

### Data set C

Each element corresponds to one dyadic state transition; they have length `N_trans`.

| Variable | Description | Type |
|---|---|---|
| `C_grp` | Social group. | integer |
| `C_dyad` | Dyad identifier. | integer |
| `C_s_from` | Current state $k$ from which the transition occurs. | integer |
| `C_s_to` | Future state $l$ to which the transition occurs. | integer |

### Data set D

Each element corresponds to one state (1, 2, 3, or 4) per dyad; they have length `N_row_D`.
They allow us to compute the naive estimates that are analogous to the simple ratio index.

| Variable | Description | Type |
|---|---|---|
| `D_dyad` | Dyad identifier. | integer |
| `D_s` | State (1, 2, 3 or 4). | integer |
| `D_exposure` | Amount of time in this state; _i.e._, sum of holding times of each dyad in each state. | numeric |
| `D_events` | Number of times each dyad was observed leaving each state. | integer |
