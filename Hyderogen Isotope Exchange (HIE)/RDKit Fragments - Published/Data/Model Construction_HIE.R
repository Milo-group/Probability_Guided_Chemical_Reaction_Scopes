"
# Install
remotes::install_github('https://github.com/barkais/rxn.cond.class.git')
"

library('rxn.cond.class')

# Hydrogen Isotope Exchange -----------------------------------------------

# Import training data from CSV and convert it into a data frame
data50 <- data.frame(data.table::fread('Training_Data.csv')) 
row.names(data50) <- data50[,1]  # Set the first column as row names
data50 <- data50[,-1]  # Remove the first column after setting row names

# Convert all columns in the data frame to factors
for(i in 1:ncol(data50)) data50[, i] <- as.factor(data50[, i])

# Add a new column 'flag' that sequentially numbers the rows
data50 <- plyr::mutate(data50, flag = seq(1,nrow(data50)))

## Sampler data

# Import sampling data from CSV and convert to data frame
sampl.dat <- data.table::fread('Sampling_Data.csv') 
sampl.dat <- data.frame(sampl.dat)
row.names(sampl.dat) <- sampl.dat[,1]  # Set the first column as row names
sampl.dat <- sampl.dat[,-1]  # Remove the first column after setting row names

# Convert the 'class' column to a factor
sampl.dat$class <- as.factor(sampl.dat$class)

# Convert the first 60 columns of sampling data to numeric values
sampl.dat[,c(1:60)] <- as.numeric(as.matrix(sampl.dat[,c(1:60)]))

# Add a 'flag' column to sequentially number the rows
sampl.dat <- plyr::mutate(sampl.dat, flag = seq(1,nrow(sampl.dat)))

# Sample from group 1 using simi.sampler, returning samples from the smallest group
one <- simi.sampler(sampl.dat, 1)
two <- simi.sampler(sampl.dat, 2)
three <- simi.sampler(sampl.dat, 3)
four <- simi.sampler(sampl.dat, 4)

# Sample from group 2 molecules that are similar to group 1
two_three <- simi.sampler(sampl.dat, 2, 3)

# Sample from group 1 molecules that are similar to group 3
one_three <- simi.sampler(sampl.dat, 1, 3)

# Sample from group 4 molecules that are similar to group 3
four_three <- simi.sampler(sampl.dat, 4, 3)

## Back to model search

# Combine similarities from the various groups sampled
similarties <- c(three,
                 union(two, two_three),
                 union(one, one_three),
                 union(four, four_three))

# Define training data from the subset of similarities
Train.data <- data50[similarties,]

# Define test data from the remaining rows not in the similarities set
Test.data <- data50[-similarties,]

# Non-ordinal multinomial logistic regression models
models <- sub_model_log(data = Train.data,
                                min = 4,
                                max = 4, 
                                ordinal = F)

test.form <- models[3,1]
# test.form <- "`class` ~ `fr_NH0` + `fr_aldehyde` + `fr_amide` + `fr_pyridine`"

test <- nnet::multinom(test.form,
                       data = Train.data,
                       maxit = 2000, 
                       trace = FALSE)

# --------- #
# --Train-- #
# --------- #

# Confusion Matrix
model.info <- mod.info(test, Train.data, FALSE, FALSE)

ct_plot(model.info$class.table, 
        plot.title = 'Training Set', 
        conformation = '3. 3rd Place')$plot



# Heatmap
prob.heatmap(test, Train.data, 
             plot.title = 'Training Set', 
             conformation = '3. 3rd Place')

# ---------- #
# -- Test -- #
# ---------- #

# Confusion Matrix
model.info <- mod.info(test, Test.data, FALSE, FALSE)

ct_plot(model.info$class.table, 
        plot.title = 'Test Set', 
        conformation = '3. 3rd Place')$plot

# Heatmap
prob.heatmap(test, Test.data, 
             plot.title = 'Test Set', 
             conformation = '3. 3rd Place')

# ----------------- #
# -- Predictions -- #
# ----------------- #

# Load and organize prediction of new substrates data
Prediction.set <- data.frame(data.table::fread('Predictions_Data.csv'), check.names = F)
RN <- Prediction.set$V1
Prediction.set <- Prediction.set[,-1]
row.names(Prediction.set) <- RN
Prediction.set <- lapply(Prediction.set, as.factor)

knitr::kable(cbind(predict(test, Prediction.set, 'probs') * 100,
                   predicted_class = predict(test, Prediction.set, 'class')))
