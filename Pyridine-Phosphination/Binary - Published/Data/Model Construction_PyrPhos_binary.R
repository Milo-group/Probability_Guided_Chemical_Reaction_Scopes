"
# Install
remotes::install_github('https://github.com/barkais/rxn.cond.class.git')
"

library('rxn.cond.class')

# Halogenation through selective phosphination ----------------------------

#----------------#
#---  Binary ----#
#----------------#

# Import training data from a CSV file
data_bin <- data.frame(data.table::fread('Training_Data_rdkit_stereoelectronic.csv'), check.names = F)

# Set the first column as row names
row.names(data_bin) <- data_bin[,1]  
data_bin <- data_bin[,-1]  # Remove the first column after setting row names

# Convert the 'class' column to a factor
data_bin$class <- as.factor(data_bin$class)

# Convert the first 50 columns to numeric values
data_bin[,c(1:50)] <- as.numeric(as.matrix(data_bin[,c(1:50)]))

# Convert columns 51 to 64 to factors after converting to numeric
for (i in 51:64) data_bin[,i] <- as.factor(as.numeric(data_bin[,i]))

# Add a 'flag' column to sequentially number the rows
data_bin <- plyr::mutate(data_bin, flag = seq(1,nrow(data_bin)))

# Sample from group 1 (first 50 columns + flag) with sample size being 75% of group 1 size
one <- simi.sampler(data_bin[, c(1:50, 65:66)], 1, sample.size = round(sum(data_bin$class == 1) * 0.75))

# Sample from group 2 (first 50 columns + flag) with sample size being 75% of group 2 size
two <- simi.sampler(data_bin[, c(1:50, 65:66)], 2, sample.size = round(sum(data_bin$class == 2) * 0.75))

# Define train and test data from the samples taken from groups 1 and 2
Train.data <- data_bin[c(one, two), ]
Test.data <- data_bin[-c(one, two), ]

# Non-ordinal multinomial logistic regression models
models <- sub_model_log(data = Train.data,
                        min = 4,
                        max = 4, 
                        ordinal = F)

test.form <- models[5,1]
# To the 5th model we added the `Avg_NPA_SM`  feature
# test.form <- "`class` ~ `Avg_NPA_SM` + `fr_Ar_N` + `fr_aryl_methyl` + `fr_benzene` + `fr_halogen`"

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
        conformation = '5. 5th Place')$plot

# Heatmap
prob.heatmap(test, Train.data, 
             plot.title = 'Training Set', 
             conformation = '5. 5th Place')

# ---------- #
# -- Test -- #
# ---------- #

# Confusion Matrix
model.info <- mod.info(test, Test.data, FALSE, FALSE)

ct_plot(model.info$class.table, 
        plot.title = 'Test Set', 
        conformation = '5. 5th Place')$plot

# Heatmap
prob.heatmap(test, Test.data, 
             plot.title = 'Test Set', 
             conformation = '5. 5th Place')
