args <- commandArgs(trailingOnly = TRUE)

# Input RData files and output file
file1 <- args[1]
file2 <- args[2]
output_file <- args[3]

# Load the first RData file
load(file1)
obj1 <- obj_db  # Replace 'obj_db' with the name of the variable in file1

# Load the second RData file
load(file2)

# Combine the objects row-wise
obj_db <- rbind(obj1, obj_db)

# Save the combined object to an RData file
rm(obj1)
save.image(file = output_file)
