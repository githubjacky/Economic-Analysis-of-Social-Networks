################################################################################
###                                                                          ###
### ---- basicRSiena.r: a basic script for the introduction to RSiena ------ ###
###                                                                          ###
###                               version: February 1, 2021                  ###
################################################################################

# An example of a basic sequence of commands
# for estimating a basic model by function siena07() of RSiena.
# With a lot of use of help pages; this can be skipped as you like.
# Note that lines starting with # are comment lines, not commands.

# What is your current working directory?
getwd()
# If you wish it to be different, change it by
setwd('/Users/jackyyeh/github/Economic-Analysis-of-Social-Networks/course/Dynamic-network-formation')
# where you can fill in your own directory between the parentheses.
# On Windows you will have to use double instead of single backslashes.

# If you have not yet installed RSiena, install it:
# install.packages("RSiena", repos="https://cloud.r-project.org")
# (if you do not specify the repos,
#  you get a choice from a list of available repositories.)
# This installation will persist on your machine, no need to repeat it later.
# To install the recent version from GitHub go to
# https://github.com/snlab-nl/rsiena/releases/
# You find instructions for installation there or in the RSiena manual,
# p. 17 of http://www.stats.ox.ac.uk/~snijders/siena/RSiena_Manual.pdf
# The GitHub version is more up to date than the version
# available from CRAN, the general R archive.
# Note that you can also download the recent versions also from the
# Siena website, http://www.stats.ox.ac.uk/~snijders/siena/siena_downloads.htm
# and then install from the local file.

# Define data sets

# If you have internet access, you can download the data
# from the Siena website ("Data sets" tab)
# http://www.stats.ox.ac.uk/~snijders/siena/s50_data.zip
# and unzip it in your working directory.
# The data description is at
# http://www.stats.ox.ac.uk/~snijders/siena/s50_data.htm
# Then you can read the data files by the commands
# (this can be replaced by using the internal data set, see below)
#        friend.data.w1 <- as.matrix(read.table("s50-network1.dat"))
#        friend.data.w2 <- as.matrix(read.table("s50-network2.dat"))
#        friend.data.w3 <- as.matrix(read.table("s50-network3.dat"))
#        drink <- as.matrix(read.table("s50-alcohol.dat"))
#        smoke <- as.matrix(read.table("s50-smoke.dat"))

# But without internet access, the data can be obtained
# from within RSiena (see below), because this is an internal data set.

library(RSiena)
# With
?RSiena
# you get to see the basic operation.
# Now we use the internally available s50 data set.
# Look at its description:
?s50
# 3 waves, 50 actors
# Look at the start and end of the first wave matrix
head(s501)
tail(s501)
# and at the alcohol variable
s50a
# Now define the objects with the same names as above
# (this step is superfluous if you read the data already).
friend.data.w1 <- s501
friend.data.w2 <- s502
friend.data.w3 <- s503
drink <- s50a
smoke <- s50s

# Now the data must be given the specific roles of variables
# in an RSiena analysis.
#
# Dependent variable
####################
?sienaDependent
# First create a 50 * 50 * 3 array composed of the 3 adjacency matrices
friendshipData <- array( c( friend.data.w1, friend.data.w2, friend.data.w3 ),
           dim = c( 50, 50, 3 ) )
# and next give this the role of the dependent variable:
friendship <- sienaDependent(friendshipData)
# What did we construct?
friendship

# We also must prepare the objects that will be the explanatory variables.
#
# Actor covariates
##################
# We use smoking for wave 1 as a constant actor covariate:
smoke1 <- coCovar( smoke[ , 1 ] )
# A variable actor covariate is defined for drinking:
alcohol <- varCovar( drink )
# (This choice is purely for the purpose of illustration here.)

# Put the variables together in the data set for analysis
?sienaDataCreate
mydata <- sienaDataCreate( friendship, smoke1, alcohol )
# Check what we have
mydata

# You can get an outline of the data set with some basic descriptives from
print01Report( mydata, modelname="s50")
# For the model specification we need to create the effects object
myeff <- getEffects( mydata )
# All the effects that are available given the structure
# of this data set can be seen from
effectsDocumentation(myeff)
# For a precise description of all effects, see Chapter 12 in the RSiena manual.
# A basic specification of the structural effects:
?includeEffects
myeff <- includeEffects( myeff, transTrip, cycle3)
# and some covariate effects:
myeff <- includeEffects( myeff, egoX, altX, simX, interaction1 = "alcohol" )
myeff <- includeEffects( myeff, simX, interaction1 = "smoke1" )
myeff

# Create object with algorithm settings
# Accept defaults but specify name for output file
# (which you may replace by any name you prefer)
?sienaAlgorithmCreate
myalgorithm <- sienaAlgorithmCreate( projname = 's50' )

# Estimate parameters
?siena07
ans <- siena07( myalgorithm, data = mydata, effects = myeff)
ans
# This gives results from a random starting point.
# To use a fixed starting point, use the "seed" parameter:
# myalgorithm <- sienaAlgorithmCreate( projname = 's50', seed=435123 )

# For checking convergence, look at the
# 'Overall maximum convergence ratio' mentioned under the parameter estimates.
# It can also be shown by requesting
ans$tconv.max
# If this is less than 0.25, convergence is good.
# If convergence is inadequate, estimate once more,
# using the result obtained as the "previous answer"
# from which estimation continues:

ans <- siena07( myalgorithm, data = mydata, effects = myeff, prevAns=ans)
ans
# If convergence is good, you can look at the estimates.
# More extensive results
summary(ans)
# Still more extensive results are given in the output file
# s50.txt in the current directory.

# Note that by putting an R command between parentheses (....),
# the result will also be printed to the screen.
# Next add the transitive reciprocated triplets effect,
# an interaction between transitive triplets and reciprocity,

(myeff <- includeEffects( myeff, transRecTrip))
(ans1 <- siena07( myalgorithm, data = mydata, effects = myeff, prevAns=ans))
# If necessary, repeat the estimation with the new result:
(ans1 <- siena07( myalgorithm, data = mydata, effects = myeff, prevAns=ans1))
# This might still not have an overall maximum convergence ratio
# less than 0.25. If not, you could go on once more.
#
# Inspect the file s50.txt in your working directory
# and understand the meaning of its contents.

# To have a joint test of the three effects of alcohol:
?Multipar.RSiena
Multipar.RSiena(ans1, 7:9)
# Focusing on alcohol similarity, the effect is significant;
# diluting the effects of alcohol by also considering ego and alter,
# the three effects simultaneously are not significant.

################################################################################
###                Assignment 1                                              ###
################################################################################

# 1a.
# Drop the effect of smoke1 similarity and estimate the model again.
# Do this by the function setEffect() using the <<include>> parameter.
# Give the changed effects object and the new answer object new names,
# such as effects1 and ans1, to distinguish them.
eff1a <- setEffect(myeff, simX, interaction1 = "smoke1", include = FALSE)
ans1a <- siena07( myalgorithm, data = mydata, effects = eff1a, verbose=FALSE)
summary(ans1a)
# 1b.
# Change the three effects of alcohol to the single effect
# of alcohol similarity, and estimate again.
eff1b <- setEffect(myeff, egoX, interaction1 = "alcohol", include = FALSE)
eff1b <- setEffect(eff1b, altX, interaction1 = "alcohol", include = FALSE)
ans1b <- siena07( myalgorithm, data = mydata, effects = eff1b, verbose=FALSE)
summary(ans1b)

################################################################################
###                Networks and behavior study                               ###
################################################################################

# Now we redefine the role of alcohol drinking
# as a dependent behaviour variable.
# Once again, look at the help file
?sienaDependent
# now paying special attention to the <<type>> parameter.
drinking <- sienaDependent( drink, type = "behavior" )

# Put the variables together in the data set for analysis
NBdata <- sienaDataCreate( friendship, smoke1, drinking )
NBdata
NBeff <- getEffects( NBdata )
effectsDocumentation(NBeff)
NBeff <- includeEffects( NBeff, transTrip, transRecTrip )
NBeff <- includeEffects( NBeff, egoX, egoSqX, altX, altSqX, diffSqX,
                         interaction1 = "drinking" )
NBeff <- includeEffects( NBeff, egoX, altX, simX, interaction1 = "smoke1" )
NBeff
# For including effects also for the dependent behaviour variable, see
# ?includeEffects
NBeff <- includeEffects( NBeff, avAlt, name="drinking",
                         interaction1 = "friendship" )
NBeff
# Define an algorithm with a new project name
myalgorithm1 <- sienaAlgorithmCreate( projname = 's50_NB' )

# Estimate again, using the second algorithm right from the start.
NBans <- siena07( myalgorithm1, data = NBdata, effects = NBeff)
# You may improve convergence (considering the overall maximum
# convergence ratio) by repeated estimation in the same way as above.

# Look at results
NBans
# Make a nicer listing of the results
siena.table(NBans, type="html", sig=TRUE)
# This produces an html file; siena.table can also produce a LaTeX file.

################################################################################
###                Assignment 2                                              ###
################################################################################

# 2a.
# Replace the average alter effect by average similarity (avSim)
# or total similarity (totSim) and estimate the model again.
eff2a <- setEffect(NBeff, avAlt, name="drinking", 
                   interaction1 = "friendship", include = FALSE)
eff2a <- includeEffects(eff2a, avSim, name="drinking",
                        interaction1 = "friendship" )
ans2a <- siena07( myalgorithm1, data = NBdata, effects = eff2a)
symmary(ans2a)
# 2b.
# Add the effect of smoking on drinking and estimate again.
eff2b <- includeEffects(NBeff, egoX, altX, simX, name="drinking", 
                        interaction1="smoke1")
ans2b <- siena07(myalgorithm1, data = NBdata, effects = eff2b,
                 verbose=FALSE, silent=TRUE)
summary(ans2b)

################################################################################

################################################################################
###                Assignment 3                                              ###
################################################################################

# Read Sections 13.3 and 13.4 of the Siena Manual, download scripts
# SelectionTables.r and InfluenceTables.r from the Siena website,
# and make plots of the selection table and influence table for drinking.

################################################################################
