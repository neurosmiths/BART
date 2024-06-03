# BART
Repository for code to process BART data

This code will allow you to demo:
(a) BART task behavior analysis
(b) Temporal Difference learning model that is applied to each individual's task behavior.

Instructions:

A. Navigate to DemoData.
   1. Here you can open the .ns2 file from a public access boxfolder. This could not be uploaded to github, due to the file size. Please download this file to run TDlearn files.
   2. You can also find the behavioral and .nev data files needed here to run the folliwng scripts.

B. Behavioral Data Example:
  1. Navigate to functions.
  2. Open BART_behavior.mat
  3. Load behavioral and neural data using the uigetfile function
  4. Run script to get task behavior analysis.
This function reproduces a figure from an example participant. We examine the differences between inflation times of passive and active balloon trials, to determine each participant's impulsivity level, using the Kullback-Liebler Divergence.

C. Temporal Difference Data Example:
  1. Navigate to functions.
  2. Open BART_behaviorTDlearn.mat
  3. Load behavioral and neural data using the uigetfile function
  4. Run script to get task behavior and TD model analysis for an example participant.

D. Look at TDlearn.mat
  1. Navigate to resources (rsrcs) folder.
  2. Open TDlearn.mat
  3. Examine code. This is the primary code for the TD models we used to calculate learning trajectories during the task. 
This function calculates each participant's TD learn outcomes: Value Expectation, Reward Prediction Error, Risk Expectation, and Risk Prediction Error. TDlearn.mat is where we highlight the temporal difference equations we use for temporal difference, unsigned, and asymmetric learning models.

