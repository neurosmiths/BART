# BART
Repository for code to process BART data

This code will allow you to demo:
(a) BART task behavior analysis
(b) Temporal Difference learning model that is applied to each individual's task behavior.

Instructions:

1. Navigate to functions.
2. Open BART_behavior.mat
3. Uncomment and load "demo" file.
4. Run script to get task behavior analysis.
This function reproduces a figure from an example participant. We examine the differences between inflation times of passive and active balloon trials, to determine each participant's impulsivity level, using the Kullback-Liebler Divergence.

1. Navigate to functions.
2. Open BART_behaviorTDlearn.mat
3. Uncomment and load "demo" file.
4. Run script to get task behavior and TD model analysis for an example participant.
This function uses TDlearn.mat to calculate each participant's TD learn outcomes: Value Expectation, Reward Prediction Error, Risk Expectation, and Risk Prediction Error. TDlearn.mat is where we highlight the temporal difference equations we use for temporal difference, unsigned, and asymmetric learning models.
