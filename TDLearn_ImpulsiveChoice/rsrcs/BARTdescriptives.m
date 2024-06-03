function [bhvStruct] = BARTdescriptives(ptID)

% Function to get sex and age descriptives of patients.
% patient details.
ptID = ['202202'];

T = readtable('D:\Data\preProcessed\BART_preprocessed\BART_Subject_Sex_Age.csv');

rowsByName = T(ptID,["Sex","Age"])
rowsByName = T(ptID,["Sex","Age"])




%sex = logical({T.Sex},'M');
if ptID
    sex = T.Sex;
    age = T.Age;
else
    fprintf('\nThis patient may not have descriptives...\n')
end

end
