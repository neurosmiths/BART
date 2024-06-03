function [] = deleteShittyBalloonOnsetStats()

[ptArray] = BARTnumbers;

for pt = 1:length(ptArray)
    
    ptID = ptArray{pt};
    
    BARTdir = sprintf('/media/user1/data4TB/data/BART/BART_EMU/%s/Data',ptID);
    fName = '*balloonOnsetHighGammaStats.mat';
    trodeFiles = dir(fullfile(BARTdir,fName));

    for fl = 1:length(trodeFiles)
        
        % deleting the old files. 
        delete(fullfile(trodeFiles(fl).folder,trodeFiles(fl).name))
        
    end
end
