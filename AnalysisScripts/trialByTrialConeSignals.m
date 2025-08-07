ii = 1;

for i = 1:size(ADT,1)
    P680 = ADT.RedPowerWattsAt680nm(i);
    P543 = ADT.GreenPowerWattsAt543nm(i);
    powers = [P680 P543];
    wvls = [680 543];
    
    rawDataIdx = strcmpi(RDT.SubjectID, ADT.SubjectID(i)) & strcmpi(RDT.Folder, ADT.Folder(i)) & ismember(RDT.ExperimentTime, ADT.ExperimentTimes{i}) & RDT.Location == ADT.Location(i);
    tempRawData = RDT(rawDataIdx,:);
     
    for j = 1:size(tempRawData,1)
        Lsig(ii,:) = powers(tempRawData.Channel(j)).*tempRawData.IntensityAU(j).*tempRawData.L(j).*L(wvl==wvls(tempRawData.Channel(j)));
        Msig(ii,:) = powers(tempRawData.Channel(j)).*tempRawData.IntensityAU(j).*tempRawData.M(j).*M(wvl==wvls(tempRawData.Channel(j)));
        LminusM(ii,:) = Lsig(ii,:)-Msig(ii,:);
        LplusM(ii,:) = Lsig(ii,:) + Msig(ii,:);
        ii = ii + 1;
    end
    
    
end