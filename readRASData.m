function RasData = readRASData(filepath)
    dataMatrix = readmatrix(filepath,"NumHeaderLines",1);
    rasAlt1 = dataMatrix(:,23);
    [RasData.Apogee,apoindex] = max(rasAlt1);
    fprintf("Maximum RASAero II Altitude: %g ft\n",RasData.Apogee);
    RasData.Altitude = dataMatrix(1:apoindex,23);
    RasData.Time = dataMatrix(1:apoindex,1);
    RasData.ApogeeTime = RasData.Time(apoindex);
    RasData.Mach = dataMatrix(1:apoindex,4);
    RasData.Velocity = dataMatrix(1:apoindex,18);
end