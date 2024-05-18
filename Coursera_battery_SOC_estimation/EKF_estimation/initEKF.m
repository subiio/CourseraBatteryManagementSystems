function ekfData = initEKF(v0,T0,SigmaX0,SigmaV,SigmaW,model)
    ir0 = 0;    ekfData.irInd = 1;
    hk0 = 0;    ekfData.hkInd = 2;
    SOC0 = SOCfromOCVtemp(v0,T0,model); ekfData.zkInd = 3;
    ekfData.xhat = [ir0 hk0 SOC0]';

    ekfData.SigmaX = SigmaX0;   ekfData.SigmaV = SigmaV;
    ekfData.SigmaW = SigmaW;    ekfData.Qbump = 5;

    ekfData.priorI = 0; ekfData.signIk = 0;

    ekfData.model = model;
end
