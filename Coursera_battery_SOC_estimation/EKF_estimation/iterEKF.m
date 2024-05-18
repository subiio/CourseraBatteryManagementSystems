function [zk,zkbnd,ekfData] = iterEKF(vk,ik,Tk,deltat,ekfData)
    model = ekfData.model;

    Q = getParamESC('QParam',Tk,model);
    G = getParamESC('GParam',Tk,model);
    M = getParamESC('MParam',Tk,model);
    M0 = getParamESC('M0Param',Tk,model);
    RC = exp(-deltat./abs(getParamESC('RCParam',Tk,model)))';
    R = getParamESC('RParam',Tk,model)';
    R0 = getParamESC('R0Param',Tk,model);
    eta = getParamESC('etaParam',Tk,model);
    if ik<0, ik = ik*eta; end;

    I = ekfData.priorI;
    SigmaX = ekfData.SigmaX;
    SigmaV = ekfData.SigmaV;
    SigmaW = ekfData.SigmaW;
    xhat = ekfData.xhat;
    irInd = ekfData.irInd;
    hkInd = ekfData.hkInd;
    zkInd = ekfData.zkInd;
    if abs(ik)>Q/100, ekfData.signIk = sign(ik); end;
    signIk = ekfData.signIk;

    nx = length(xhat); Ahat = zeros(nx,nx); Bhat = zeros(nx,1);
    Ahat(zkInd,zkInd) = 1; Bhat(zkInd) = -deltat/(3600*Q);
    Ahat(irInd,irInd) = diag(RC); Bhat(irInd) = 1-RC(:);
    Ah = exp(-abs(I*G*deltat/(3600*Q)));
    Ahat(hkInd,hkInd) = Ah;
    B = [Bhat, 0*Bhat];
    Bhat(hkInd) = -abs(G*deltat/(3600*Q))*Ah*(1+sign(I)*xhat(hkInd));
    B(hkInd,2) = Ah - 1;
    xhat = Ahat*xhat + B*[I; sign(I)];

    SigmaX = Ahat*SigmaX*Ahat' + Bhat*SigmaW*Bhat';
    
    yhat  = OCVfromSOCtemp(xhat(zkInd),Tk,model) + M0*signIk + ...
        M*xhat(hkInd) - R*xhat(irInd) - R0*ik;

    Chat = zeros(1,nx);
    Chat(zkInd) = dOCVfromSOCtemp(xhat(zkInd),Tk,model);
    Chat(hkInd) = M;
    Chat(irInd) = -R;
    Dhat = 1;
    SigmaY = Chat*SigmaX*Chat' + Dhat*SigmaV*Dhat';
    L = SigmaX*Chat'/SigmaY;

    r = vk - yhat;
    if r^2 > 100*SigmaY, L(:)= 0.0; end
    xhat = xhat + L*r;
    xhat(hkInd) = min(1,max(-1,xhat(hkInd)));
    xhat(zkInd) = min(1.05,max(-0.05,xhat(zkInd)));

    SigmaX = SigmaX - L*SigmaY*L';
    if r^2 > 4*SigmaY,
        fprintf('Bumping SigmaX\n');
        SigmaX(zkInd,zkInd) = SigmaX(zkInd,zkInd)*ekfData.Qbump;
    end
    [~,S,V] =svd(SigmaX);
    HH = V*S*V';
    SigmaX = (SigmaX + SigmaX' + HH + HH')/4;
    ekfData.priorI = ik;
    ekfData.SigmaX = SigmaX;
    ekfData.xhat = xhat;
    zk = xhat(zkInd);
    zkbnd = 3*sqrt(SigmaX(zkInd,zkInd));
end

