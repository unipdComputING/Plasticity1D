clear all; clc;
% perfect plasticity
%mat = PLASTICITY1D('fy', 200.0, 'E', 210000.0);
% iso hardening plasticity
%mat = PLASTICITY1D('fy', 200.0, 'E', 210000.0, 'K', 1000.0);
% kine hardening plasticity
%mat = PLASTICITY1D('fy', 200.0, 'E', 210000.0, 'H', 1500.0);
% iso + kine hardening plasticity
mat = PLASTICITY1D('fy', 200.0, 'E', 210000.0, 'K', 1000.0, 'H', 1500.0);

%--------------------------------------------------------------------------
strain = 0.0;
dstrain = 0.0001;
totInc = 300;
loadUnload = [100; 200; 300;];
%--------------------------------------------------------------------------
statev  = zeros(3, 1);
STRESS  = zeros(totInc, 1);
STRAIN  = zeros(totInc, 1);
PSTRAIN = zeros(totInc, 1);
ALPHA   = zeros(totInc, 1);
BACK    = zeros(totInc, 1);
DEP     = zeros(totInc, 1);
iLU = 1;
totLU = size(loadUnload, 1);
for i = 1: totInc
  [STRESS(i), DEP(i), statev] = mat.getConstMat(strain, dstrain,statev);
  strain = strain + dstrain;
  STRAIN(i) = strain;
  PSTRAIN(i) = statev(1);
  ALPHA(i) = statev(2);
  BACK(i) = statev(3);
  if ((iLU <= totLU) && (i > loadUnload(iLU)))
    dstrain = -1 * dstrain;
    iLU = iLU + 1;
  end
end

plot(STRAIN, STRESS, 'k','LineWidth', 2);

