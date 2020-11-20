classdef PLASTICITY1D
  %------------------------------------------------------------------------
  properties
    fy; %yield stress
    E; %young mod.
    K; %iso plastic mod;
    H; %kin plastic mod;
  end
  %------------------------------------------------------------------------
  methods
    %----------------------------------------------------------------------
    function this = PLASTICITY1D(varargin)
      this.fy = 0;
      this.E = 0;
      this.K = 0;
      this.H = 0;
      if (nargin ~= 0)
        cont = 1;
        while cont <= nargin
          switch varargin{cont}
            case 'fy'
              cont = cont + 1;
              this.fy = varargin{cont};
            case 'E'
              cont = cont + 1;
              this.E = varargin{cont};
            case 'K'
              cont = cont + 1;
              this.K = varargin{cont};
            case 'H'
              cont = cont + 1;
              this.H = varargin{cont};
          end
          cont = cont + 1;
        end
      end
    end
    %----------------------------------------------------------------------
    function [f] = yield(this, stress, alpha, q)
      f = abs(stress - q) -(this.fy + this.K * alpha);
    end
    %----------------------------------------------------------------------
    function [df] = dyield(this, stress, q)
      csi = stress - q;
      df = 0;
      if (abs(csi) == 0)
        return
      end
      df = csi / abs(csi);
    end
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    function [stress, Dep, statev] = getConstMat(this, strain, dstrain, ...
                                                                    statev)
      strain = strain + dstrain;
      pstrain = statev(1);
      alpha = statev(2);
      q = statev(3);
      Dep = this.E;
      %Tr
      stress = Dep * (strain - pstrain);
      f = this.yield(stress, alpha, q);
      if (f > 0) %plasticity
        %return mapping
        gamma = f / (this.E + this.K + this.H);
        df = this.dyield(stress, q);
        stress = stress - gamma * Dep * df;
        pstrain = pstrain + gamma * df;
        if (this.K ~= 0)
          alpha = alpha + gamma;
        end
        if (this.H ~= 0)
          q = q + gamma * this.H;
        end
        Dep = this.E * (this.K + this.H) / (this.E + this.K + this.H);
      end
      statev = [pstrain; alpha; q;];
    end
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
  end
  %------------------------------------------------------------------------
end