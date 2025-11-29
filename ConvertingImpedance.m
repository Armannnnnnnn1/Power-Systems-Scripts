function [Xpu_new, X_SI] = ConvertingImpedance(Vb_new, Sbnew, Xpu_old, Sb_old, Vb_old)
%CONVERTINGIMPEDANCE Converts per-unit impedance when changing power base
%   Inputs:
%   Vbnew - New voltage base (kV)
%   Sbnew - New apparent power base (MVA)  
%   Zb_old - Old base impedance (ohms)
%   Xpu_old - Per-unit impedance with old base
%
%   Outputs:
%   Xpu_new - Per-unit impedance with new base
%   Zb_new - New base impedance (ohms)
%   X_SI - Actual impedance in ohms

Zb_new = (Vb_new*Vb_new) / Sbnew;
Zb_old = (Vb_old^2) / Sb_old;
X_SI = Xpu_old * Zb_old;
Xpu_new = X_SI / Zb_new;
end
