
% CELL_VOLTAGE  Calculates a single-cell voltage [V]

function V = cell_voltage(i, p)
%   i [A/cm^2]   : current density
%   p            : parameter struct from params_hess

  R     = p.R;                % J/(mol·K)
  F     = p.F;                % C/mol
  T_el  = p.el.T;             % K (operating electrolyzer temperature)


  %% ---- Reversible voltage (Nernst equation) ----
E0   = (p.el.DG0_R) / (p.z*p.F);      % Eq. (30)
Vrev = E0 + (R*T_el/(p.z*F)) * log ((p.el.P_H2/p.std.P_cat_ref_bar) * sqrt(p.el.P_O2/p.std.P_an_ref_bar));

%% ---- Membrane conductivity (Bernardi–Verbrugge) ----
  sigma_mem = (F^2 * p.el.C_Hp_mol_cm3 * p.el.D_Hp_cm2_s) / (R * T_el); % [S/cm]
  sigma_mem = max(sigma_mem, 1e-12);

  %% ---- Membrane ohmic overpotential ----
  Vohm = (p.el.delta_mem_cm / sigma_mem) * i; % δ_mem/σ_mem * i [V]

  %% ---- Roughness factor γ_M for anode and cathode ----
  % γ_M = φ_I * m_M / (6 * ρ_M * d_M)  [cm^2/cm^2]
gamma_M_an = (p.el.phiI_an * p.el.mM_an) * (6 / (p.el.rhoM_an * p.el.dM_an)); % anode
gamma_M_ca = (p.el.phiI_ca * p.el.mM_ca) * (6 / (p.el.rhoM_ca * p.el.dM_ca)); % cathode

  %% ---- Temperature-dependent exchange current density (Arrhenius) ----
  % i0* = i0_ref * exp(-Ea/R * (1/T - 1/Tref))  [A/cm^2]
  i0_star_an = p.el.i0_ref_an * exp(-p.el.Ea_an / R * (1/T_el - 1/p.el.Tref)); % anode
  i0_star_ca = p.el.i0_ref_ca * exp(-p.el.Ea_ca / R * (1/T_el - 1/p.el.Tref)); % cathode

  %% ---- Exchange current density including roughness ----
  % i0 = γ_M * i0*
  i0_an = gamma_M_an * i0_star_an; % anode exchange current density [A/cm^2]
  i0_ca = gamma_M_ca * i0_star_ca; % cathode exchange current density [A/cm^2]

  %% ---- Activation overpotential (replace i_p/i_n with i0_an/i0_ca) ----
  Vact = R*T_el/(p.el.alpha_an*F) * asinh(i/(2*i0_an)) + ...
         R*T_el/(p.el.alpha_ca*F) * asinh(i/(2*i0_ca));

  %% ---- Total cell voltage ----
  V = Vrev + Vohm + Vact;

end

