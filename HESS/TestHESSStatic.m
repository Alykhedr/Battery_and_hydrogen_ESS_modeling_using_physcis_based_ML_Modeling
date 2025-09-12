classdef TestHESSStatic < matlab.unittest.TestCase
  methods (Test)
    function testLoadProfilesStructure(tc)
      [t, Pel, Pfc] = load_profiles();

      % 1) Check endpoints: 0 to 170 hr
      tc.verifyEqual(t(1), 0);
      tc.verifyEqual(t(end), 170);

      % 2) Check uniform 1-hr steps
      expectedDiff = ones(1, numel(t)-1);
      tc.verifyEqual(diff(t), expectedDiff);

      % 3) Mutual exclusivity: never both > 0 simultaneously
      tc.verifyTrue(all(~(Pel>0 & Pfc>0)));
    end

    function testParamsDefaults(tc)
      t = 0:1:10;
      p = params(t);
      % constants
      tc.verifyEqual(p.LHV_H2,      33.33);
      tc.verifyEqual(p.eta_el,      0.60);
      tc.verifyEqual(p.eta_fc,      0.49);
      tc.verifyEqual(p.max_mass_H2, 350);
      % initial tank level should be 100 kg
      tc.verifyEqual(p.init.mass_H2(1), 100);
    end

    function testPlotHESSRuns(tc)
      % ensure plot_hess runs without error; supply all required args
      t        = (0:1:5)';
      Pel      = zeros(size(t));
      Pfc      = zeros(size(t));
      mH2      = zeros(size(t));
      m_dot_el = zeros(size(t));
      m_dot_fc = zeros(size(t));
      tc.verifyWarningFree(@() plot_hess(t, Pel, Pfc, mH2, m_dot_el, m_dot_fc));
    end
  end
end
