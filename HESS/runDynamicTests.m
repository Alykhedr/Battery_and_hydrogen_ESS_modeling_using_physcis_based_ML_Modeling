function runDynamicTests()
  % Define test profiles
  profiles = { ...
    struct('P_el',100,'P_fc',0,'desc','Constant Charge'); ...
    struct('P_el',0,'P_fc',100,'desc','Constant Discharge'); ...
    struct('P_el',0,'P_fc',0,'desc','Zero Power'); ...
  };

  % Loop through each profile
  for i=1:numel(profiles)
    time = 0:1:10;   % [0…10] hr
    P_el = profiles{i}.P_el * ones(size(time));
    P_fc = profiles{i}.P_fc * ones(size(time));

    % Simulate Model 1
    p         = params(time);
    m_dot_el  = p.init.m_dot_el;
    m_dot_fc  = p.init.m_dot_fc;
    mass_H2   = p.init.mass_H2;
    dt        = 1;  % hr

    for k = 2:numel(time)
      if P_el(k)>0
        m_dot_el(k) = p.eta_el * P_el(k) / p.LHV_H2;
      end
      if P_fc(k)>0
        m_dot_fc(k) = P_fc(k) / (p.eta_fc * p.LHV_H2);
      end
      mass_H2(k) = mass_H2(k-1) + (m_dot_el(k)-m_dot_fc(k)) * dt;
      mass_H2(k) = min(max(mass_H2(k),0), p.max_mass_H2);
    end

    % Compute expected final mass
    switch profiles{i}.desc
      case 'Constant Charge'
        expected = p.init.mass_H2(1) + ...
                   (p.eta_el * profiles{i}.P_el / p.LHV_H2) * time(end);
      case 'Constant Discharge'
        expected = p.init.mass_H2(1) - ...
                   (profiles{i}.P_fc / (p.eta_fc * p.LHV_H2)) * time(end);
        expected = max(expected, 0);
      case 'Zero Power'
        expected = p.init.mass_H2(1);
    end

    % Check error
    actual = mass_H2(end);
    err    = abs(actual - expected);

    % Report
    if err < 1e-6
      status = 'PASS';
    else
      status = 'FAIL';
    end

    fprintf('%-18s Test: actual=%6.2f kg, expected=%6.2f kg, err=%.2e → %s\n', ...
      profiles{i}.desc, actual, expected, err, status);
  end
end
