function D_PhiM_pm = get_rotation_over_position(c)

    PhiM = c.OM;
    M0 = c.M0;
    D_PhiM_pm = cell(1,c.NM);
    for m = 1:c.NM
        % w.r.t xm ***********
        t_temp = [-sind(PhiM(1))*cosd(PhiM(2)); ...
                    -sind(PhiM(1))*sind(PhiM(2))*sind(PhiM(3)) - cosd(PhiM(1))*cosd(PhiM(3)); ...
                    cosd(PhiM(1))*sind(PhiM(3))-sind(PhiM(1))*sind(PhiM(2))*cosd(PhiM(3))];
        D_alphaM_xm = M0(:,m)'*t_temp;
        t_temp = [-cosd(PhiM(1))*sind(PhiM(2)); cosd(PhiM(1))*cosd(PhiM(2))*sind(PhiM(3)); ...
                    cosd(PhiM(1))*cosd(PhiM(2))*cosd(PhiM(3))];
        D_betaM_xm  = M0(:,m)'*t_temp;
        t_temp = [0; cosd(PhiM(1))*sind(PhiM(2))*cosd(PhiM(3)) + sind(PhiM(1))*sind(PhiM(3)); ...
                    sind(PhiM(1))*cosd(PhiM(3)) - cosd(PhiM(1))*sind(PhiM(2))*sind(PhiM(3))];
        D_gammaM_xm = M0(:,m)'*t_temp;

        % w.r.t ym ************
        t_temp = [cosd(PhiM(1))*cosd(PhiM(2)); ...
                    -cosd(PhiM(3))*sind(PhiM(1)) + cosd(PhiM(1))*sind(PhiM(2))*sind(PhiM(3)); ...
                    cosd(PhiM(3))*cosd(PhiM(1))*sind(PhiM(2)) + sind(PhiM(1))*sind(PhiM(3))];
        D_alphaM_ym = M0(:,m)'*t_temp;
        t_temp = [-sind(PhiM(1))*sind(PhiM(2)); sind(PhiM(1))*cosd(PhiM(2))*sind(PhiM(3)); ...
                    cosd(PhiM(3))*sind(PhiM(1))*cosd(PhiM(2))];
        D_betaM_ym = M0(:,m)'*t_temp;
        t_temp = [0; -sind(PhiM(3))*cosd(PhiM(1)) + sind(PhiM(1))*sind(PhiM(2))*cosd(PhiM(3)); ...
                    -sind(PhiM(3))*cosd(PhiM(1))*sind(PhiM(2)) - cosd(PhiM(1))*cosd(PhiM(3))];
        D_gammaM_ym = M0(:,m)'*t_temp;

        % w.r.t zm ************
        D_alphaM_zm = 0;
        t_temp = [-cosd(PhiM(2)); -sind(PhiM(2))*sind(PhiM(3)); -sind(PhiM(2))*cosd(PhiM(3))];
        D_betaM_zm = M0(:,m)'*t_temp;
        t_temp = [0; cosd(PhiM(2))*cosd(PhiM(3)); -cosd(PhiM(2))*sind(PhiM(3))];
        D_gammaM_zm = M0(:,m)'*t_temp;

        D_PhiM_xm = [D_alphaM_xm D_betaM_xm D_gammaM_xm];
        D_PhiM_ym = [D_alphaM_ym D_betaM_ym D_gammaM_ym];
        D_PhiM_zm = [D_alphaM_zm D_betaM_zm D_gammaM_zm];

        D_PhiM_pm{m} = [D_PhiM_xm; D_PhiM_ym; D_PhiM_zm];
    end

end