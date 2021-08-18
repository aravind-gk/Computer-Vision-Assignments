% Q1

m_pt = importdata('matches.txt');
m_pt_1 = m_pt(:,1:2);
m_pt_2 = m_pt(:,3:4);
K1 = importdata('calibration_cam1.txt');
K2 = importdata('calibration_cam2.txt');
motion = importdata('motion.txt');
Rg = motion(:,1:3);
tg = motion(:,4);
tg = tg/norm(tg);

[E,HE,F,HF] = eightPoint(m_pt_1, m_pt_2, K1, K2);

[Re,te] = Rotation_and_trans(E);

[h_Re , h_te] = Rotation_and_trans(HE);

[theta,n] = theta_and_axis(((Rg)\Re));
[h_theta,h_n] = theta_and_axis(((Rg)\h_Re));

c1 = dot(tg,te);
c2 = dot(tg,h_te);






