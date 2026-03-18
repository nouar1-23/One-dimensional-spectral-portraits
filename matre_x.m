function[A0,A1,A3]=matre_x()
main_diag = -1 : -0.1 : -1.9;

upper_diag0 = 2* ones(1, 9);
upper_diag1 = 10* ones(1, 9);
upper_diag2 = 50 * ones(1, 9);
upper_diag3 = [3^8, 3^6, 3^3, 3^2, 3^2, 3, 3, 3^-1, 3^-1];
D0 = diag(main_diag) + diag(upper_diag0, 1);
D1 = diag(main_diag) + diag(upper_diag1, 1);
D2 = diag(main_diag) + diag(upper_diag2, 1);

D3 = diag(main_diag) + diag(upper_diag3, 1);
L_main = ones(1, 10);

L1_sub = -1 * ones(1, 9);
L2_sub = -2 * ones(1, 9);

L1 = diag(L_main) + diag(L1_sub, -1);
L2 = diag(L_main) + diag(L2_sub, -1);
in_v_L1=inv(L1);
in_v_L2=inv(L2);


A0 = L1 * D0 * in_v_L1;
A1 = L1 * D1 * in_v_L1;
A2 = L1 * D2 * in_v_L1;
A3 = L2 * D1 * in_v_L2;
A4 = L2 * D2 * in_v_L2;
A5 = L1 * D3 * in_v_L1;
