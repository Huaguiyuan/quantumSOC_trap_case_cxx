%% plot data
clear all
clc
load photon_number.data
load photon_fluctu.data
i = 1;
for N = 1:10
    for Q = 1:10
        a(N,Q) = photon_number(i);
        b(N,Q) = photon_fluctu(i);
        i = i+1;
    end
end
[aa,bb]=meshgrid(1:10,1:10);
figure(1)
surf(aa,bb,a)
figure(2)
surf(aa,bb,b)