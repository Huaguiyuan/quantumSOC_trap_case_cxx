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
