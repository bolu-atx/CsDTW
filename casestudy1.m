% case study 1.
clc;clear all;close all;
load testdata;

[M,i] = max(size(ref));
if i ~= 1;ref = ref';end

[N,i] = max(size(target));
if i ~= 1;target = target';end

subplot(311);ylim([0 0.1])
plot(ref,'-','LineWidth',1.5);hold on;plot(target,':r','LineWidth',0.5);
% simple interpolation
target_s = interp1(linspace(1,N,N),target,linspace(1,N,M)');
plot(target_s,':g','LineWidth',1.5)
xlabel('Time');
ylabel('Amplitude');
legend('Reference','Target','Warped');
title('Interpolation');
ylim([0 0.1]);xlim([50 300]);

subplot(312);
plot(ref,'-','LineWidth',1.5);hold on;plot(target,':r','LineWidth',0.5);

[Dist, D,k,w,target_s] = dtw(ref,target,0);
plot(target_s,':g','LineWidth',1.5)
xlabel('Time');
ylabel('Amplitude');
title('Dynamic Time Warping (unconstrained)');
ylim([0 0.1]);xlim([50 300]);

subplot(313);ylim([0 0.1])
plot(ref,'-','LineWidth',1.5);hold on;plot(target,':r','LineWidth',0.5);
[Dist,D,d,k,w,target_s] = LCdtw(ref,target,ones(numel(ref),numel(target)),2,0);
plot(target_s,':g','LineWidth',1.5)
xlabel('Time');
ylabel('Amplitude');
title('Dynamic Time Warping (with slope constraint (0.5,2))');
ylim([0 0.1]);xlim([50 300]);
%%
figure;
subplot(311);ylim([0 0.1])
plot(ref,'-','LineWidth',1.5);hold on;plot(target,':r','LineWidth',0.5);
% apply my warping method
warped_result = mywarping(ref,target,[150],[250],0);
target_s = warped_result.tww;

plot(target_s,':g','LineWidth',1.5)
xlabel('Time');
ylabel('Amplitude');
legend('Reference','Target','Warped');
title('Proposed Warping Method');
ylim([0 0.1]);xlim([50 300]);

subplot(312);
plot(ref,'-','LineWidth',1.5);hold on;plot(target,':r','LineWidth',0.5);
[Warping, target_s] = cow(ref',target',10,5,[0 1 1 0 0]);
[Warping, target_s2] = cow(ref',target',20,10,[0 1 1 0 0]);
plot(target_s,':g','LineWidth',1.5)
xlabel('Time');
ylabel('Amplitude');
title('Correlation Optimized Warping (Segments = 10, slack = 5)');
ylim([0 0.1]);xlim([50 300]);
subplot(313);
plot(ref,'-','LineWidth',1.5);hold on;plot(target,':r','LineWidth',0.5);
plot(target_s2,':g','LineWidth',1.5)
xlabel('Time');
ylabel('Amplitude');
title('Correlation Optimized Warping (Segments = 20, slack = 10)');
ylim([0 0.1]);xlim([50 300]);