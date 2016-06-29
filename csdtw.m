% Constrained Selective Dynamic Time Warping
% function dtw = CsDTW(ref,target,segs_start,segs_end,pflag,use_correlation,stretch_first)
% Constrained Selective Dynamic Time Warping
% Inputs:
%   ref - reference tarjectory, column vector
%   target - target trajectory, column vector
%   segs_start - starting indices of the relevant segments in the reference
%        trajectory
%   segs_end   - terminating indices of the relevant segments in the
%        reference trajectory
%   pflag - plotting flag - set to 1 enable to show diagnostic plots
%   use_correlation - instead of solving subsequence DTW, we can apply
%         correlation matching
%
%   stretch_first -
%      if reference / target vectors do not  differ by a large amount in length
%        then enable this flag to stretch the target trajectory through interpolation
% Outputs:
%   The output of the CsDTW is stored as a structure file:
%
%   dtw - dynamic time warping results data structure with several fields:
%       dtw.Dist = Dist;  total accumulated distance of the optimal path
%       dtw.D = D;        cumulative distance matrix
%       dtw.w = w;        symmetric warping path
%       dtw.t_warped = t_warped; synchronized warped target trajectory
%       dtw.k = k;        
%       dtw.segs_start = segs_start;  starting point of relevant segment 
%       dtw.segs_end = segs_end;      termination point of relevant segment
%       dtw.new_segs_start      identified starting point in target
%       dtw.new_segs_end        identified terminating point in target

function dtw = mywarping(ref,target,segs_start,segs_end,pflag,use_correlation,stretch_first)
if nargin < 7;stretch_first = 0;end
if nargin < 6;use_correlation = 0;end
if nargin < 5;pflag = 1;end
l = 10; % search domain, within range of 10 shifts
contype = 2;
% define the relevant segments based on ref
% get the dimensions of N and M, check whether it's column or row vector
[M,i] = max(size(ref));
if i ~= 1;ref = ref';end

[N,i] = max(size(target));
if i ~= 1;target = target';end

% by default, we do not stretch the dataset to the same length
if (N ~= M) && (stretch_first);
    target_s = interp1(linspace(1,N,N),target,linspace(1,N,M)');
else
    target_s = target;
end

new_segs_start = [];
new_segs_end = [];

% let's identify the new segment positions based either on correlation or
% subsequence DTW identification.
if use_correlation == 1;
    for i = numel(segs_start):-1:1;
        segment = segs_start(i):segs_end(i);
        [rho, seg_shift] = maxcorr(ref(segment),target_s(segment),l);
        % apply the shift to the target trajectory
        shifted_segment = segment + seg_shift;                
        % Now we know the coordinates of the new segment.
        new_segs_start = [shifted_segment(1) new_segs_start];
        new_segs_end = [shifted_segment(end) new_segs_end];
    end
else
    for i = numel(segs_start):-1:1;
        % apply subsequence DTW matching to identify the target coordinates
        results = subsequenceDTW(target,ref(segs_start(i):segs_end(i)),0);
        l = segs_end(i)-segs_start(i);
        startpoint = round((results.a + results.b - l)/2);
        endpoint = startpoint + l;
        new_segs_start = [startpoint new_segs_start];
        new_segs_end = [endpoint new_segs_end];
    end
end

% Generate bandconstraint + penalty matrix (BC) for warping
BC = zeros(M,N);
% generate diagonal lines to force 1-1 mapping during VIP sensitive regions
for i = 1:numel(segs_start)
    a = [segs_start(i),new_segs_start(i)];
    b = [segs_end(i),new_segs_end(i)];
    start=a(1):b(1);
    finish=a(2):b(2);
    path = [start' finish'];
    
    for j = 1:size(path,1);
        BC(path(j,1),path(j,2)) = 1;
    end
end

% generate the rest of the allowed bands.
ref_node_locs = sort([segs_start segs_end 1 M]);
target_node_locs = sort([new_segs_start new_segs_end 1 M]);

for i = 1:numel(ref_node_locs)-1;
    if any(segs_start == ref_node_locs(i))
        continue;
    end
    rows = ref_node_locs(i):ref_node_locs(i+1);
    cols = target_node_locs(i):target_node_locs(i+1);
    BC(rows,cols) = 1;
end

% Apply bandconstrained DTW warping
[Dist,D,d,k,w,t_warped] = LCdtw(ref,target,BC,contype,pflag);

% save all the useful information here
dtw.Dist = Dist;
dtw.D = D;
dtw.w = w;
dtw.t_warped = t_warped;
dtw.k = k;
dtw.segs_start = segs_start;
dtw.segs_end = segs_end;
dtw.new_segs_start = new_segs_start;
dtw.new_segs_end = new_segs_end;

end

% find the shift that maximizes the correlation between segments
% positive shift means shift the target trajectory to the left to match the ref
% negative shift means shift to the right.
function [rho, shift] = maxcorr(r,t,L)
if L > numel(r);error('The search slack parameter L exceeds the segment length.');end;
% forward
k = 1;
for i = -L:L;
    % backward lag
    if i < 0
        rho(k) = corr(t(1:end+i),r(-i+1:end));    
    % forward lag
    else
        rho(k) = corr(r(1:end-i),t(i+1:end));
    end
    k = k + 1;
end

[rho, k_max] = max(rho);
shift = k_max - L -1;
end