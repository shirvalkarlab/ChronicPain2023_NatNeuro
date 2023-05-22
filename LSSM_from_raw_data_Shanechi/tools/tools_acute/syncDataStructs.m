% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Omid Sani, Maryam Shanechi
%   Shanechi Lab, University of Southern California, 2018
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%syncDataStructs Syncs two data structs by shedding samples from each that
%do not have a corresponding sample in the other, while ensuring one
%samples form the first stay in the lead.
%   Inputs:
%     - (1) leaderIn: Leader data struct. Closeness will be determined based
%               on the Fs of this struct. Only samples from the leaderIn
%               that have a corresponding sample from followerIn within
%               1/Fs after then will be kept.
%     - (2) followerIn: Follower data struct. Only samples from followerIn
%               that are within 1/Fs after some kept sample from leaderIn
%               will be kept (if there are more than onee, only the first
%               will be kept)
%     - (3) Ts (optional, default: 1/leaderIn.Fs): sampling period of the 
%               leader. 
%     
%   Outputs:
%     - (1) leaderOut: synced output of leader
%     - (2) followerOut: synced output of follower
%     - (3) leaderKeepInd: Kept indices of leaderIn
%     - (4) followerKeepInd: Kept indices of followerIn
%   Usage example:
%       % Sync with joints and save
%       [featsSync2, jDataSync, keepIndFeat2, keptIndJ2] = syncDataStructs(feats, jData);

function [leaderOut, followerOut, leaderKeepInd, followerKeepInd] = syncDataStructs(leaderIn, followerIn, Ts)

if nargin < 3, Ts = 1/leaderIn.Fs; end

delta1 = [0 1]*Ts;
followerKeepInd = [];
for i = 1:length(leaderIn.time)
    tPer = leaderIn.time(i) + delta1;
    ind = find(followerIn.time >= tPer(1) & followerIn.time < tPer(2), 1, 'first');
    if i == 1 || ~isequal(ind, followerKeepInd(end))
        followerKeepInd = cat(1, followerKeepInd, ind);
    end
end

followerOut = followerIn;
followerOut.data = followerOut.data(followerKeepInd, :);
followerOut.time = followerOut.time(followerKeepInd,:);

delta2 = [-1 0]*Ts;
leaderKeepInd = [];
for i = 1:length(followerOut.time)
    tPer = followerOut.time(i) + delta2;
    ind = find(leaderIn.time > tPer(1) & leaderIn.time <= tPer(2), 1, 'last');
    leaderKeepInd = cat(1, leaderKeepInd, ind);
end

leaderOut = leaderIn;
leaderOut.data = leaderOut.data(leaderKeepInd, :);
leaderOut.time = leaderOut.time(leaderKeepInd, :);

end