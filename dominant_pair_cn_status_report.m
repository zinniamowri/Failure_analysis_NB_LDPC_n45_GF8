%% DOMINANT VN-PAIR AND CN-STATUS REPORT
clear;
clc;

matFile = 'failed_frame_log_16p0dB.mat';

if ~isfile(matFile)
    error('Cannot find MAT file: %s', matFile);
end

S = load(matFile);

requiredFields = {'fail_log','CN_lst','Eb_No_db'};
for k = 1:numel(requiredFields)
    if ~isfield(S, requiredFields{k})
        error('Missing required field: %s', requiredFields{k});
    end
end

% Only these two dominant erroneous VN pairs are reported.
dominantPairs = [26 31; 22 37];

snrText = strrep(sprintf('%.1f', S.Eb_No_db), '.', 'p');
reportFile = sprintf('dominant_pair_cn_status_%sdB.txt', snrText);

fid = fopen(reportFile, 'w');
if fid == -1
    error('Could not open report file: %s', reportFile);
end
cleanupObj = onCleanup(@() fclose(fid));

fprintf(fid, 'DOMINANT VN-PAIR CN-STATUS REPORT\n');
fprintf(fid, 'Eb/N0 = %.1f dB\n', S.Eb_No_db);
fprintf(fid, '============================================================\n\n');

for p = 1:size(dominantPairs,1)

    pair = dominantPairs(p,:);

    % Find CNs connected to either VN in the pair.
    connectedCNs = [];
    for c = 1:numel(S.CN_lst)
        cnEntry = S.CN_lst(c);
        vnList = double(cnEntry{1});
        vnList = vnList(:).';

        if any(ismember(pair, vnList))
            connectedCNs(end+1) = c; 
        end
    end

    connectedCNs = unique(connectedCNs, 'stable');

    % Keep only failures whose final erroneous VN set is exactly this pair.
    matchingRecords = [];
    for r = 1:numel(S.fail_log)
        errVNs = sort(double(S.fail_log(r).err_pos(:).'));
        if isequal(errVNs, pair)
            matchingRecords(end+1) = r; %#ok<SAGROW>
        end
    end

    fprintf(fid, 'Dominant pair: (v_%d, v_%d)\n', pair(1), pair(2));
    fprintf(fid, 'Connected CNs: %s\n', mat2str(connectedCNs));
    fprintf(fid, 'Number of failed frames: %d\n\n', numel(matchingRecords));

    % Count how often each connected CN is satisfied/unsatisfied.
    unsatCount = zeros(size(connectedCNs));
    satCount   = zeros(size(connectedCNs));

    fprintf(fid, 'Frame-by-frame CN status:\n');
    fprintf(fid, '------------------------------------------------------------\n');

    for q = 1:numel(matchingRecords)
        rec = S.fail_log(matchingRecords(q));
        syndrome = double(rec.final_syn(:).');

        fprintf(fid, 'Frame %d: ', rec.frame);

        for j = 1:numel(connectedCNs)
            c = connectedCNs(j);

            if syndrome(c) ~= 0
                statusText = 'UNSATISFIED';
                unsatCount(j) = unsatCount(j) + 1;
            else
                statusText = 'satisfied';
                satCount(j) = satCount(j) + 1;
            end

            fprintf(fid, 'CN %d = %s', c, statusText);
            if j < numel(connectedCNs)
                fprintf(fid, ', ');
            end
        end
        fprintf(fid, '\n');
    end

    fprintf(fid, '\nCN-status summary:\n');
    fprintf(fid, '------------------------------------------------------------\n');

    for j = 1:numel(connectedCNs)
        fprintf(fid, ['CN %d: UNSATISFIED in %d frames, ' ...
                      'satisfied in %d frames\n'], ...
                      connectedCNs(j), unsatCount(j), satCount(j));
    end

    fprintf(fid, '\n============================================================\n\n');
end

fprintf('Report created successfully:\n%s\n', reportFile);
